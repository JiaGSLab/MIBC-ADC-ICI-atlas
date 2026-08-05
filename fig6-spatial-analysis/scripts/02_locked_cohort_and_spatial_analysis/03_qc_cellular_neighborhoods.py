#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import stats


FIG6_ROOT = Path("/data/home/chenwenhui/mIHC/fig6_pm")
OUTPUT_DIR = Path("/data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh/neighborhood_radius200_B0")
SAMPLE_VALIDATION = FIG6_ROOT / "deliverables/WP0.3_sample_validation.csv"
DELIVERABLES = FIG6_ROOT / "deliverables"

EXPECTED_ROWS = 23_087_820
EXPECTED_SAMPLES = 72
EXPECTED_K = [10, 11]
LOW_CELL_THRESHOLD = 30_000


def format_int(value: int | float) -> str:
    if pd.isna(value):
        return ""
    return f"{int(value):,}"


def format_float(value: float, digits: int = 4) -> str:
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return f"{float(value):.{digits}f}"


def format_pct(value: float, digits: int = 2) -> str:
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return f"{float(value) * 100:.{digits}f}%"


def markdown_table(df: pd.DataFrame, columns: Iterable[str] | None = None) -> str:
    if columns is not None:
        df = df.loc[:, list(columns)]
    if df.empty:
        return "_None._"

    render = df.copy()
    for column in render.columns:
        if pd.api.types.is_float_dtype(render[column]):
            render[column] = render[column].map(lambda x: format_float(x, 4))
        elif pd.api.types.is_integer_dtype(render[column]):
            render[column] = render[column].map(lambda x: format_int(x))
        else:
            render[column] = render[column].astype(str).replace({"nan": ""})

    headers = list(render.columns)
    rows = render.astype(str).values.tolist()
    out = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    out.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(out)


def phenotype_prop_columns(df: pd.DataFrame) -> list[str]:
    columns = [column for column in df.columns if column.startswith("neigh_prop__")]
    if not columns:
        raise ValueError("No neighborhood proportion columns found.")
    return columns


def strip_prop_prefix(column: str) -> str:
    return column.replace("neigh_prop__", "", 1)


def cn_sort_key(cn: str) -> tuple[int, str]:
    match = re.match(r"CN(\d+)$", str(cn))
    if match:
        return (int(match.group(1)), str(cn))
    return (10**9, str(cn))


def top_dominant_celltypes(centers: pd.DataFrame, top_n: int = 3) -> pd.DataFrame:
    prop_cols = phenotype_prop_columns(centers)
    rows = []
    for _, row in centers.iterrows():
        values = row[prop_cols].astype(float).sort_values(ascending=False)
        top = values.head(top_n)
        rows.append(
            {
                "k": int(row["k"]),
                "CN": str(row["CN"]),
                "label_column": str(row.get("label_column", "")),
                "n_cells": int(row["n_cells"]),
                "dominant_celltype": strip_prop_prefix(top.index[0]),
                "dominant_fraction": float(top.iloc[0]),
                "top_celltypes": "; ".join(
                    f"{strip_prop_prefix(column)} {float(value) * 100:.1f}%"
                    for column, value in top.items()
                ),
            }
        )
    out = pd.DataFrame(rows)
    out = out.sort_values(["k", "CN"], key=lambda series: series.map(cn_sort_key) if series.name == "CN" else series)
    return out.reset_index(drop=True)


def sample_qc_table(
    sample_validation: pd.DataFrame,
    count_summary: pd.DataFrame,
    low_cell_threshold: int = LOW_CELL_THRESHOLD,
) -> pd.DataFrame:
    left = sample_validation.copy()
    right = count_summary[
        [
            "sample",
            "n_cells",
            "neighbor_count_min",
            "neighbor_count_q01",
            "neighbor_count_median",
            "neighbor_count_mean",
            "n_valid_cells",
            "valid_fraction",
        ]
    ].copy()
    merged = left.merge(right, on="sample", how="outer", validate="one_to_one", indicator=True)
    if not merged["_merge"].eq("both").all():
        missing = merged.loc[~merged["_merge"].eq("both"), ["sample", "_merge"]]
        raise ValueError(f"Sample mismatch between WP0.3 validation and neighborhood summary:\n{missing}")
    merged = merged.drop(columns=["_merge"])

    mismatch = merged["b0_defined_rows"].astype(int) != merged["n_cells"].astype(int)
    if mismatch.any():
        bad = merged.loc[mismatch, ["sample", "b0_defined_rows", "n_cells"]]
        raise ValueError(f"B0 defined counts do not match neighborhood rows:\n{bad}")

    merged["cn_participant_cells"] = merged["n_valid_cells"].astype(int)
    merged["cn_filtered_cells"] = merged["b0_defined_rows"].astype(int) - merged["cn_participant_cells"]
    merged["cn_participant_fraction"] = merged["cn_participant_cells"] / merged["b0_defined_rows"]
    merged["low_cell_lt30000"] = merged["b0_defined_rows"].astype(int) < int(low_cell_threshold)

    columns = [
        "sample",
        "Time_short",
        "Outcome_short",
        "Group",
        "union_rows",
        "b0_defined_rows",
        "defined_fraction",
        "cn_participant_cells",
        "cn_filtered_cells",
        "cn_participant_fraction",
        "low_cell_lt30000",
        "neighbor_count_min",
        "neighbor_count_q01",
        "neighbor_count_median",
        "neighbor_count_mean",
    ]
    return merged.loc[:, columns].sort_values("sample").reset_index(drop=True)


def defined_fraction_group_tables(qc: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    group_stats = (
        qc.groupby("Outcome_short", observed=False)["defined_fraction"]
        .agg(n="count", mean="mean", median="median", sd="std", min="min", max="max")
        .reset_index()
        .sort_values("Outcome_short")
    )

    cr = qc.loc[qc["Outcome_short"].eq("CR"), "defined_fraction"].dropna().astype(float).to_numpy()
    nr = qc.loc[qc["Outcome_short"].eq("NR"), "defined_fraction"].dropna().astype(float).to_numpy()
    if len(cr) < 2 or len(nr) < 2:
        comparison = pd.DataFrame(
            [
                {
                    "comparison": "NR-CR",
                    "n_CR": len(cr),
                    "n_NR": len(nr),
                    "mean_diff_NR_minus_CR": np.nan,
                    "mean_diff_95ci_low": np.nan,
                    "mean_diff_95ci_high": np.nan,
                    "cohens_d": np.nan,
                    "welch_t_p": np.nan,
                    "welch_t_fdr_q": np.nan,
                    "mannwhitney_p": np.nan,
                }
            ]
        )
        return group_stats, comparison

    diff = float(np.mean(nr) - np.mean(cr))
    var_cr = float(np.var(cr, ddof=1))
    var_nr = float(np.var(nr, ddof=1))
    se = math.sqrt(var_cr / len(cr) + var_nr / len(nr))
    df_num = (var_cr / len(cr) + var_nr / len(nr)) ** 2
    df_den = (var_cr / len(cr)) ** 2 / (len(cr) - 1) + (var_nr / len(nr)) ** 2 / (len(nr) - 1)
    dof = df_num / df_den if df_den else np.nan
    tcrit = stats.t.ppf(0.975, dof) if math.isfinite(dof) else np.nan
    ci_low = diff - tcrit * se if math.isfinite(tcrit) else np.nan
    ci_high = diff + tcrit * se if math.isfinite(tcrit) else np.nan

    pooled = math.sqrt(((len(cr) - 1) * var_cr + (len(nr) - 1) * var_nr) / (len(cr) + len(nr) - 2))
    cohens_d = diff / pooled if pooled else np.nan
    welch = stats.ttest_ind(nr, cr, equal_var=False, nan_policy="omit")
    mw = stats.mannwhitneyu(nr, cr, alternative="two-sided")

    comparison = pd.DataFrame(
        [
            {
                "comparison": "NR-CR",
                "n_CR": len(cr),
                "n_NR": len(nr),
                "mean_CR": float(np.mean(cr)),
                "mean_NR": float(np.mean(nr)),
                "mean_diff_NR_minus_CR": diff,
                "mean_diff_95ci_low": float(ci_low),
                "mean_diff_95ci_high": float(ci_high),
                "cohens_d": float(cohens_d),
                "welch_t_p": float(welch.pvalue),
                "welch_t_fdr_q": float(welch.pvalue),
                "mannwhitney_p": float(mw.pvalue),
            }
        ]
    )
    return group_stats, comparison


def dominant_summary_sentence(top_df: pd.DataFrame, k: int) -> str:
    sub = top_df.loc[top_df["k"].eq(k)].copy()
    parts = []
    for celltype, rows in sub.groupby("dominant_celltype", sort=True):
        cns = sorted(rows["CN"].tolist(), key=cn_sort_key)
        parts.append(f"{celltype} ({', '.join(cns)})")
    return f"k={k}: " + "; ".join(parts)


def read_inputs(output_dir: Path) -> dict[str, pd.DataFrame]:
    return {
        "run_params": pd.read_csv(output_dir / "run_params.csv"),
        "k_scores": pd.read_csv(output_dir / "k_selection_scores.csv"),
        "centers": pd.read_csv(output_dir / "neighborhood_cluster_centers.csv"),
        "count_summary": pd.read_csv(output_dir / "neighborhood_count_summary.csv"),
        "manifest": pd.read_csv(output_dir / "sample_neighborhood_manifest.csv"),
        "sample_validation": pd.read_csv(SAMPLE_VALIDATION),
    }


def render_report(
    output_dir: Path,
    deliverables: Path,
    low_cell_threshold: int = LOW_CELL_THRESHOLD,
) -> dict[str, Path]:
    data = read_inputs(output_dir)
    parquet_path = output_dir / "neighborhood_composition.parquet"
    parquet_file = pq.ParquetFile(parquet_path)
    total_rows = int(parquet_file.metadata.num_rows)

    run_params = data["run_params"]
    k_scores = data["k_scores"]
    centers = data["centers"]
    count_summary = data["count_summary"]
    manifest = data["manifest"]

    qc = sample_qc_table(data["sample_validation"], count_summary, low_cell_threshold=low_cell_threshold)
    top = top_dominant_celltypes(centers, top_n=3)
    group_stats, group_comparison = defined_fraction_group_tables(qc)

    qc_path = deliverables / "WP1.1_cn_qc_by_sample.csv"
    top_path = deliverables / "WP1.1_cn_top_dominant_celltypes.csv"
    group_stats_path = deliverables / "WP1.1_defined_fraction_group_stats.csv"
    group_comparison_path = deliverables / "WP1.1_defined_fraction_group_comparison.csv"
    metadata_path = deliverables / "WP1.1_run_metadata.json"
    report_path = deliverables / "WP1.1_cn_on_b0_report.md"

    qc.to_csv(qc_path, index=False)
    top.to_csv(top_path, index=False)
    group_stats.to_csv(group_stats_path, index=False)
    group_comparison.to_csv(group_comparison_path, index=False)

    params = run_params.iloc[0].to_dict()
    valid_cells = int(params["n_valid_cells"])
    filtered_cells = total_rows - valid_cells
    low_samples = qc.loc[qc["low_cell_lt30000"]].copy()
    comparison = group_comparison.iloc[0]

    command = (
        "source ~/miniconda3/etc/profile.d/conda.sh && conda activate scimap && "
        "python /data/home/chenwenhui/mIHC/scripts/compute_radius_neighborhood_composition_and_cluster.py "
        "--input_table /data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh/all_samples_cell_phenotypes_B0.parquet "
        "--output_dir /data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh/neighborhood_radius200_B0 "
        "--radius 200 --min_neighbors 10 --k_values 10 11 --primary_k 10 "
        "--sample_col sample --cellid_col cellid --x_col spatial_x --y_col spatial_y "
        "--phenotype_col celltype --workers 20 --random_state 0 --max_cells_for_k_selection 100000 "
        "--prefer_output parquet"
    )

    metadata = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_table": params["input_table"],
        "output_dir": str(output_dir),
        "neighborhood_composition_rows": total_rows,
        "neighborhood_composition_columns": len(parquet_file.schema_arrow.names),
        "valid_cells": valid_cells,
        "filtered_cells": filtered_cells,
        "samples": int(qc["sample"].nunique()),
        "k_values": EXPECTED_K,
        "command": command,
        "deliverables": {
            "report": str(report_path),
            "sample_qc": str(qc_path),
            "top_dominant_celltypes": str(top_path),
            "group_stats": str(group_stats_path),
            "group_comparison": str(group_comparison_path),
        },
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=True), encoding="utf-8")

    systemic_text = (
        "No clear sample-level CR/NR systematic difference was detected"
        if comparison["welch_t_p"] >= 0.05 and comparison["mean_diff_95ci_low"] <= 0 <= comparison["mean_diff_95ci_high"]
        else "A sample-level CR/NR difference was detected"
    )

    low_columns = [
        "sample",
        "Time_short",
        "Outcome_short",
        "b0_defined_rows",
        "defined_fraction",
        "cn_participant_cells",
        "cn_filtered_cells",
        "cn_participant_fraction",
    ]
    top_columns = ["k", "CN", "n_cells", "dominant_celltype", "dominant_fraction", "top_celltypes"]

    lines = [
        "# WP1.1 CN on B0 Report",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Scope",
        "",
        "This delivery rebuilds radius200 neighborhood composition and KMeans CN labels on the raw B0 table only. "
        "It does not perform CN annotation, downstream celltype/CN proportion association, or adjacency analysis.",
        "",
        "## Inputs and Parameters",
        "",
        f"- Input B0 table: `{params['input_table']}`",
        f"- Output directory: `{output_dir}`",
        f"- Radius: {format_float(float(params['radius']), 1)}",
        f"- min_neighbors: {int(params['min_neighbors'])}",
        f"- k_values: {params['k_values']}; primary_k: {int(params['primary_k'])}",
        f"- columns: sample=`{params['sample_col']}`, cellid=`{params['cellid_col']}`, x=`{params['x_col']}`, y=`{params['y_col']}`, phenotype=`{params['phenotype_col']}`",
        f"- random_state: 0; workers: {int(params['workers'])}; max_cells_for_k_selection: 100,000",
        "",
        "Command:",
        "",
        f"```bash\n{command}\n```",
        "",
        "## Output Inventory",
        "",
        f"- Cell-level neighborhood/CN table: `{parquet_path}` ({format_int(total_rows)} rows, {len(parquet_file.schema_arrow.names)} columns)",
        f"- Per-sample neighborhood intermediates: `{output_dir / 'sample_neighborhoods'}` ({int(manifest['status'].eq('ok').sum())} ok samples)",
        f"- Cluster centers: `{output_dir / 'neighborhood_cluster_centers.csv'}`",
        f"- K selection scores: `{output_dir / 'k_selection_scores.csv'}`",
        f"- Run params: `{output_dir / 'run_params.csv'}`",
        f"- QC CSV: `{qc_path}`",
        "",
        "## Effective Cells",
        "",
        f"- Total B0 cells processed: {format_int(total_rows)}",
        f"- CN participant cells (`neighbor_count >= 10`): {format_int(valid_cells)} ({format_pct(valid_cells / total_rows)})",
        f"- Filtered cells (`neighbor_count < 10`): {format_int(filtered_cells)} ({format_pct(filtered_cells / total_rows)})",
        "",
        "## K Selection",
        "",
        markdown_table(k_scores),
        "",
        "## CN Dominant Celltypes Preview",
        "",
        "The table reports the top center proportions only; these are previews for WP1.2 and are not formal CN annotations.",
        "",
        markdown_table(top, columns=top_columns),
        "",
        "Dominant-celltype summary:",
        "",
        f"- {dominant_summary_sentence(top, 10)}",
        f"- {dominant_summary_sentence(top, 11)}",
        "",
        "## Phase 1 QC",
        "",
        "### Defined Fraction by CR/NR",
        "",
        markdown_table(group_stats),
        "",
        markdown_table(group_comparison),
        "",
        f"{systemic_text}: NR-CR mean difference = {format_pct(float(comparison['mean_diff_NR_minus_CR']))} "
        f"(95% CI {format_pct(float(comparison['mean_diff_95ci_low']))} to {format_pct(float(comparison['mean_diff_95ci_high']))}), "
        f"Cohen's d = {format_float(float(comparison['cohens_d']), 3)}, "
        f"Welch p/q = {format_float(float(comparison['welch_t_p']), 4)}/{format_float(float(comparison['welch_t_fdr_q']), 4)}, "
        f"Mann-Whitney p = {format_float(float(comparison['mannwhitney_p']), 4)}.",
        "",
        "### Low-Cell Samples",
        "",
        f"Threshold: B0 defined cells < {format_int(low_cell_threshold)}. These samples are flagged for downstream per-sample weighting/exclusion decisions only; none were removed in WP1.1.",
        "",
        markdown_table(low_samples[low_columns]),
        "",
        "## Self-Check",
        "",
        f"- Full B0 row count matches expectation: {total_rows == EXPECTED_ROWS} ({format_int(total_rows)} rows).",
        f"- Sample count matches expectation: {qc['sample'].nunique() == EXPECTED_SAMPLES} ({qc['sample'].nunique()} samples).",
        f"- All sample neighborhood jobs completed: {manifest['status'].eq('ok').all()} ({manifest['status'].eq('ok').sum()}/72 ok).",
        f"- CN columns available: `CN_k10` and `CN_k11` in `neighborhood_composition.parquet`.",
        "- No annotation, downstream proportion association, or adjacency outputs were generated.",
    ]
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return {
        "report": report_path,
        "qc": qc_path,
        "top": top_path,
        "group_stats": group_stats_path,
        "group_comparison": group_comparison_path,
        "metadata": metadata_path,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build WP1.1 report and QC deliverables.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--deliverables", type=Path, default=DELIVERABLES)
    parser.add_argument("--low-cell-threshold", type=int, default=LOW_CELL_THRESHOLD)
    args = parser.parse_args()

    paths = render_report(args.output_dir, args.deliverables, low_cell_threshold=args.low_cell_threshold)
    for name, path in paths.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
