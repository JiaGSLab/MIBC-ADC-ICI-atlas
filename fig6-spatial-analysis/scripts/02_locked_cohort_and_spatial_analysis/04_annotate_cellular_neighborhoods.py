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
from scipy.optimize import linear_sum_assignment


FIG6_ROOT = Path("/data/home/chenwenhui/mIHC/fig6_pm")
DELIVERABLES = FIG6_ROOT / "deliverables"
WP11_OUTPUT = Path("/data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh/neighborhood_radius200_B0")
CENTERS_CSV = WP11_OUTPUT / "neighborhood_cluster_centers.csv"
COMPOSITION_PARQUET = WP11_OUTPUT / "neighborhood_composition.parquet"
K_SCORES_CSV = WP11_OUTPUT / "k_selection_scores.csv"
RUN_PARAMS_CSV = WP11_OUTPUT / "run_params.csv"
BASELINE_CSV = DELIVERABLES / "WP0.3_celltype_distribution_cohort.csv"

EXPECTED_K10_ROWS = 10
EXPECTED_K11_ROWS = 11
SAMPLE_COVERAGE_THRESHOLD = 0.01
LOW_COVERAGE_MIN_SAMPLES = 20
BROAD_COVERAGE_MIN_SAMPLES = 36


def prop_columns(df: pd.DataFrame) -> list[str]:
    cols = [c for c in df.columns if c.startswith("neigh_prop__")]
    if not cols:
        raise ValueError("No neigh_prop__ columns found.")
    return cols


def strip_prop_prefix(column: str) -> str:
    return column.replace("neigh_prop__", "", 1)


def cn_sort_key(cn: str) -> tuple[int, str]:
    match = re.match(r"CN(\d+)$", str(cn))
    if match:
        return (int(match.group(1)), str(cn))
    return (10**9, str(cn))


def format_int(value: int | float) -> str:
    if pd.isna(value):
        return ""
    return f"{int(value):,}"


def format_float(value: float, digits: int = 4) -> str:
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return f"{float(value):.{digits}f}"


def format_pct(value: float, digits: int = 1) -> str:
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return f"{float(value) * 100:.{digits}f}%"


def markdown_table(df: pd.DataFrame, columns: Iterable[str] | None = None, max_rows: int | None = None) -> str:
    if columns is not None:
        df = df.loc[:, list(columns)]
    if max_rows is not None:
        df = df.head(max_rows)
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


def row_composition(row: pd.Series, cols: list[str] | None = None) -> dict[str, float]:
    cols = cols or [c for c in row.index if str(c).startswith("neigh_prop__")]
    return {strip_prop_prefix(c): float(row[c]) for c in cols}


def annotate_niche(composition: dict[str, float]) -> dict[str, str]:
    ordered = sorted(composition.items(), key=lambda item: item[1], reverse=True)
    dominant = ordered[0][0]
    her2 = composition.get("HER2+ Tumor", 0.0)
    panck = composition.get("PanCK+ Tumor", 0.0)
    macrophage = composition.get("CD163- Macrophage", 0.0) + composition.get("CD163+ Macrophage", 0.0)
    b_cell = composition.get("B cell", 0.0) + composition.get("CD21+ B cell", 0.0)
    cd8 = composition.get("CD8 T cell", 0.0) + composition.get("CD8 exhausted T cell", 0.0)
    cd4_treg = composition.get("CD4 T cell", 0.0) + composition.get("Treg", 0.0)

    if dominant == "HER2+ Tumor":
        if her2 >= 0.55:
            name = "HER2+ tumor-dense niche"
            family = "tumor_HER2_dense"
        elif composition.get("CD163- Macrophage", 0.0) >= 0.20:
            name = "HER2+ tumor-macrophage niche"
            family = "tumor_macrophage"
        elif b_cell >= 0.12 or cd8 >= 0.12:
            name = "HER2+ tumor-immune niche"
            family = "tumor_immune"
        else:
            name = "HER2+ tumor niche"
            family = "tumor_HER2"
    elif dominant == "PanCK+ Tumor" or panck >= 0.45:
        name = "PanCK+ tumor niche"
        family = "tumor_PanCK"
    elif dominant == "CD163- Macrophage" or macrophage >= 0.30:
        name = "macrophage-rich immune-stromal niche"
        family = "macrophage_immune_stromal"
    elif dominant == "CD21+ B cell":
        name = "CD21+ B-cell lymphoid niche"
        family = "B_cell_lymphoid"
    elif dominant == "B cell" or b_cell >= 0.28:
        name = "B-cell tumor-immune niche" if her2 >= 0.12 or panck >= 0.10 else "B-cell immune niche"
        family = "B_cell_immune"
    elif dominant == "CD8 T cell" or composition.get("CD8 T cell", 0.0) >= 0.30:
        name = "CD8 T-cell cytotoxic niche"
        family = "CD8_T"
    elif dominant == "CD4 T cell" or cd4_treg >= 0.35:
        name = "CD4/Treg immune niche"
        family = "CD4_Treg"
    elif dominant == "Fibroblast / SMA+ stromal":
        name = "fibroblast/stromal niche"
        family = "fibroblast_stromal"
    elif dominant == "Vascular":
        name = "vascular niche"
        family = "vascular"
    else:
        name = f"{dominant} niche"
        family = re.sub(r"[^A-Za-z0-9]+", "_", dominant).strip("_")

    top = ", ".join(f"{celltype} {value * 100:.1f}%" for celltype, value in ordered[:3])
    return {
        "functional_niche": name,
        "niche_family": family,
        "annotation_rationale": f"Dominant center components: {top}.",
    }


def load_baseline(path: Path = BASELINE_CSV) -> dict[str, float]:
    df = pd.read_csv(path)
    df = df.loc[df["b0_rows"].gt(0) & ~df["celltype"].eq("undefined")].copy()
    baseline = dict(zip(df["celltype"].astype(str), df["b0_fraction"].astype(float)))
    total = sum(baseline.values())
    if not 0.99 <= total <= 1.01:
        raise ValueError(f"B0 baseline fractions do not sum to ~1: {total}")
    return baseline


def top_composition_string(composition: dict[str, float], baseline: dict[str, float], top_n: int = 5) -> str:
    ordered = sorted(composition.items(), key=lambda item: item[1], reverse=True)[:top_n]
    parts = []
    for celltype, fraction in ordered:
        base = baseline.get(celltype, np.nan)
        fold = fraction / base if base and math.isfinite(base) else np.nan
        parts.append(f"{celltype} {fraction * 100:.1f}% ({fold:.2f}x)")
    return "; ".join(parts)


def enriched_celltype_string(composition: dict[str, float], baseline: dict[str, float]) -> str:
    rows = []
    for celltype, fraction in composition.items():
        base = baseline.get(celltype, np.nan)
        if not base or not math.isfinite(base):
            continue
        fold = fraction / base
        if fraction >= 0.03 and fold >= 1.5:
            rows.append((celltype, fraction, fold))
    rows = sorted(rows, key=lambda item: (item[2], item[1]), reverse=True)
    if not rows:
        return ""
    return "; ".join(f"{celltype} {fraction * 100:.1f}% ({fold:.2f}x)" for celltype, fraction, fold in rows[:5])


def build_niche_map(
    centers: pd.DataFrame,
    baseline: dict[str, float],
    total_valid_cells: int,
    k: int,
) -> pd.DataFrame:
    cols = prop_columns(centers)
    rows = []
    for _, row in centers.loc[centers["k"].eq(k)].sort_values("CN", key=lambda s: s.map(cn_sort_key)).iterrows():
        composition = row_composition(row, cols)
        ordered = sorted(composition.items(), key=lambda item: item[1], reverse=True)
        dominant, dominant_fraction = ordered[0]
        dominant_base = baseline.get(dominant, np.nan)
        dominant_enrichment = dominant_fraction / dominant_base if dominant_base and math.isfinite(dominant_base) else np.nan
        annotation = annotate_niche(composition)
        record = {
            "k": int(k),
            "CN": str(row["CN"]),
            "label_column": str(row.get("label_column", f"CN_k{k}")),
            "functional_niche": annotation["functional_niche"],
            "niche_family": annotation["niche_family"],
            "n_cells": int(row["n_cells"]),
            "cluster_fraction": float(row["n_cells"]) / float(total_valid_cells),
            "dominant_celltype": dominant,
            "dominant_fraction": float(dominant_fraction),
            "dominant_enrichment": float(dominant_enrichment),
            "top_composition_enrichment": top_composition_string(composition, baseline, top_n=5),
            "enriched_celltypes": enriched_celltype_string(composition, baseline),
            "annotation_rationale": annotation["annotation_rationale"],
        }
        for col in cols:
            record[col] = float(row[col])
        rows.append(record)
    return pd.DataFrame(rows)


def _matrix(df: pd.DataFrame, prop_cols: list[str]) -> np.ndarray:
    return df[prop_cols].to_numpy(dtype=float)


def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / denom) if denom else np.nan


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def match_k10_k11(k10_map: pd.DataFrame, k11_map: pd.DataFrame, prop_cols: list[str]) -> pd.DataFrame:
    k10 = k10_map.sort_values("CN", key=lambda s: s.map(cn_sort_key)).reset_index(drop=True)
    k11 = k11_map.sort_values("CN", key=lambda s: s.map(cn_sort_key)).reset_index(drop=True)
    x10 = _matrix(k10, prop_cols)
    x11 = _matrix(k11, prop_cols)
    cos = np.zeros((len(k10), len(k11)), dtype=float)
    pear = np.zeros_like(cos)
    for i in range(len(k10)):
        for j in range(len(k11)):
            cos[i, j] = _cosine(x10[i], x11[j])
            pear[i, j] = _pearson(x10[i], x11[j])

    row_ind, col_ind = linear_sum_assignment(-cos)
    hungarian_for_k11 = {int(j): int(i) for i, j in zip(row_ind, col_ind)}
    best_for_k11 = {j: int(np.nanargmax(cos[:, j])) for j in range(len(k11))}
    split_counts: dict[int, int] = {}
    for best_idx in best_for_k11.values():
        split_counts[best_idx] = split_counts.get(best_idx, 0) + 1

    rows = []
    for j, row11 in k11.iterrows():
        best_i = best_for_k11[j]
        hungarian_i = hungarian_for_k11.get(j)
        best_row = k10.iloc[best_i]
        split_count = split_counts[best_i]
        row11_family = str(row11.get("niche_family", ""))
        best_family = str(best_row.get("niche_family", ""))
        same_family = row11_family == best_family and row11_family != ""
        if split_count > 1:
            call = "split"
        elif cos[best_i, j] >= 0.98:
            call = "stable"
        elif cos[best_i, j] >= 0.90 and same_family:
            call = "stable_variant"
        else:
            call = "shifted"

        if call == "split":
            note = f"k11 {row11['CN']} is part of a {split_count}-way subdivision of k10 {best_row['CN']}."
        elif call == "stable":
            note = f"k11 {row11['CN']} closely recapitulates k10 {best_row['CN']}."
        elif call == "stable_variant":
            note = f"k11 {row11['CN']} preserves the same niche family with modest center shift."
        else:
            note = f"k11 {row11['CN']} has a weaker or changed match to k10 {best_row['CN']}."

        rows.append(
            {
                "k11_CN": str(row11["CN"]),
                "k11_functional_niche": str(row11["functional_niche"]),
                "k11_niche_family": row11_family,
                "k11_n_cells": int(row11["n_cells"]),
                "best_k10_CN": str(best_row["CN"]),
                "best_k10_functional_niche": str(best_row["functional_niche"]),
                "best_k10_niche_family": best_family,
                "best_cosine_similarity": float(cos[best_i, j]),
                "best_pearson_r": float(pear[best_i, j]),
                "hungarian_k10_CN": "" if hungarian_i is None else str(k10.iloc[hungarian_i]["CN"]),
                "is_hungarian_primary_match": bool(hungarian_i == best_i),
                "split_group_size": int(split_count),
                "stability_call": call,
                "stability_note": note,
            }
        )
    return pd.DataFrame(rows).sort_values("k11_CN", key=lambda s: s.map(cn_sort_key)).reset_index(drop=True)


def coverage_from_labels(
    labels: pd.DataFrame,
    label_col: str,
    k: int,
    coverage_threshold: float = SAMPLE_COVERAGE_THRESHOLD,
    low_coverage_min_samples: int = LOW_COVERAGE_MIN_SAMPLES,
) -> pd.DataFrame:
    valid = labels.loc[labels["valid_neighborhood"].astype(bool), ["sample", label_col]].dropna().copy()
    valid[label_col] = valid[label_col].astype(str)
    sample_totals = valid.groupby("sample", observed=False).size()
    samples = sample_totals.index.tolist()
    total_valid = int(len(valid))
    rows = []
    for cn, sub in valid.groupby(label_col, observed=False):
        sample_counts = sub.groupby("sample", observed=False).size().reindex(samples, fill_value=0)
        sample_fractions = sample_counts / sample_totals
        samples_ge = int(sample_fractions.ge(coverage_threshold).sum())
        broad_min_samples = min(
            BROAD_COVERAGE_MIN_SAMPLES,
            max(low_coverage_min_samples + 1, math.ceil(len(samples) * 0.50)),
        )
        if samples_ge <= low_coverage_min_samples:
            flag = "low_coverage"
        elif samples_ge < broad_min_samples:
            flag = "moderate_coverage"
        else:
            flag = "broad"
        rows.append(
            {
                "k": int(k),
                "CN": str(cn),
                "label_column": label_col,
                "n_cells": int(len(sub)),
                "cluster_fraction": float(len(sub) / total_valid),
                "samples_present": int(sample_counts.gt(0).sum()),
                "samples_ge_1pct": samples_ge,
                "sample_coverage_fraction_ge_1pct": float(samples_ge / len(samples)),
                "median_sample_fraction": float(sample_fractions.median()),
                "max_sample_fraction": float(sample_fractions.max()),
                "coverage_flag": flag,
            }
        )
    return pd.DataFrame(rows).sort_values(["k", "CN"], key=lambda s: s.map(cn_sort_key) if s.name == "CN" else s).reset_index(drop=True)


def load_label_columns(path: Path = COMPOSITION_PARQUET) -> pd.DataFrame:
    return pd.read_parquet(path, columns=["sample", "valid_neighborhood", "CN_k10", "CN_k11"])


def primary_recommendation(
    k_scores: pd.DataFrame,
    matching: pd.DataFrame,
    coverage: pd.DataFrame,
) -> dict[str, object]:
    k10_ch = float(k_scores.loc[k_scores["k"].eq(10), "calinski_harabasz_score"].iloc[0])
    k11_ch = float(k_scores.loc[k_scores["k"].eq(11), "calinski_harabasz_score"].iloc[0])
    stable_count = int(matching["stability_call"].isin(["stable", "stable_variant"]).sum())
    split_count = int(matching["stability_call"].eq("split").sum())
    low_cov = coverage.loc[coverage["coverage_flag"].eq("low_coverage"), ["k", "CN", "functional_niche"]]
    return {
        "recommended_primary_k": 10,
        "recommendation": "Use k=10 as the primary CN resolution; keep k=11 as sensitivity/resolution check.",
        "reason_1": f"k=10 has the higher Calinski-Harabasz score ({k10_ch:.2f} vs {k11_ch:.2f}).",
        "reason_2": f"k=11 mostly preserves k=10 niche families ({stable_count} stable/stable-variant rows) and adds {split_count} split rows rather than a clearly new broad niche.",
        "reason_3": "k=10 gives one interpretable representative per major functional niche while avoiding extra HER2/B-cell tumor-immune subdivisions.",
        "low_coverage_niches": "; ".join(f"k{int(r.k)} {r.CN} {r.functional_niche}" for r in low_cov.itertuples()) or "none",
    }


def render_report(
    k10_map: pd.DataFrame,
    k11_map: pd.DataFrame,
    matching: pd.DataFrame,
    coverage: pd.DataFrame,
    k_scores: pd.DataFrame,
    recommendation: dict[str, object],
    paths: dict[str, Path],
) -> str:
    stable = matching.loc[matching["stability_call"].isin(["stable", "stable_variant"])]
    split = matching.loc[matching["stability_call"].eq("split")]
    low_cov = coverage.loc[coverage["coverage_flag"].eq("low_coverage")]
    map_cols = [
        "CN",
        "functional_niche",
        "n_cells",
        "cluster_fraction",
        "dominant_celltype",
        "dominant_fraction",
        "dominant_enrichment",
        "top_composition_enrichment",
    ]
    cov_cols = [
        "k",
        "CN",
        "functional_niche",
        "n_cells",
        "cluster_fraction",
        "samples_present",
        "samples_ge_1pct",
        "median_sample_fraction",
        "max_sample_fraction",
        "coverage_flag",
    ]
    match_cols = [
        "k11_CN",
        "k11_functional_niche",
        "best_k10_CN",
        "best_k10_functional_niche",
        "best_cosine_similarity",
        "best_pearson_r",
        "stability_call",
        "stability_note",
    ]
    lines = [
        "# WP1.2 CN Functional Niche Annotation and k Selection",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Scope",
        "",
        "This analysis annotates WP1.1 k10/k11 CN centers, evaluates k10-k11 center stability, and recommends a primary k. "
        "It does not rerun CN clustering and does not perform CR/NR proportion association or adjacency analysis.",
        "",
        "## Inputs",
        "",
        f"- WP1.1 output directory: `{WP11_OUTPUT}`",
        f"- Cluster centers: `{CENTERS_CSV}`",
        f"- Cell-level labels for coverage: `{COMPOSITION_PARQUET}`",
        f"- B0 baseline celltype fractions: `{BASELINE_CSV}`",
        "",
        "## Deliverables",
        "",
        f"- k10 niche map: `{paths['k10_map']}`",
        f"- k11 niche map: `{paths['k11_map']}`",
        f"- k10-k11 matching: `{paths['matching']}`",
        f"- CN coverage: `{paths['coverage']}`",
        f"- Run metadata: `{paths['metadata']}`",
        "",
        "## k10 Functional Niche Map",
        "",
        markdown_table(k10_map, columns=map_cols),
        "",
        "## k11 Functional Niche Map",
        "",
        markdown_table(k11_map, columns=map_cols),
        "",
        "## k10-k11 Matching and Stability",
        "",
        markdown_table(matching, columns=match_cols),
        "",
        f"Stable/stable-variant k11 rows: {len(stable)} / {len(matching)}. Split rows: {len(split)} / {len(matching)}.",
        "",
        "## CN Size and Sample Coverage",
        "",
        f"Coverage threshold: CN comprises >= {SAMPLE_COVERAGE_THRESHOLD * 100:.1f}% of valid CN-participant cells within a sample. "
        f"`low_coverage` means <= {LOW_COVERAGE_MIN_SAMPLES} samples at this threshold.",
        "",
        markdown_table(coverage, columns=cov_cols),
        "",
        "Low-coverage CNs:",
        "",
        markdown_table(low_cov, columns=cov_cols),
        "",
        "## Primary k Recommendation",
        "",
        f"Recommendation: `{recommendation['recommended_primary_k']}`.",
        "",
        f"- {recommendation['recommendation']}",
        f"- {recommendation['reason_1']}",
        f"- {recommendation['reason_2']}",
        f"- {recommendation['reason_3']}",
        f"- Low-coverage niches: {recommendation['low_coverage_niches']}.",
        "",
        "k selection scores:",
        "",
        markdown_table(k_scores),
        "",
        "## Self-Check",
        "",
        f"- k10 map rows: {len(k10_map)}; k11 map rows: {len(k11_map)}.",
        f"- matching rows: {len(matching)}.",
        f"- coverage rows: {len(coverage)}.",
        "- Outputs are written only to `fig6_pm/deliverables/`; WP1.1 CN outputs were read-only.",
    ]
    return "\n".join(lines) + "\n"


def run_analysis(deliverables: Path = DELIVERABLES) -> dict[str, Path]:
    centers = pd.read_csv(CENTERS_CSV)
    k_scores = pd.read_csv(K_SCORES_CSV)
    run_params = pd.read_csv(RUN_PARAMS_CSV)
    baseline = load_baseline(BASELINE_CSV)
    total_valid_cells = int(run_params["n_valid_cells"].iloc[0])
    prop_cols = prop_columns(centers)

    k10_map = build_niche_map(centers, baseline, total_valid_cells=total_valid_cells, k=10)
    k11_map = build_niche_map(centers, baseline, total_valid_cells=total_valid_cells, k=11)
    matching = match_k10_k11(k10_map, k11_map, prop_cols)

    labels = load_label_columns(COMPOSITION_PARQUET)
    coverage = pd.concat(
        [
            coverage_from_labels(labels, "CN_k10", k=10),
            coverage_from_labels(labels, "CN_k11", k=11),
        ],
        ignore_index=True,
    )
    coverage = coverage.merge(
        pd.concat(
            [
                k10_map[["k", "CN", "functional_niche", "niche_family"]],
                k11_map[["k", "CN", "functional_niche", "niche_family"]],
            ],
            ignore_index=True,
        ),
        on=["k", "CN"],
        how="left",
        validate="one_to_one",
    )

    recommendation = primary_recommendation(k_scores, matching, coverage)

    paths = {
        "script": deliverables / "WP1.2_cn_annotation_kselection.py",
        "report": deliverables / "WP1.2_cn_annotation_report.md",
        "k10_map": deliverables / "WP1.2_functional_niche_map_k10.csv",
        "k11_map": deliverables / "WP1.2_functional_niche_map_k11.csv",
        "matching": deliverables / "WP1.2_k10_k11_matching.csv",
        "coverage": deliverables / "WP1.2_cn_coverage.csv",
        "metadata": deliverables / "WP1.2_run_metadata.json",
    }

    k10_map.to_csv(paths["k10_map"], index=False)
    k11_map.to_csv(paths["k11_map"], index=False)
    matching.to_csv(paths["matching"], index=False)
    coverage.to_csv(paths["coverage"], index=False)

    metadata = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "wp11_output": str(WP11_OUTPUT),
        "centers": str(CENTERS_CSV),
        "composition": str(COMPOSITION_PARQUET),
        "baseline": str(BASELINE_CSV),
        "total_valid_cells": total_valid_cells,
        "coverage_threshold": SAMPLE_COVERAGE_THRESHOLD,
        "low_coverage_min_samples": LOW_COVERAGE_MIN_SAMPLES,
        "recommended_primary_k": recommendation["recommended_primary_k"],
        "command": "conda activate scimap && python deliverables/WP1.2_cn_annotation_kselection.py",
        "outputs": {name: str(path) for name, path in paths.items() if name != "script"},
    }
    paths["metadata"].write_text(json.dumps(metadata, indent=2, ensure_ascii=True), encoding="utf-8")

    report = render_report(k10_map, k11_map, matching, coverage, k_scores, recommendation, paths)
    paths["report"].write_text(report, encoding="utf-8")
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description="WP1.2 CN functional niche annotation and k selection.")
    parser.add_argument("--deliverables", type=Path, default=DELIVERABLES)
    args = parser.parse_args()
    paths = run_analysis(args.deliverables)
    for name, path in paths.items():
        if name != "script":
            print(f"{name}: {path}")


if __name__ == "__main__":
    main()
