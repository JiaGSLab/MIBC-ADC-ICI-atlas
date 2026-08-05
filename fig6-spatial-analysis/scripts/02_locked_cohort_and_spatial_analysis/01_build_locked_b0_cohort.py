#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import time
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from scipy.spatial import cKDTree


FIG6_ROOT = Path("/data/home/chenwenhui/mIHC/fig6_pm")
INPUT_ROOT = Path("/data/home/chenwenhui/mIHC/CellInfo/BLCA_PixieCrosstalk_noSOM_fig6_refresh")
OUTPUT_ROOT = Path("/data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh")
INPUT_PARQUET = INPUT_ROOT / "all_samples_cell_phenotypes.parquet"
SAMPLE_INFO = INPUT_ROOT / "CN_pathology_k_selection_fig6_refresh/sample_info_with_analysis_id.csv"
OUTPUT_PARQUET = OUTPUT_ROOT / "all_samples_cell_phenotypes_B0.parquet"
DELIVERABLES = FIG6_ROOT / "deliverables"
SCRATCH = FIG6_ROOT / "scratch_wp0.3"

EXPECTED_ACTIVE_SAMPLES = 72
NEAR_OVERLAP_RADIUS_PX = 30.0
UM_PER_PX = 0.4994729534954469
BATCH_SIZE = 262_144

METADATA_COLUMNS = [
    "Sample_ID",
    "Therapy",
    "Time_point",
    "Time_short",
    "Outcome",
    "Outcome_short",
    "Group",
    "patient_id_for_matching",
    "Sex",
    "Age",
    "Grade",
    "Invasive",
    "Differentiation",
]


def json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        val = float(obj)
        return None if not math.isfinite(val) else val
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)


def filter_defined_table(table: pa.Table) -> pa.Table:
    celltype = table["celltype"]
    mask = pc.and_(pc.is_valid(celltype), pc.not_equal(celltype, pa.scalar("undefined")))
    return table.filter(mask)


def _metadata_value(value, column: str):
    if pd.isna(value):
        return None
    if column == "Age":
        return float(value)
    return str(value)


def append_metadata_columns(table: pa.Table, metadata: pd.DataFrame, columns: list[str]) -> pa.Table:
    if table.num_rows == 0:
        return table
    meta = metadata.drop_duplicates("sample", keep="first").set_index("sample")
    sample_values = table["sample"].to_pylist()
    arrays = []
    names = []
    for column in columns:
        if column not in metadata.columns:
            continue
        mapping = {sample: _metadata_value(value, column) for sample, value in meta[column].items()}
        values = [mapping.get(sample) for sample in sample_values]
        pa_type = pa.float64() if column == "Age" else pa.string()
        arrays.append(pa.array(values, type=pa_type))
        names.append(column)
    if not arrays:
        return table
    return _append_many(table, names, arrays)


def _append_many(table: pa.Table, names: list[str], arrays: list[pa.Array]) -> pa.Table:
    out = table
    for name, array in zip(names, arrays):
        out = out.append_column(name, array)
    return out


def panel_contribution_table(counts: dict[tuple[str, str], int] | Counter) -> pd.DataFrame:
    sample_totals: Counter = Counter()
    for (sample, _panel), count in counts.items():
        sample_totals[sample] += int(count)
    rows = []
    for (sample, panel), count in sorted(counts.items()):
        total = sample_totals[sample]
        rows.append(
            {
                "sample": sample,
                "panel": panel,
                "b0_defined_rows": int(count),
                "panel_fraction_of_b0": float(count / total) if total else np.nan,
            }
        )
    return pd.DataFrame(rows)


def near_overlap_summary_for_sample(
    sample: str,
    coords: np.ndarray,
    panels: np.ndarray,
    radius_px: float = NEAR_OVERLAP_RADIUS_PX,
) -> tuple[dict, pd.DataFrame]:
    n_rows = int(len(coords))
    if n_rows == 0:
        return {
            "sample": sample,
            "defined_rows": 0,
            "rows_with_cross_panel_neighbor_lt_radius": 0,
            "fraction_rows_with_cross_panel_neighbor_lt_radius": np.nan,
            "radius_px": float(radius_px),
            "radius_um": float(radius_px * UM_PER_PX),
        }, pd.DataFrame()

    panels = panels.astype(str)
    unique_panels = sorted(np.unique(panels).tolist())
    row_has_neighbor = np.zeros(n_rows, dtype=bool)
    min_dist = np.full(n_rows, np.nan, dtype=np.float32)
    pair_rows = []

    for panel_a in unique_panels:
        idx_a = np.flatnonzero(panels == panel_a)
        coords_a = coords[idx_a]
        for panel_b in unique_panels:
            if panel_a == panel_b:
                continue
            idx_b = np.flatnonzero(panels == panel_b)
            if len(idx_a) == 0 or len(idx_b) == 0:
                pair_rows.append(
                    {
                        "sample": sample,
                        "panel_a": panel_a,
                        "panel_b": panel_b,
                        "n_a": int(len(idx_a)),
                        "n_b": int(len(idx_b)),
                        "rows_with_neighbor_lt_radius": 0,
                        "fraction_a_with_neighbor_lt_radius": np.nan,
                        "median_nn_distance_px": np.nan,
                    }
                )
                continue
            tree_b = cKDTree(coords[idx_b])
            dist, _idx = tree_b.query(coords_a, k=1, workers=-1)
            hit = dist < radius_px
            row_has_neighbor[idx_a[hit]] = True
            current = min_dist[idx_a]
            current = np.where(np.isnan(current), dist.astype(np.float32), np.minimum(current, dist.astype(np.float32)))
            min_dist[idx_a] = current
            pair_rows.append(
                {
                    "sample": sample,
                    "panel_a": panel_a,
                    "panel_b": panel_b,
                    "n_a": int(len(idx_a)),
                    "n_b": int(len(idx_b)),
                    "rows_with_neighbor_lt_radius": int(hit.sum()),
                    "fraction_a_with_neighbor_lt_radius": float(hit.mean()) if len(hit) else np.nan,
                    "median_nn_distance_px": float(np.median(dist)) if len(dist) else np.nan,
                    "p10_nn_distance_px": float(np.percentile(dist, 10)) if len(dist) else np.nan,
                    "p90_nn_distance_px": float(np.percentile(dist, 90)) if len(dist) else np.nan,
                    "radius_px": float(radius_px),
                    "radius_um": float(radius_px * UM_PER_PX),
                }
            )

    finite_min = min_dist[np.isfinite(min_dist)]
    summary = {
        "sample": sample,
        "defined_rows": n_rows,
        "rows_with_cross_panel_neighbor_lt_radius": int(row_has_neighbor.sum()),
        "fraction_rows_with_cross_panel_neighbor_lt_radius": float(row_has_neighbor.mean()),
        "median_min_cross_panel_nn_distance_px": float(np.median(finite_min)) if len(finite_min) else np.nan,
        "p10_min_cross_panel_nn_distance_px": float(np.percentile(finite_min, 10)) if len(finite_min) else np.nan,
        "p90_min_cross_panel_nn_distance_px": float(np.percentile(finite_min, 90)) if len(finite_min) else np.nan,
        "radius_px": float(radius_px),
        "radius_um": float(radius_px * UM_PER_PX),
    }
    return summary, pd.DataFrame(pair_rows)


def add_counts(counter: Counter, series: pd.Series) -> None:
    counter.update({key: int(value) for key, value in series.value_counts(dropna=False).items()})


def add_group_counts(counter: Counter, df: pd.DataFrame, columns: list[str]) -> None:
    if df.empty:
        return
    grouped = df.groupby(columns, dropna=False).size()
    counter.update({tuple(map(str, key if isinstance(key, tuple) else (key,))): int(value) for key, value in grouped.items()})


def update_batch_summaries(
    batch_df: pd.DataFrame,
    b0_df: pd.DataFrame,
    counters: dict[str, Counter],
) -> None:
    batch_df = batch_df.copy()
    b0_df = b0_df.copy()
    batch_df["celltype"] = batch_df["celltype"].astype(object).where(batch_df["celltype"].notna(), "__NULL__")
    b0_df["celltype"] = b0_df["celltype"].astype(object).where(b0_df["celltype"].notna(), "__NULL__")

    add_counts(counters["union_sample"], batch_df["sample"].astype(str))
    undefined_mask = batch_df["celltype"].astype(str).eq("undefined") | batch_df["celltype"].astype(str).eq("__NULL__")
    add_counts(counters["undefined_sample"], batch_df.loc[undefined_mask, "sample"].astype(str))
    add_counts(counters["b0_sample"], b0_df["sample"].astype(str))
    add_group_counts(counters["union_sample_celltype"], batch_df[["sample", "celltype"]], ["sample", "celltype"])
    add_group_counts(counters["b0_sample_celltype"], b0_df[["sample", "celltype"]], ["sample", "celltype"])
    add_group_counts(counters["union_celltype"], batch_df[["celltype"]], ["celltype"])
    add_group_counts(counters["b0_celltype"], b0_df[["celltype"]], ["celltype"])
    add_group_counts(counters["b0_panel"], b0_df[["sample", "panel"]], ["sample", "panel"])
    if not b0_df.empty:
        coord_bad = b0_df["spatial_x"].isna() | b0_df["spatial_y"].isna()
        add_counts(counters["b0_bad_coord_sample"], b0_df.loc[coord_bad, "sample"].astype(str))


def collect_coords(b0_df: pd.DataFrame, coord_chunks: dict[str, list[np.ndarray]], panel_chunks: dict[str, list[np.ndarray]]) -> None:
    if b0_df.empty:
        return
    for sample, part in b0_df.groupby("sample", sort=False):
        coords = part[["spatial_x", "spatial_y"]].to_numpy(dtype=np.float32, copy=True)
        panels = part["panel"].astype(str).to_numpy(copy=True)
        coord_chunks[str(sample)].append(coords)
        panel_chunks[str(sample)].append(panels)


def counter_to_celltype_df(union_counter: Counter, b0_counter: Counter, level: str) -> pd.DataFrame:
    keys = sorted(set(union_counter) | set(b0_counter))
    rows = []
    if level == "cohort":
        union_total = sum(union_counter.values())
        b0_total = sum(b0_counter.values())
        for key in keys:
            celltype = key[0] if isinstance(key, tuple) else key
            union_n = int(union_counter.get(key, 0))
            b0_n = int(b0_counter.get(key, 0))
            rows.append(
                {
                    "celltype": celltype,
                    "union_rows": union_n,
                    "b0_rows": b0_n,
                    "dropped_rows": union_n - b0_n,
                    "union_fraction": float(union_n / union_total) if union_total else np.nan,
                    "b0_fraction": float(b0_n / b0_total) if b0_total else np.nan,
                }
            )
    else:
        sample_totals_union: Counter = Counter()
        sample_totals_b0: Counter = Counter()
        for sample, _celltype in union_counter:
            sample_totals_union[sample] += union_counter[(sample, _celltype)]
        for sample, _celltype in b0_counter:
            sample_totals_b0[sample] += b0_counter[(sample, _celltype)]
        for sample, celltype in keys:
            union_n = int(union_counter.get((sample, celltype), 0))
            b0_n = int(b0_counter.get((sample, celltype), 0))
            rows.append(
                {
                    "sample": sample,
                    "celltype": celltype,
                    "union_rows": union_n,
                    "b0_rows": b0_n,
                    "dropped_rows": union_n - b0_n,
                    "union_fraction": float(union_n / sample_totals_union[sample]) if sample_totals_union[sample] else np.nan,
                    "b0_fraction": float(b0_n / sample_totals_b0[sample]) if sample_totals_b0[sample] else np.nan,
                }
            )
    return pd.DataFrame(rows)


def sample_validation_table(counters: dict[str, Counter], sample_info: pd.DataFrame) -> pd.DataFrame:
    samples = sorted(set(counters["union_sample"]) | set(counters["b0_sample"]))
    meta = sample_info.drop_duplicates("sample", keep="first").set_index("sample")
    rows = []
    for sample in samples:
        union_n = int(counters["union_sample"].get(sample, 0))
        b0_n = int(counters["b0_sample"].get(sample, 0))
        undefined_n = int(counters["undefined_sample"].get(sample, 0))
        bad_coord = int(counters["b0_bad_coord_sample"].get(sample, 0))
        rows.append(
            {
                "sample": sample,
                "Time_short": meta.at[sample, "Time_short"] if sample in meta.index and "Time_short" in meta.columns else None,
                "Outcome_short": meta.at[sample, "Outcome_short"] if sample in meta.index and "Outcome_short" in meta.columns else None,
                "Group": meta.at[sample, "Group"] if sample in meta.index and "Group" in meta.columns else None,
                "union_rows": union_n,
                "b0_defined_rows": b0_n,
                "dropped_undefined_rows": undefined_n,
                "defined_fraction": float(b0_n / union_n) if union_n else np.nan,
                "undefined_fraction": float(undefined_n / union_n) if union_n else np.nan,
                "b0_bad_coordinate_rows": bad_coord,
            }
        )
    return pd.DataFrame(rows)


def validate_output_parquet(path: Path, expected_rows: int) -> dict:
    pf = pq.ParquetFile(path)
    rows = int(pf.metadata.num_rows)
    undefined_rows = 0
    bad_coord_rows = 0
    samples: set[str] = set()
    panels: set[str] = set()
    for batch in pf.iter_batches(columns=["sample", "panel", "celltype", "spatial_x", "spatial_y"], batch_size=BATCH_SIZE):
        table = pa.Table.from_batches([batch])
        celltype = table["celltype"]
        undefined_rows += int(pc.sum(pc.equal(celltype, pa.scalar("undefined"))).as_py() or 0)
        bad_x = pc.invert(pc.is_valid(table["spatial_x"]))
        bad_y = pc.invert(pc.is_valid(table["spatial_y"]))
        bad_coord_rows += int(pc.sum(pc.or_(bad_x, bad_y)).as_py() or 0)
        samples.update(table["sample"].unique().to_pylist())
        panels.update(table["panel"].unique().to_pylist())
    return {
        "path": str(path),
        "rows": rows,
        "expected_rows": int(expected_rows),
        "row_count_matches_expected": rows == expected_rows,
        "undefined_rows": undefined_rows,
        "bad_coordinate_rows": bad_coord_rows,
        "sample_count": len(samples),
        "panel_count": len(panels),
        "samples": sorted(samples),
        "panels": sorted(panels),
    }


def format_number(value) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"
    if isinstance(value, (float, np.floating)):
        if abs(float(value)) < 1 and float(value) != 0:
            return f"{float(value):.4f}"
        return f"{float(value):,.2f}"
    return str(value)


def markdown_table(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "_No rows._"
    view = df.head(max_rows).copy()
    columns = [str(col) for col in view.columns]
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for _idx, row in view.iterrows():
        lines.append("| " + " | ".join(format_number(row[col]).replace("|", "\\|") for col in view.columns) + " |")
    return "\n".join(lines)


def make_report(
    output_parquet: Path,
    input_rows: int,
    output_verify: dict,
    sample_validation: pd.DataFrame,
    cohort_celltype: pd.DataFrame,
    sample_celltype: pd.DataFrame,
    panel_contrib: pd.DataFrame,
    near_summary: pd.DataFrame,
    run_metadata: dict,
) -> str:
    total_b0 = int(sample_validation["b0_defined_rows"].sum())
    total_dropped = int(sample_validation["dropped_undefined_rows"].sum())
    b0_fraction = total_b0 / input_rows if input_rows else np.nan
    panel_check = (
        panel_contrib.groupby("sample")["panel_fraction_of_b0"].sum().round(12).eq(1.0).all()
        if not panel_contrib.empty
        else False
    )
    example_samples = ["pt16-pre", "pt20-post", "pt69-pre"]
    examples = sample_celltype[sample_celltype["sample"].isin(example_samples)].copy()
    examples = examples.sort_values(["sample", "b0_rows"], ascending=[True, False])
    examples = examples.groupby("sample", group_keys=False).head(8)
    top_cohort = cohort_celltype.sort_values("b0_rows", ascending=False).head(20)
    near_view = near_summary.sort_values("fraction_rows_with_cross_panel_neighbor_lt_radius", ascending=False)

    lines = [
        "# WP0.3 · Full-cohort raw B0 cell table build",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.",
        "",
        "## Method",
        "",
        "Locked method: raw B0. The input `all_samples_cell_phenotypes.parquet` was filtered to drop rows with `celltype == \"undefined\"`; all remaining defined cells were pooled in the existing registered coordinate system. No cross-panel fusion, deduplication, panel normalization, CN, neighborhood, or adjacency calculation was run.",
        "",
        f"Input parquet: `{INPUT_PARQUET}`",
        f"Output parquet: `{output_parquet}`",
        "",
        "## Cohort Totals",
        "",
        markdown_table(
            pd.DataFrame(
                [
                    {
                        "union_rows": input_rows,
                        "b0_defined_rows": total_b0,
                        "dropped_undefined_rows": total_dropped,
                        "b0_over_union": b0_fraction,
                        "output_undefined_rows": output_verify["undefined_rows"],
                        "output_bad_coordinate_rows": output_verify["bad_coordinate_rows"],
                        "active_samples": output_verify["sample_count"],
                        "panel_fraction_sum_ok": panel_check,
                    }
                ]
            )
        ),
        "",
        "## Per-sample Validation",
        "",
        markdown_table(sample_validation, max_rows=90),
        "",
        "Full CSV: `deliverables/WP0.3_sample_validation.csv`.",
        "",
        "## Celltype Distribution: B0 vs Union",
        "",
        "Cohort-level celltype comparison:",
        "",
        markdown_table(top_cohort, max_rows=25),
        "",
        "Example sample-level changes (top B0 celltypes):",
        "",
        markdown_table(examples, max_rows=30),
        "",
        "Full CSVs: `deliverables/WP0.3_celltype_distribution_cohort.csv` and `deliverables/WP0.3_celltype_distribution_by_sample.csv`.",
        "",
        "## Panel Contribution",
        "",
        "Panel fractions are computed after dropping undefined rows. Each sample sums to 1.0 across panels; this confirms B0 is a single pooled defined-cell table rather than a 3x analytic denominator.",
        "",
        markdown_table(panel_contrib, max_rows=90),
        "",
        "Full CSV: `deliverables/WP0.3_panel_contribution_by_sample.csv`.",
        "",
        "## Residual Cross-panel Near-overlap",
        "",
        f"Reported only; no rows were deleted. A row is counted if it has at least one defined cell from another panel within {NEAR_OVERLAP_RADIUS_PX:.0f}px ({NEAR_OVERLAP_RADIUS_PX * UM_PER_PX:.2f} um). These are retained as distinct real cells under raw B0.",
        "",
        markdown_table(near_view, max_rows=90),
        "",
        "Full CSVs: `deliverables/WP0.3_near_overlap_by_sample.csv` and `deliverables/WP0.3_near_overlap_by_sample_pair.csv`.",
        "",
        "## Output Files",
        "",
        f"- `{output_parquet}`",
        "- `/data/home/chenwenhui/mIHC/fig6_pm/deliverables/WP0.3_B0_build_report.md`",
        "- `/data/home/chenwenhui/mIHC/fig6_pm/deliverables/WP0.3_sample_validation.csv`",
        "- `/data/home/chenwenhui/mIHC/fig6_pm/deliverables/WP0.3_celltype_distribution_cohort.csv`",
        "- `/data/home/chenwenhui/mIHC/fig6_pm/deliverables/WP0.3_celltype_distribution_by_sample.csv`",
        "- `/data/home/chenwenhui/mIHC/fig6_pm/deliverables/WP0.3_panel_contribution_by_sample.csv`",
        "- `/data/home/chenwenhui/mIHC/fig6_pm/deliverables/WP0.3_near_overlap_by_sample.csv`",
        "- `/data/home/chenwenhui/mIHC/fig6_pm/deliverables/WP0.3_near_overlap_by_sample_pair.csv`",
        "",
        "## Reproducibility",
        "",
        "- Conda env: `scimap`.",
        "- Command: `python deliverables/WP0.3_build_b0_cohort.py`.",
        f"- Batch size: {run_metadata['batch_size']:,}.",
        f"- Runtime seconds: {run_metadata['runtime_seconds']:.1f}.",
        "- Existing input parquet was read only; the new B0 parquet was written to a new output directory.",
    ]
    return "\n".join(lines) + "\n"


def run(args: argparse.Namespace) -> None:
    started = time.time()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    DELIVERABLES.mkdir(parents=True, exist_ok=True)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    if OUTPUT_PARQUET.exists() and not args.allow_existing_output:
        raise FileExistsError(f"Refusing to overwrite existing output: {OUTPUT_PARQUET}")

    metadata = pd.read_csv(SAMPLE_INFO)
    metadata_columns = [col for col in METADATA_COLUMNS if col in metadata.columns]
    pf = pq.ParquetFile(INPUT_PARQUET)
    input_rows = int(pf.metadata.num_rows)
    temp_output = OUTPUT_PARQUET.with_name(f"{OUTPUT_PARQUET.name}.tmp.{datetime.now().strftime('%Y%m%dT%H%M%S')}")

    counters = {
        "union_sample": Counter(),
        "undefined_sample": Counter(),
        "b0_sample": Counter(),
        "union_sample_celltype": Counter(),
        "b0_sample_celltype": Counter(),
        "union_celltype": Counter(),
        "b0_celltype": Counter(),
        "b0_panel": Counter(),
        "b0_bad_coord_sample": Counter(),
    }
    coord_chunks: dict[str, list[np.ndarray]] = defaultdict(list)
    panel_chunks: dict[str, list[np.ndarray]] = defaultdict(list)
    writer = None
    output_rows = 0

    try:
        for batch_i, batch in enumerate(pf.iter_batches(batch_size=args.batch_size), start=1):
            table = pa.Table.from_batches([batch])
            summary_table = table.select(["sample", "panel", "spatial_x", "spatial_y", "celltype"])
            b0_table = filter_defined_table(table)
            b0_summary_table = b0_table.select(["sample", "panel", "spatial_x", "spatial_y", "celltype"])
            summary_df = summary_table.to_pandas()
            b0_summary_df = b0_summary_table.to_pandas()
            update_batch_summaries(summary_df, b0_summary_df, counters)
            collect_coords(b0_summary_df, coord_chunks, panel_chunks)

            if b0_table.num_rows:
                b0_table = append_metadata_columns(b0_table, metadata, metadata_columns)
                if writer is None:
                    writer = pq.ParquetWriter(
                        temp_output,
                        b0_table.schema,
                        compression="zstd",
                        use_dictionary=True,
                        version="2.6",
                    )
                writer.write_table(b0_table)
                output_rows += b0_table.num_rows
            if batch_i % 20 == 0:
                print(f"processed_batches={batch_i} output_rows={output_rows}", flush=True)
    finally:
        if writer is not None:
            writer.close()

    sample_validation = sample_validation_table(counters, metadata)
    active_samples = sample_validation["sample"].nunique()
    if active_samples != EXPECTED_ACTIVE_SAMPLES:
        raise ValueError(f"Expected {EXPECTED_ACTIVE_SAMPLES} active samples, observed {active_samples}")
    if int(sample_validation["b0_bad_coordinate_rows"].sum()) != 0:
        raise ValueError("B0 output has bad coordinates before write verification")

    if OUTPUT_PARQUET.exists() and not args.allow_existing_output:
        raise FileExistsError(f"Refusing to overwrite existing output before rename: {OUTPUT_PARQUET}")
    if OUTPUT_PARQUET.exists() and args.allow_existing_output:
        OUTPUT_PARQUET.unlink()
    os.rename(temp_output, OUTPUT_PARQUET)

    output_verify = validate_output_parquet(OUTPUT_PARQUET, expected_rows=output_rows)
    if not output_verify["row_count_matches_expected"]:
        raise ValueError("Output parquet row count does not match build counter")
    if output_verify["undefined_rows"] != 0:
        raise ValueError("Output parquet still contains undefined rows")
    if output_verify["bad_coordinate_rows"] != 0:
        raise ValueError("Output parquet has bad coordinates")
    if output_verify["sample_count"] != EXPECTED_ACTIVE_SAMPLES:
        raise ValueError(f"Output sample count is {output_verify['sample_count']}, expected {EXPECTED_ACTIVE_SAMPLES}")

    sample_celltype = counter_to_celltype_df(counters["union_sample_celltype"], counters["b0_sample_celltype"], level="sample")
    cohort_celltype = counter_to_celltype_df(counters["union_celltype"], counters["b0_celltype"], level="cohort")
    panel_contrib = panel_contribution_table(counters["b0_panel"])
    panel_sums = panel_contrib.groupby("sample")["panel_fraction_of_b0"].sum().round(12)
    if not panel_sums.eq(1.0).all():
        raise ValueError("Panel contribution fractions do not sum to 1 for every sample")

    near_rows = []
    near_pair_parts = []
    for sample in sorted(coord_chunks):
        coords = np.concatenate(coord_chunks[sample], axis=0)
        panels = np.concatenate(panel_chunks[sample], axis=0)
        summary, pair_df = near_overlap_summary_for_sample(sample, coords, panels, radius_px=NEAR_OVERLAP_RADIUS_PX)
        near_rows.append(summary)
        near_pair_parts.append(pair_df)
    near_summary = pd.DataFrame(near_rows)
    near_pair = pd.concat(near_pair_parts, ignore_index=True) if near_pair_parts else pd.DataFrame()

    run_metadata = {
        "input_parquet": str(INPUT_PARQUET),
        "output_parquet": str(OUTPUT_PARQUET),
        "input_rows": input_rows,
        "output_rows": output_rows,
        "expected_active_samples": EXPECTED_ACTIVE_SAMPLES,
        "batch_size": args.batch_size,
        "near_overlap_radius_px": NEAR_OVERLAP_RADIUS_PX,
        "um_per_px": UM_PER_PX,
        "runtime_seconds": time.time() - started,
        "metadata_columns_added": metadata_columns,
        "output_verify": output_verify,
    }

    sample_validation.to_csv(DELIVERABLES / "WP0.3_sample_validation.csv", index=False)
    cohort_celltype.sort_values("b0_rows", ascending=False).to_csv(
        DELIVERABLES / "WP0.3_celltype_distribution_cohort.csv", index=False
    )
    sample_celltype.sort_values(["sample", "b0_rows"], ascending=[True, False]).to_csv(
        DELIVERABLES / "WP0.3_celltype_distribution_by_sample.csv", index=False
    )
    panel_contrib.to_csv(DELIVERABLES / "WP0.3_panel_contribution_by_sample.csv", index=False)
    near_summary.to_csv(DELIVERABLES / "WP0.3_near_overlap_by_sample.csv", index=False)
    near_pair.to_csv(DELIVERABLES / "WP0.3_near_overlap_by_sample_pair.csv", index=False)
    (DELIVERABLES / "WP0.3_run_metadata.json").write_text(
        json.dumps(run_metadata, indent=2, ensure_ascii=False, default=json_default),
        encoding="utf-8",
    )

    report = make_report(
        output_parquet=OUTPUT_PARQUET,
        input_rows=input_rows,
        output_verify=output_verify,
        sample_validation=sample_validation,
        cohort_celltype=cohort_celltype,
        sample_celltype=sample_celltype,
        panel_contrib=panel_contrib,
        near_summary=near_summary,
        run_metadata=run_metadata,
    )
    (DELIVERABLES / "WP0.3_B0_build_report.md").write_text(report, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build full-cohort raw B0 cell table for Fig. 6.")
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    parser.add_argument("--allow-existing-output", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
