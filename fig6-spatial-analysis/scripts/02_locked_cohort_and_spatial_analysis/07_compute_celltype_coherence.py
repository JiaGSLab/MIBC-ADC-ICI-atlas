#!/usr/bin/env python3
"""WP1.5 celltype proximity and coherence analysis."""

from __future__ import annotations

import argparse
import json
import math
import re
import time
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from scipy.spatial import cKDTree


PROJECT_DIR = Path("/data/home/chenwenhui/mIHC/fig6_pm")
CELLINFO_DIR = Path("/data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh")
DELIVERABLE_DIR = PROJECT_DIR / "deliverables"
SCRATCH_DIR = PROJECT_DIR / "scratch_wp1.5"

B0_TABLE = CELLINFO_DIR / "all_samples_cell_phenotypes_B0.parquet"
SAMPLE_VALIDATION = DELIVERABLE_DIR / "WP0.3_sample_validation.csv"
CN_QC = DELIVERABLE_DIR / "WP1.1_cn_qc_by_sample.csv"

OUT_PROXIMITY = DELIVERABLE_DIR / "WP1.5_proximity_sample_metrics.csv"
OUT_COHERENCE = DELIVERABLE_DIR / "WP1.5_coherence_sample_metrics.csv"
OUT_STATS = DELIVERABLE_DIR / "WP1.5_CRvsNR_stats.csv"
OUT_REPORT = DELIVERABLE_DIR / "WP1.5_proximity_coherence_report.md"

PRIMARY_LOW_CELL_THRESHOLD = 30000
SENSITIVITY_LOW_CELL_THRESHOLD = 20000
PRESENCE_MIN_PER_GROUP = 5
DEFAULT_BOOTSTRAP = 5000
DEFAULT_PERMUTATIONS = 25
DEFAULT_PERM_SOURCE_CAP = 20000
RADII = (20, 30, 50)
KNN_NEIGHBORS = 6

V1_LIKE_SAMPLES = {
    "pt16-pre",
    "pt17-pre",
    "pt19-pre",
    "pt20-pre",
    "pt21-pre",
    "pt23-post",
    "pt23-pre",
}
V1_LAYOUT = "v1_like_missing_SMA_CD3_PDL1"
V2_LAYOUT = "v2_like"

CELLTYPE_SETS = {
    "HER2+ Tumor": ("HER2+ Tumor",),
    "PanCK+ Tumor": ("PanCK+ Tumor",),
    "Vascular": ("Vascular",),
    "CD8 T cell": ("CD8 T cell",),
    "Treg": ("Treg",),
    "Macrophage": ("CD163+ Macrophage", "CD163- Macrophage"),
    "CD21+ B cell": ("CD21+ B cell",),
    "CD163+ Macrophage": ("CD163+ Macrophage",),
    "CD163- Macrophage": ("CD163- Macrophage",),
}


@dataclass(frozen=True)
class PairDef:
    pair_id: str
    source_label: str
    target_label: str
    pair_group: str
    same_panel_only: bool = False
    panel_context: str = "pooled_3D_collapsed"


PAIR_DEFS = [
    PairDef(
        "HER2_Tumor_to_Vascular_same_panel",
        "HER2+ Tumor",
        "Vascular",
        "headline_same_panel",
        same_panel_only=True,
        panel_context="same_panel_single_section",
    ),
    PairDef(
        "PanCK_Tumor_to_Vascular_same_panel",
        "PanCK+ Tumor",
        "Vascular",
        "headline_same_panel",
        same_panel_only=True,
        panel_context="same_panel_single_section",
    ),
    PairDef(
        "HER2_Tumor_to_Vascular_3D_all",
        "HER2+ Tumor",
        "Vascular",
        "headline_3D_sensitivity",
        same_panel_only=False,
    ),
    PairDef(
        "PanCK_Tumor_to_Vascular_3D_all",
        "PanCK+ Tumor",
        "Vascular",
        "headline_3D_sensitivity",
        same_panel_only=False,
    ),
    PairDef("HER2_Tumor_to_CD8", "HER2+ Tumor", "CD8 T cell", "general_proximity"),
    PairDef("PanCK_Tumor_to_CD8", "PanCK+ Tumor", "CD8 T cell", "general_proximity"),
    PairDef("HER2_Tumor_to_Treg", "HER2+ Tumor", "Treg", "general_proximity"),
    PairDef("PanCK_Tumor_to_Treg", "PanCK+ Tumor", "Treg", "general_proximity"),
    PairDef("Macrophage_to_CD8", "Macrophage", "CD8 T cell", "general_proximity"),
    PairDef("CD8_to_Treg", "CD8 T cell", "Treg", "general_proximity"),
]

COHERENCE_CELLTYPES = [
    "CD8 T cell",
    "CD21+ B cell",
    "CD163+ Macrophage",
    "CD163- Macrophage",
    "HER2+ Tumor",
    "PanCK+ Tumor",
    "Vascular",
]


@dataclass(frozen=True)
class AnalysisConfig:
    name: str
    low_cell_threshold: int | None
    exclude_v1_all_features: bool = False


ANALYSIS_CONFIGS = [
    AnalysisConfig("primary_threshold30k_panelaware", PRIMARY_LOW_CELL_THRESHOLD),
    AnalysisConfig("sensitivity_include_lowcells_panelaware", None),
    AnalysisConfig("sensitivity_threshold20k_panelaware", SENSITIVITY_LOW_CELL_THRESHOLD),
    AnalysisConfig("sensitivity_threshold30k_exclude_v1_all", PRIMARY_LOW_CELL_THRESHOLD, True),
]


def log(message: str) -> None:
    print(message, flush=True)


def celltype_mask(celltypes: pd.Series | np.ndarray, label: str) -> np.ndarray:
    values = np.asarray(celltypes).astype(object)
    allowed = CELLTYPE_SETS.get(label, (label,))
    return np.isin(values, allowed)


def nearest_neighbor_summary(source_coords: np.ndarray, target_coords: np.ndarray, radii: Iterable[int] = RADII) -> dict:
    source_coords = np.asarray(source_coords, dtype=float)
    target_coords = np.asarray(target_coords, dtype=float)
    out = {
        "source_count": int(len(source_coords)),
        "target_count": int(len(target_coords)),
        "median_nn_distance": np.nan,
        "mean_nn_distance": np.nan,
    }
    for radius in radii:
        out[f"cooccurrence_prop_r{int(radius)}"] = np.nan
    if len(source_coords) == 0 or len(target_coords) == 0:
        return out
    tree = cKDTree(target_coords)
    distances, _ = tree.query(source_coords, k=1, workers=-1)
    distances = np.asarray(distances, dtype=float)
    out["median_nn_distance"] = float(np.nanmedian(distances))
    out["mean_nn_distance"] = float(np.nanmean(distances))
    for radius in radii:
        out[f"cooccurrence_prop_r{int(radius)}"] = float(np.mean(distances <= float(radius)))
    return out


def permutation_cooccurrence_null(
    source_coords: np.ndarray,
    universe_coords: np.ndarray,
    target_count: int,
    radii: Iterable[int],
    *,
    n_permutations: int,
    random_seed: int,
    source_cap: int,
) -> dict:
    radii = tuple(int(r) for r in radii)
    out = {}
    for radius in radii:
        out[f"null_mean_cooccurrence_prop_r{radius}"] = np.nan
        out[f"null_sd_cooccurrence_prop_r{radius}"] = np.nan
        out[f"cooccurrence_enrichment_r{radius}"] = np.nan
        out[f"cooccurrence_z_r{radius}"] = np.nan
    source_coords = np.asarray(source_coords, dtype=float)
    universe_coords = np.asarray(universe_coords, dtype=float)
    if (
        n_permutations <= 0
        or len(source_coords) == 0
        or len(universe_coords) == 0
        or target_count <= 0
        or target_count > len(universe_coords)
    ):
        return out
    rng = np.random.default_rng(random_seed)
    if len(source_coords) > int(source_cap):
        source_coords = source_coords[rng.choice(len(source_coords), size=int(source_cap), replace=False)]
    null = np.empty((int(n_permutations), len(radii)), dtype=float)
    for i in range(int(n_permutations)):
        idx = rng.choice(len(universe_coords), size=int(target_count), replace=False)
        tree = cKDTree(universe_coords[idx])
        distances, _ = tree.query(source_coords, k=1, workers=-1)
        for j, radius in enumerate(radii):
            null[i, j] = float(np.mean(distances <= radius))
    for j, radius in enumerate(radii):
        vals = null[:, j]
        out[f"null_mean_cooccurrence_prop_r{radius}"] = float(np.mean(vals))
        out[f"null_sd_cooccurrence_prop_r{radius}"] = float(np.std(vals, ddof=1)) if len(vals) > 1 else np.nan
    return out


def _combine_panel_null(null_parts: list[tuple[int, dict]], radii: Iterable[int]) -> dict:
    total_source = sum(n for n, _ in null_parts)
    out = {}
    for radius in radii:
        means = [n * d.get(f"null_mean_cooccurrence_prop_r{radius}", np.nan) for n, d in null_parts]
        sds = [d.get(f"null_sd_cooccurrence_prop_r{radius}", np.nan) for _, d in null_parts]
        if total_source and np.isfinite(means).any():
            out[f"null_mean_cooccurrence_prop_r{radius}"] = float(np.nansum(means) / total_source)
        else:
            out[f"null_mean_cooccurrence_prop_r{radius}"] = np.nan
        finite_sds = [x for x in sds if np.isfinite(x)]
        out[f"null_sd_cooccurrence_prop_r{radius}"] = float(np.mean(finite_sds)) if finite_sds else np.nan
    return out


def compute_pair_metrics_for_sample(
    sample_df: pd.DataFrame,
    pair: PairDef,
    radii: Iterable[int] = RADII,
    *,
    n_permutations: int,
    random_seed: int,
    perm_source_cap: int = DEFAULT_PERM_SOURCE_CAP,
) -> dict:
    radii = tuple(int(r) for r in radii)
    df = sample_df
    source_mask = celltype_mask(df["celltype"], pair.source_label)
    target_mask = celltype_mask(df["celltype"], pair.target_label)
    all_source_count = int(source_mask.sum())
    all_target_count = int(target_mask.sum())
    record = {
        "sample": str(df["sample"].iloc[0]) if len(df) else "",
        "pair_id": pair.pair_id,
        "source_label": pair.source_label,
        "target_label": pair.target_label,
        "pair_group": pair.pair_group,
        "same_panel_only": bool(pair.same_panel_only),
        "panel_context": pair.panel_context,
        "all_source_count": all_source_count,
        "all_target_count": all_target_count,
        "shared_panel_count": 0,
        "panels_used": "",
        "status": "ok",
    }
    if pair.same_panel_only:
        source_panels = set(df.loc[source_mask, "panel"].astype(str))
        target_panels = set(df.loc[target_mask, "panel"].astype(str))
        shared_panels = sorted(source_panels & target_panels)
        record["shared_panel_count"] = int(len(shared_panels))
        record["panels_used"] = ";".join(shared_panels)
        if not shared_panels:
            record.update(nearest_neighbor_summary(np.empty((0, 2)), np.empty((0, 2)), radii))
            record["all_source_count"] = all_source_count
            record["all_target_count"] = all_target_count
            record["status"] = "no_shared_panel"
            for radius in radii:
                record[f"null_mean_cooccurrence_prop_r{radius}"] = np.nan
                record[f"null_sd_cooccurrence_prop_r{radius}"] = np.nan
                record[f"cooccurrence_enrichment_r{radius}"] = np.nan
                record[f"cooccurrence_z_r{radius}"] = np.nan
            return record
        source_coords_parts = []
        target_coords_parts = []
        null_parts = []
        for panel in shared_panels:
            panel_mask = df["panel"].astype(str).eq(panel)
            src = df.loc[panel_mask & source_mask, ["spatial_x", "spatial_y"]].to_numpy(dtype=float)
            tgt = df.loc[panel_mask & target_mask, ["spatial_x", "spatial_y"]].to_numpy(dtype=float)
            uni = df.loc[panel_mask, ["spatial_x", "spatial_y"]].to_numpy(dtype=float)
            if len(src) and len(tgt):
                source_coords_parts.append(src)
                target_coords_parts.append(tgt)
                null = permutation_cooccurrence_null(
                    src,
                    uni,
                    len(tgt),
                    radii,
                    n_permutations=n_permutations,
                    random_seed=random_seed + zlib.crc32(panel.encode("utf-8")),
                    source_cap=perm_source_cap,
                )
                null_parts.append((len(src), null))
        source_coords = np.vstack(source_coords_parts) if source_coords_parts else np.empty((0, 2))
        target_coords = np.vstack(target_coords_parts) if target_coords_parts else np.empty((0, 2))
        universe_null = _combine_panel_null(null_parts, radii)
    else:
        source_coords = df.loc[source_mask, ["spatial_x", "spatial_y"]].to_numpy(dtype=float)
        target_coords = df.loc[target_mask, ["spatial_x", "spatial_y"]].to_numpy(dtype=float)
        universe_coords = df[["spatial_x", "spatial_y"]].to_numpy(dtype=float)
        universe_null = permutation_cooccurrence_null(
            source_coords,
            universe_coords,
            len(target_coords),
            radii,
            n_permutations=n_permutations,
            random_seed=random_seed,
            source_cap=perm_source_cap,
        )
    record.update(nearest_neighbor_summary(source_coords, target_coords, radii))
    record["source_count"] = int(len(source_coords))
    record["target_count"] = int(len(target_coords))
    record.update(universe_null)
    if len(source_coords) == 0 or len(target_coords) == 0:
        record["status"] = "missing_source_or_target"
    for radius in radii:
        obs = record.get(f"cooccurrence_prop_r{radius}", np.nan)
        null_mean = record.get(f"null_mean_cooccurrence_prop_r{radius}", np.nan)
        null_sd = record.get(f"null_sd_cooccurrence_prop_r{radius}", np.nan)
        record[f"cooccurrence_enrichment_r{radius}"] = obs / null_mean if null_mean and np.isfinite(null_mean) else np.nan
        record[f"cooccurrence_z_r{radius}"] = (
            (obs - null_mean) / null_sd if null_sd and np.isfinite(null_sd) and null_sd > 0 else np.nan
        )
    return record


def make_tumor_vascular_contrasts(proximity_df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for sample, sdf in proximity_df.groupby("sample"):
        her2 = sdf[sdf["pair_id"].eq("HER2_Tumor_to_Vascular_same_panel")]
        panck = sdf[sdf["pair_id"].eq("PanCK_Tumor_to_Vascular_same_panel")]
        if her2.empty or panck.empty:
            continue
        her2 = her2.iloc[0]
        panck = panck.iloc[0]
        delta_dist = her2.get("median_nn_distance", np.nan) - panck.get("median_nn_distance", np.nan)
        delta_r50 = her2.get("cooccurrence_prop_r50", np.nan) - panck.get("cooccurrence_prop_r50", np.nan)
        closer = "tie_or_NA"
        if np.isfinite(delta_dist):
            closer = "HER2+ Tumor" if delta_dist < 0 else "PanCK+ Tumor" if delta_dist > 0 else "tie"
        records.append(
            {
                "sample": sample,
                "pair_id": "HER2_minus_PanCK_to_Vascular_same_panel",
                "source_label": "HER2+ Tumor minus PanCK+ Tumor",
                "target_label": "Vascular",
                "pair_group": "headline_contrast",
                "same_panel_only": True,
                "panel_context": "same_panel_single_section",
                "delta_median_nn_distance_HER2_minus_PanCK": delta_dist,
                "delta_cooccurrence_prop_r50_HER2_minus_PanCK": delta_r50,
                "closer_to_vascular": closer,
                "status": "ok" if np.isfinite(delta_dist) else "missing_strict_same_panel_metric",
            }
        )
    return pd.DataFrame(records)


def compute_knn_coherence(
    sample: str,
    coords: np.ndarray,
    labels: np.ndarray,
    target_labels: Iterable[str],
    *,
    k_neighbors: int = KNN_NEIGHBORS,
) -> pd.DataFrame:
    coords = np.asarray(coords, dtype=float)
    labels = np.asarray(labels).astype(object)
    n = len(coords)
    if n <= 1:
        return pd.DataFrame()
    k_eff = min(int(k_neighbors), n - 1)
    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k_eff + 1, workers=-1)
    neighbor_idx = idx[:, 1:]
    unique_labels, label_codes = np.unique(labels.astype(str), return_inverse=True)
    code_by_label = {label: i for i, label in enumerate(unique_labels)}
    neighbor_codes = label_codes[neighbor_idx]
    records = []
    w = float(n * k_eff)
    for label in target_labels:
        if label not in code_by_label:
            records.append(
                {
                    "sample": sample,
                    "celltype": label,
                    "cell_count": 0,
                    "k_neighbors": k_eff,
                    "same_type_knn_fraction": np.nan,
                    "same_type_knn_expected": np.nan,
                    "same_type_knn_enrichment": np.nan,
                    "binary_moran_i_knn": np.nan,
                }
            )
            continue
        code = code_by_label[label]
        mask = label_codes == code
        count = int(mask.sum())
        same_fraction = float(np.mean(neighbor_codes[mask] == code)) if count else np.nan
        expected = (count - 1) / (n - 1) if n > 1 else np.nan
        x = mask.astype(float)
        x_mean = float(x.mean())
        centered = x - x_mean
        denom = float(np.sum(centered * centered))
        moran = np.nan
        if denom > 0 and w > 0:
            moran_num = float(np.sum(centered[:, None] * centered[neighbor_idx]))
            moran = (n / w) * (moran_num / denom)
        records.append(
            {
                "sample": sample,
                "celltype": label,
                "cell_count": count,
                "k_neighbors": k_eff,
                "same_type_knn_fraction": same_fraction,
                "same_type_knn_expected": expected,
                "same_type_knn_enrichment": same_fraction / expected if expected and np.isfinite(expected) else np.nan,
                "binary_moran_i_knn": moran,
            }
        )
    return pd.DataFrame(records)


def read_sample_meta() -> pd.DataFrame:
    validation = pd.read_csv(SAMPLE_VALIDATION)
    qc = pd.read_csv(CN_QC)
    qc_cols = [
        "sample",
        "cn_participant_cells",
        "cn_filtered_cells",
        "cn_participant_fraction",
        "low_cell_lt30000",
    ]
    meta = validation.merge(qc[qc_cols], on="sample", how="left", validate="one_to_one")
    meta["panel_layout"] = np.where(meta["sample"].isin(V1_LIKE_SAMPLES), V1_LAYOUT, V2_LAYOUT)
    meta["is_v1_like"] = meta["sample"].isin(V1_LIKE_SAMPLES)
    meta["passes_lowcell_30000"] = meta["cn_participant_cells"].ge(PRIMARY_LOW_CELL_THRESHOLD)
    meta["passes_lowcell_20000"] = meta["cn_participant_cells"].ge(SENSITIVITY_LOW_CELL_THRESHOLD)
    return meta


def process_one_sample(
    sample: str,
    sample_df: pd.DataFrame,
    *,
    n_permutations: int,
    random_seed: int,
    perm_source_cap: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_seed = (random_seed + zlib.crc32(sample.encode("utf-8"))) & 0xFFFFFFFF
    prox_records = []
    for i, pair in enumerate(PAIR_DEFS):
        prox_records.append(
            compute_pair_metrics_for_sample(
                sample_df,
                pair,
                RADII,
                n_permutations=n_permutations,
                random_seed=sample_seed + i * 1009,
                perm_source_cap=perm_source_cap,
            )
        )
    proximity = pd.DataFrame(prox_records)
    contrasts = make_tumor_vascular_contrasts(proximity)
    if not contrasts.empty:
        proximity = pd.concat([proximity, contrasts], ignore_index=True, sort=False)
    coords = sample_df[["spatial_x", "spatial_y"]].to_numpy(dtype=float)
    labels = sample_df["celltype"].astype(str).to_numpy()
    coherence = compute_knn_coherence(sample, coords, labels, COHERENCE_CELLTYPES, k_neighbors=KNN_NEIGHBORS)
    return proximity, coherence


def cohens_d_nr_minus_cr(cr_values: Iterable[float], nr_values: Iterable[float]) -> float:
    cr = np.asarray(list(cr_values), dtype=float)
    nr = np.asarray(list(nr_values), dtype=float)
    cr = cr[np.isfinite(cr)]
    nr = nr[np.isfinite(nr)]
    if len(cr) == 0 or len(nr) == 0:
        return np.nan
    diff = float(np.mean(nr) - np.mean(cr))
    if len(cr) + len(nr) < 3:
        return np.nan
    var_cr = np.var(cr, ddof=1) if len(cr) > 1 else 0.0
    var_nr = np.var(nr, ddof=1) if len(nr) > 1 else 0.0
    pooled_var = ((len(cr) - 1) * var_cr + (len(nr) - 1) * var_nr) / (len(cr) + len(nr) - 2)
    if pooled_var <= 0 or not np.isfinite(pooled_var):
        return 0.0 if diff == 0 else np.nan
    return diff / math.sqrt(pooled_var)


def bootstrap_cohens_d_ci(cr_values, nr_values, *, n_boot: int, seed: int) -> tuple[float, float]:
    cr = np.asarray(list(cr_values), dtype=float)
    nr = np.asarray(list(nr_values), dtype=float)
    cr = cr[np.isfinite(cr)]
    nr = nr[np.isfinite(nr)]
    if len(cr) == 0 or len(nr) == 0:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    boot = np.empty(int(n_boot), dtype=float)
    for i in range(int(n_boot)):
        boot[i] = cohens_d_nr_minus_cr(
            rng.choice(cr, size=len(cr), replace=True),
            rng.choice(nr, size=len(nr), replace=True),
        )
    boot = boot[np.isfinite(boot)]
    if len(boot) == 0:
        return (np.nan, np.nan)
    return (float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975)))


def mannwhitney_p(cr_values, nr_values) -> float:
    cr = np.asarray(list(cr_values), dtype=float)
    nr = np.asarray(list(nr_values), dtype=float)
    cr = cr[np.isfinite(cr)]
    nr = nr[np.isfinite(nr)]
    if len(cr) == 0 or len(nr) == 0:
        return np.nan
    try:
        return float(scipy_stats.mannwhitneyu(nr, cr, alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def bh_fdr(p_values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(p_values), dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return q
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    q[valid] = restored
    return q


def proximity_feature_table(proximity: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    pair_metric_cols = [
        "median_nn_distance",
        "mean_nn_distance",
        "cooccurrence_prop_r20",
        "cooccurrence_prop_r30",
        "cooccurrence_prop_r50",
        "cooccurrence_enrichment_r20",
        "cooccurrence_enrichment_r30",
        "cooccurrence_enrichment_r50",
        "cooccurrence_z_r20",
        "cooccurrence_z_r30",
        "cooccurrence_z_r50",
    ]
    contrast_metric_cols = [
        "delta_median_nn_distance_HER2_minus_PanCK",
        "delta_cooccurrence_prop_r50_HER2_minus_PanCK",
    ]
    records = []
    for _, row in proximity.iterrows():
        metric_cols = contrast_metric_cols if row.get("pair_group", "") == "headline_contrast" else pair_metric_cols
        for metric in metric_cols:
            if metric not in proximity.columns:
                continue
            value = row.get(metric, np.nan)
            records.append(
                {
                    "sample": row["sample"],
                    "feature_class": "proximity",
                    "metric": metric,
                    "feature_id": row["pair_id"],
                    "feature_name": f"{row.get('source_label', '')} -> {row.get('target_label', '')}",
                    "pair_group": row.get("pair_group", ""),
                    "panel_context": row.get("panel_context", ""),
                    "same_panel_only": bool(row.get("same_panel_only", False)),
                    "status": row.get("status", ""),
                    "value": value,
                }
            )
    features = pd.DataFrame(records)
    return features.merge(meta, on="sample", how="left", validate="many_to_one")


def coherence_feature_table(coherence: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    metric_cols = ["same_type_knn_fraction", "same_type_knn_enrichment", "binary_moran_i_knn"]
    records = []
    for _, row in coherence.iterrows():
        for metric in metric_cols:
            records.append(
                {
                    "sample": row["sample"],
                    "feature_class": "coherence",
                    "metric": metric,
                    "feature_id": row["celltype"],
                    "feature_name": row["celltype"],
                    "pair_group": "coherence",
                    "panel_context": "pooled_3D_collapsed_knn",
                    "same_panel_only": False,
                    "status": "ok" if row.get("cell_count", 0) > 0 else "missing_celltype",
                    "value": row.get(metric, np.nan),
                }
            )
    features = pd.DataFrame(records)
    return features.merge(meta, on="sample", how="left", validate="many_to_one")


def _filter_for_config(features: pd.DataFrame, config: AnalysisConfig) -> pd.Series:
    mask = features["Time_short"].eq("Pre") & features["Outcome_short"].isin(["CR", "NR"])
    if config.low_cell_threshold is not None:
        mask &= features["cn_participant_cells"].ge(config.low_cell_threshold)
    if config.exclude_v1_all_features:
        mask &= features["panel_layout"].ne(V1_LAYOUT)
    return mask


def compute_crnr_stats(
    features: pd.DataFrame,
    *,
    bootstrap_n: int = DEFAULT_BOOTSTRAP,
    min_group_presence: int = PRESENCE_MIN_PER_GROUP,
    random_seed: int = 0,
) -> pd.DataFrame:
    records = []
    group_cols = [
        "feature_class",
        "metric",
        "feature_id",
        "feature_name",
        "pair_group",
        "panel_context",
        "same_panel_only",
    ]
    for config in ANALYSIS_CONFIGS:
        analysis = features.loc[_filter_for_config(features, config)].copy()
        for key, sub in analysis.groupby(group_cols, dropna=False, sort=False):
            feature_class, metric, feature_id, feature_name, pair_group, panel_context, same_panel_only = key
            values = sub[["Outcome_short", "value", "status", "panel_layout"]].copy()
            finite = np.isfinite(pd.to_numeric(values["value"], errors="coerce"))
            values = values.loc[finite]
            cr = values.loc[values["Outcome_short"].eq("CR"), "value"].astype(float).to_numpy()
            nr = values.loc[values["Outcome_short"].eq("NR"), "value"].astype(float).to_numpy()
            presence_pass = len(cr) >= int(min_group_presence) and len(nr) >= int(min_group_presence)
            seed = (
                random_seed
                + zlib.crc32(f"{config.name}|{feature_class}|{metric}|{feature_id}".encode("utf-8"))
            ) & 0xFFFFFFFF
            d = cohens_d_nr_minus_cr(cr, nr)
            ci_low, ci_high = bootstrap_cohens_d_ci(cr, nr, n_boot=bootstrap_n, seed=seed)
            p_value = mannwhitney_p(cr, nr)
            if np.isfinite(d):
                direction = "higher_in_NR" if d > 0 else "higher_in_CR" if d < 0 else "no_difference"
            else:
                diff = np.nanmean(nr) - np.nanmean(cr) if len(cr) and len(nr) else np.nan
                direction = "higher_in_NR" if diff > 0 else "higher_in_CR" if diff < 0 else "no_difference"
            if config.exclude_v1_all_features:
                panel_handling = "excluded_v1_all_features"
            elif same_panel_only:
                panel_handling = "same_panel_strict_nonshared_panels_as_missing"
            else:
                panel_handling = "included_all_panel_layouts"
            records.append(
                {
                    "analysis_set": config.name,
                    "low_cell_threshold": config.low_cell_threshold if config.low_cell_threshold is not None else "none",
                    "panel_handling": panel_handling,
                    "feature_class": feature_class,
                    "metric": metric,
                    "feature_id": feature_id,
                    "feature_name": feature_name,
                    "pair_group": pair_group,
                    "panel_context": panel_context,
                    "same_panel_only": bool(same_panel_only),
                    "presence_n_CR": int(len(cr)),
                    "presence_n_NR": int(len(nr)),
                    "presence_filter_min_per_group": int(min_group_presence),
                    "presence_filter_pass": bool(presence_pass),
                    "mean_CR": float(np.nanmean(cr)) if len(cr) else np.nan,
                    "mean_NR": float(np.nanmean(nr)) if len(nr) else np.nan,
                    "median_CR": float(np.nanmedian(cr)) if len(cr) else np.nan,
                    "median_NR": float(np.nanmedian(nr)) if len(nr) else np.nan,
                    "diff_NR_minus_CR": float(np.nanmean(nr) - np.nanmean(cr)) if len(cr) and len(nr) else np.nan,
                    "cohens_d_NR_minus_CR": d,
                    "cohens_d_bootstrap95_low": ci_low,
                    "cohens_d_bootstrap95_high": ci_high,
                    "wilcoxon_mannwhitney_p": p_value,
                    "direction": direction,
                }
            )
    stats_df = pd.DataFrame(records)
    stats_df["FDR_BH_within_analysis_set"] = np.nan
    for _, idx in stats_df.groupby("analysis_set").groups.items():
        idx = np.array(list(idx))
        pass_mask = stats_df.loc[idx, "presence_filter_pass"].astype(bool).to_numpy()
        idx_pass = idx[pass_mask]
        if len(idx_pass):
            stats_df.loc[idx_pass, "FDR_BH_within_analysis_set"] = bh_fdr(
                stats_df.loc[idx_pass, "wilcoxon_mannwhitney_p"].to_numpy()
            )
    stats_df["passes_FDR_0.10"] = stats_df["FDR_BH_within_analysis_set"].lt(0.10)
    return stats_df.sort_values(
        ["analysis_set", "FDR_BH_within_analysis_set", "feature_class", "metric", "feature_id"],
        na_position="last",
    ).reset_index(drop=True)


def fmt_float(value: float, digits: int = 3) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def fmt_p(value: float) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    if value < 0.001:
        return f"{value:.2e}"
    return f"{value:.3f}"


def markdown_table(df: pd.DataFrame, columns: list[str], max_rows: int = 12) -> str:
    if df.empty:
        return "None."
    view = df.loc[:, columns].head(max_rows).copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            if "p" in col.lower() or "fdr" in col.lower():
                view[col] = view[col].map(fmt_p)
            else:
                view[col] = view[col].map(lambda x: fmt_float(x, 3))
    view = view.fillna("NA").astype(str).map(lambda x: x.replace("|", "\\|").replace("\n", " "))
    return "\n".join(
        [
            "| " + " | ".join(view.columns) + " |",
            "| " + " | ".join(["---"] * len(view.columns)) + " |",
            *["| " + " | ".join(row) + " |" for row in view.to_numpy()],
        ]
    )


def panel_audit(proximity: pd.DataFrame) -> dict:
    her2 = proximity[proximity["pair_id"].eq("HER2_Tumor_to_Vascular_same_panel")]
    panck = proximity[proximity["pair_id"].eq("PanCK_Tumor_to_Vascular_same_panel")]
    return {
        "HER2_Vascular_shared_panel_samples": int(her2["shared_panel_count"].gt(0).sum()) if not her2.empty else 0,
        "PanCK_Vascular_shared_panel_samples": int(panck["shared_panel_count"].gt(0).sum()) if not panck.empty else 0,
        "HER2_no_shared_panel_samples": her2.loc[her2["shared_panel_count"].eq(0), "sample"].tolist() if not her2.empty else [],
    }


def build_report(proximity: pd.DataFrame, coherence: pd.DataFrame, stats_df: pd.DataFrame) -> str:
    audit = panel_audit(proximity)
    primary = stats_df[stats_df["analysis_set"].eq("primary_threshold30k_panelaware")].copy()
    primary["abs_d"] = primary["cohens_d_NR_minus_CR"].abs()
    headline = primary[
        primary["pair_group"].isin(["headline_same_panel", "headline_contrast", "headline_3D_sensitivity"])
        & primary["metric"].isin(
            [
                "median_nn_distance",
                "cooccurrence_prop_r50",
                "cooccurrence_z_r50",
                "delta_median_nn_distance_HER2_minus_PanCK",
                "delta_cooccurrence_prop_r50_HER2_minus_PanCK",
            ]
        )
    ].sort_values(["pair_group", "metric", "feature_id"])
    sig = primary[primary["FDR_BH_within_analysis_set"].lt(0.10)]
    top_general = primary[
        primary["pair_group"].eq("general_proximity")
        & primary["metric"].isin(["median_nn_distance", "cooccurrence_z_r50", "cooccurrence_prop_r50"])
    ].sort_values("abs_d", ascending=False)
    top_coherence = primary[primary["feature_class"].eq("coherence")].sort_values("abs_d", ascending=False)
    summary = (
        stats_df.groupby("analysis_set")
        .agg(
            n_features=("feature_id", "size"),
            n_features_in_FDR=("presence_filter_pass", "sum"),
            n_FDR_lt_0_10=("passes_FDR_0.10", "sum"),
            min_FDR=("FDR_BH_within_analysis_set", "min"),
        )
        .reset_index()
        .sort_values("analysis_set")
    )
    lines = [
        "# WP1.5 Celltype Proximity and Coherence",
        "",
        "## Panel Audit",
        "",
        "The intended same-panel tumor-to-vascular headline is only fully true for PanCK+ Tumor -> Vascular. "
        f"PanCK+ Tumor and Vascular share a panel in {audit['PanCK_Vascular_shared_panel_samples']}/72 samples. "
        f"HER2+ Tumor and Vascular share a panel in {audit['HER2_Vascular_shared_panel_samples']}/72 samples; "
        f"the non-shared samples are {', '.join(audit['HER2_no_shared_panel_samples']) or 'none'}. "
        "Therefore the primary HER2+ Tumor -> Vascular single-section metric treats non-shared-panel samples as missing, "
        "and an all-sample pooled 3D sensitivity is reported separately.",
        "",
        "## Primary Headline Results",
        "",
        "No primary WP1.5 feature passed BH-FDR < 0.10." if sig.empty else f"{len(sig)} primary features passed BH-FDR < 0.10.",
        "",
        markdown_table(
            headline,
            [
                "pair_group",
                "metric",
                "feature_id",
                "presence_n_CR",
                "presence_n_NR",
                "mean_CR",
                "mean_NR",
                "cohens_d_NR_minus_CR",
                "wilcoxon_mannwhitney_p",
                "FDR_BH_within_analysis_set",
                "direction",
                "panel_handling",
            ],
            30,
        ),
        "",
        "## General Proximity",
        "",
        markdown_table(
            top_general,
            [
                "metric",
                "feature_id",
                "presence_n_CR",
                "presence_n_NR",
                "mean_CR",
                "mean_NR",
                "cohens_d_NR_minus_CR",
                "wilcoxon_mannwhitney_p",
                "FDR_BH_within_analysis_set",
                "direction",
            ],
            15,
        ),
        "",
        "## Coherence",
        "",
        markdown_table(
            top_coherence,
            [
                "metric",
                "feature_id",
                "presence_n_CR",
                "presence_n_NR",
                "mean_CR",
                "mean_NR",
                "cohens_d_NR_minus_CR",
                "wilcoxon_mannwhitney_p",
                "FDR_BH_within_analysis_set",
                "direction",
            ],
            15,
        ),
        "",
        "## Sensitivity Summary",
        "",
        markdown_table(summary, ["analysis_set", "n_features", "n_features_in_FDR", "n_FDR_lt_0_10", "min_FDR"], 20),
        "",
        "## Outputs",
        "",
        f"- `{OUT_PROXIMITY}`",
        f"- `{OUT_COHERENCE}`",
        f"- `{OUT_STATS}`",
        f"- `{OUT_REPORT}`",
        "",
    ]
    return "\n".join(lines)


def write_metadata(path: Path, input_path: Path, command: str) -> None:
    meta = {
        "output": str(path),
        "input": str(input_path),
        "generated_at_unix": time.time(),
        "command": command,
    }
    if path.suffix == ".csv":
        df = pd.read_csv(path)
        meta["row_count"] = int(df.shape[0])
        meta["column_count"] = int(df.shape[1])
        meta["columns"] = list(map(str, df.columns))
    else:
        text = path.read_text(encoding="utf-8")
        meta["byte_count"] = int(path.stat().st_size)
        meta["line_count"] = int(text.count("\n") + 1)
    path.with_suffix(path.suffix + ".metadata.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="WP1.5 celltype proximity and coherence analysis")
    parser.add_argument("--n-permutations", type=int, default=DEFAULT_PERMUTATIONS)
    parser.add_argument("--bootstrap-n", type=int, default=DEFAULT_BOOTSTRAP)
    parser.add_argument("--random-seed", type=int, default=15)
    parser.add_argument("--perm-source-cap", type=int, default=DEFAULT_PERM_SOURCE_CAP)
    parser.add_argument("--samples", nargs="*", default=None)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    DELIVERABLE_DIR.mkdir(parents=True, exist_ok=True)
    intermediate_dir = SCRATCH_DIR / "sample_intermediates"
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    meta = read_sample_meta()

    columns = ["sample", "panel", "celltype", "spatial_x", "spatial_y", "Time_short", "Outcome_short", "Group"]
    data = pd.read_parquet(B0_TABLE, columns=columns)
    data = data.dropna(subset=["sample", "panel", "celltype", "spatial_x", "spatial_y"]).copy()
    sample_order = meta["sample"].tolist()
    if args.samples:
        requested = set(args.samples)
        sample_order = [s for s in sample_order if s in requested]
        data = data[data["sample"].isin(sample_order)].copy()

    proximity_parts = []
    coherence_parts = []
    for i, sample in enumerate(sample_order, start=1):
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample)
        prox_path = intermediate_dir / f"{safe}.proximity.csv"
        coh_path = intermediate_dir / f"{safe}.coherence.csv"
        if args.resume and prox_path.exists() and coh_path.exists():
            log(f"[{i}/{len(sample_order)}] reuse {sample}")
            proximity_parts.append(pd.read_csv(prox_path))
            coherence_parts.append(pd.read_csv(coh_path))
            continue
        log(f"[{i}/{len(sample_order)}] process {sample}")
        sdf = data[data["sample"].eq(sample)].copy()
        prox, coh = process_one_sample(
            sample,
            sdf,
            n_permutations=args.n_permutations,
            random_seed=args.random_seed,
            perm_source_cap=args.perm_source_cap,
        )
        prox.to_csv(prox_path, index=False)
        coh.to_csv(coh_path, index=False)
        proximity_parts.append(prox)
        coherence_parts.append(coh)

    proximity = pd.concat(proximity_parts, ignore_index=True)
    coherence = pd.concat(coherence_parts, ignore_index=True)
    proximity = proximity.merge(meta, on="sample", how="left", validate="many_to_one")
    coherence = coherence.merge(meta, on="sample", how="left", validate="many_to_one")

    features = pd.concat(
        [proximity_feature_table(proximity, meta), coherence_feature_table(coherence, meta)],
        ignore_index=True,
    )
    stats_df = compute_crnr_stats(features, bootstrap_n=args.bootstrap_n, random_seed=args.random_seed)
    report = build_report(proximity, coherence, stats_df)

    proximity.to_csv(OUT_PROXIMITY, index=False)
    coherence.to_csv(OUT_COHERENCE, index=False)
    stats_df.to_csv(OUT_STATS, index=False)
    OUT_REPORT.write_text(report, encoding="utf-8")
    command = (
        "python deliverables/WP1.5_proximity_coherence.py "
        f"--n-permutations {args.n_permutations} --bootstrap-n {args.bootstrap_n} "
        f"--random-seed {args.random_seed} --perm-source-cap {args.perm_source_cap}"
    )
    if args.resume:
        command += " --resume"
    for path in [OUT_PROXIMITY, OUT_COHERENCE, OUT_STATS, OUT_REPORT]:
        write_metadata(path, B0_TABLE, command)
    log(f"Wrote {OUT_PROXIMITY}")
    log(f"Wrote {OUT_COHERENCE}")
    log(f"Wrote {OUT_STATS}")
    log(f"Wrote {OUT_REPORT}")


if __name__ == "__main__":
    main()
