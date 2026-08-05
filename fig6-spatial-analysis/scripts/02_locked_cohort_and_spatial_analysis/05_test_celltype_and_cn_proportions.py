#!/usr/bin/env python3
"""WP1.3 baseline CR-vs-NR sample proportion association.

This script does not recluster cells. It reads WP1.1 CN assignments and WP1.2
functional niche labels, aggregates valid CN-participant cells to sample-level
proportions, and tests baseline Pre CR-vs-NR associations.
"""

from __future__ import annotations

import json
import math
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats


PROJECT_DIR = Path("/data/home/chenwenhui/mIHC/fig6_pm")
CELLINFO_DIR = Path("/data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh")
CN_DIR = CELLINFO_DIR / "neighborhood_radius200_B0"
DELIVERABLE_DIR = PROJECT_DIR / "deliverables"

CN_TABLE = CN_DIR / "neighborhood_composition.parquet"
SAMPLE_VALIDATION = DELIVERABLE_DIR / "WP0.3_sample_validation.csv"
CN_QC = DELIVERABLE_DIR / "WP1.1_cn_qc_by_sample.csv"
NICHE_MAP_K10 = DELIVERABLE_DIR / "WP1.2_functional_niche_map_k10.csv"
NICHE_MAP_K11 = DELIVERABLE_DIR / "WP1.2_functional_niche_map_k11.csv"

OUT_LONG = DELIVERABLE_DIR / "WP1.3_sample_proportion_long.csv"
OUT_STATS = DELIVERABLE_DIR / "WP1.3_CRvsNR_stats.csv"
OUT_QC_PANEL = DELIVERABLE_DIR / "WP1.3_qc_panel_confound.csv"
OUT_CD21 = DELIVERABLE_DIR / "WP1.3_cd21_zero_inflation.csv"
OUT_REPORT = DELIVERABLE_DIR / "WP1.3_association_report.md"

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
PANEL_CELLTYPE = "Fibroblast / SMA+ stromal"
PRIMARY_LOW_CELL_THRESHOLD = 30000
SENSITIVITY_LOW_CELL_THRESHOLD = 20000


@dataclass(frozen=True)
class AnalysisConfig:
    name: str
    include_celltype: bool
    cn_k: str | None
    low_cell_threshold: int | None
    exclude_v1_all_features: bool = False


ANALYSIS_CONFIGS = [
    AnalysisConfig(
        name="primary_k10_threshold30k_panelaware",
        include_celltype=True,
        cn_k="10",
        low_cell_threshold=PRIMARY_LOW_CELL_THRESHOLD,
    ),
    AnalysisConfig(
        name="sensitivity_k10_include_lowcells_panelaware",
        include_celltype=True,
        cn_k="10",
        low_cell_threshold=None,
    ),
    AnalysisConfig(
        name="sensitivity_k10_threshold20k_panelaware",
        include_celltype=True,
        cn_k="10",
        low_cell_threshold=SENSITIVITY_LOW_CELL_THRESHOLD,
    ),
    AnalysisConfig(
        name="sensitivity_k10_threshold30k_exclude_v1_all",
        include_celltype=True,
        cn_k="10",
        low_cell_threshold=PRIMARY_LOW_CELL_THRESHOLD,
        exclude_v1_all_features=True,
    ),
    AnalysisConfig(
        name="sensitivity_k11_threshold30k_panelaware",
        include_celltype=False,
        cn_k="11",
        low_cell_threshold=PRIMARY_LOW_CELL_THRESHOLD,
    ),
]


def normalize_cn_label(value) -> str:
    """Return labels as CN0/CN1/... regardless of string or numeric input."""
    if pd.isna(value):
        return ""
    text = str(value)
    if text.startswith("CN"):
        return text
    try:
        return f"CN{int(float(text))}"
    except ValueError:
        return text


def read_sample_meta() -> pd.DataFrame:
    validation = pd.read_csv(SAMPLE_VALIDATION)
    qc = pd.read_csv(CN_QC)
    qc_cols = [
        "sample",
        "cn_participant_cells",
        "cn_filtered_cells",
        "cn_participant_fraction",
        "low_cell_lt30000",
        "neighbor_count_median",
        "neighbor_count_mean",
    ]
    meta = validation.merge(qc[qc_cols], on="sample", how="left", validate="one_to_one")
    meta["panel_layout"] = np.where(meta["sample"].isin(V1_LIKE_SAMPLES), V1_LAYOUT, V2_LAYOUT)
    meta["is_v1_like"] = meta["sample"].isin(V1_LIKE_SAMPLES)
    meta["is_baseline_pre"] = meta["Time_short"].eq("Pre")
    meta["passes_lowcell_30000"] = meta["cn_participant_cells"].ge(PRIMARY_LOW_CELL_THRESHOLD)
    meta["passes_lowcell_20000"] = meta["cn_participant_cells"].ge(SENSITIVITY_LOW_CELL_THRESHOLD)
    return meta


def read_cn_labels() -> pd.DataFrame:
    columns = ["sample", "celltype", "valid_neighborhood", "CN_k10", "CN_k11"]
    labels = pd.read_parquet(CN_TABLE, columns=columns)
    labels = labels.loc[labels["valid_neighborhood"].astype(bool), columns].copy()
    for col in ["sample", "celltype", "CN_k10", "CN_k11"]:
        labels[col] = labels[col].astype("category")
    return labels


def prepare_niche_map(path: Path) -> pd.DataFrame:
    niche_map = pd.read_csv(path)
    niche_map = niche_map.copy()
    niche_map["CN"] = niche_map["CN"].map(normalize_cn_label)
    return niche_map


def _complete_feature_counts(
    counts: pd.DataFrame,
    samples: list[str],
    features: list[str],
    denominator: pd.Series,
    feature_domain: str,
    k: str,
    feature_name_map: dict[str, str],
    family_map: dict[str, str],
    v1_like_samples: set[str],
) -> pd.DataFrame:
    index = pd.MultiIndex.from_product([samples, features], names=["sample", "feature_id"])
    out = counts.set_index(["sample", "feature_id"]).reindex(index, fill_value=0).reset_index()
    out["feature_domain"] = feature_domain
    out["k"] = k
    out["feature_name"] = out["feature_id"].map(feature_name_map).fillna(out["feature_id"])
    out["niche_family"] = out["feature_id"].map(family_map).fillna("")
    out["count"] = out["count"].astype("int64")
    out["denominator"] = out["sample"].map(denominator).astype("int64")
    out["proportion"] = np.where(out["denominator"] > 0, out["count"] / out["denominator"], np.nan)
    out["panel_confounded_direct"] = is_direct_panel_confounded(
        out["feature_domain"], out["feature_name"], out["niche_family"]
    )
    out["is_v1_like_sample"] = out["sample"].isin(v1_like_samples)
    return out


def is_direct_panel_confounded(
    feature_domain: pd.Series | str,
    feature_name: pd.Series | str,
    niche_family: pd.Series | str,
) -> pd.Series | bool:
    if isinstance(feature_domain, pd.Series):
        return (feature_domain.eq("celltype") & feature_name.eq(PANEL_CELLTYPE)) | (
            feature_domain.eq("CN") & niche_family.eq("fibroblast_stromal")
        )
    return (feature_domain == "celltype" and feature_name == PANEL_CELLTYPE) or (
        feature_domain == "CN" and niche_family == "fibroblast_stromal"
    )


def build_sample_proportion_long(
    labels: pd.DataFrame,
    sample_meta: pd.DataFrame,
    niche_map_k10: pd.DataFrame,
    niche_map_k11: pd.DataFrame,
    v1_like_samples: set[str],
) -> pd.DataFrame:
    sample_meta = sample_meta.copy()
    if "panel_layout" not in sample_meta.columns:
        sample_meta["panel_layout"] = np.where(sample_meta["sample"].isin(v1_like_samples), V1_LAYOUT, V2_LAYOUT)
    if "is_v1_like" not in sample_meta.columns:
        sample_meta["is_v1_like"] = sample_meta["sample"].isin(v1_like_samples)
    if "passes_lowcell_30000" not in sample_meta.columns:
        sample_meta["passes_lowcell_30000"] = sample_meta["cn_participant_cells"].ge(PRIMARY_LOW_CELL_THRESHOLD)
    if "passes_lowcell_20000" not in sample_meta.columns:
        sample_meta["passes_lowcell_20000"] = sample_meta["cn_participant_cells"].ge(SENSITIVITY_LOW_CELL_THRESHOLD)
    labels = labels.loc[labels["valid_neighborhood"].astype(bool)].copy()
    sample_order = sample_meta["sample"].tolist()
    denominator = labels.groupby("sample", observed=True).size().astype("int64")
    denominator.index = denominator.index.astype(str)

    celltypes = sorted(labels["celltype"].astype(str).dropna().unique())
    cell_counts = (
        labels.groupby(["sample", "celltype"], observed=True)
        .size()
        .rename("count")
        .reset_index()
    )
    cell_counts["sample"] = cell_counts["sample"].astype(str)
    cell_counts["feature_id"] = cell_counts["celltype"].astype(str)
    cell_counts = cell_counts.groupby(["sample", "feature_id"], as_index=False, observed=True)["count"].sum()
    cell_long = _complete_feature_counts(
        counts=cell_counts,
        samples=sample_order,
        features=celltypes,
        denominator=denominator,
        feature_domain="celltype",
        k="celltype",
        feature_name_map={feature: feature for feature in celltypes},
        family_map={feature: "" for feature in celltypes},
        v1_like_samples=v1_like_samples,
    )

    cn_longs = []
    for k, col, niche_map in [("10", "CN_k10", niche_map_k10), ("11", "CN_k11", niche_map_k11)]:
        features = niche_map["CN"].map(normalize_cn_label).tolist()
        counts = labels.groupby(["sample", col], observed=True).size().rename("count").reset_index()
        counts["sample"] = counts["sample"].astype(str)
        counts["feature_id"] = counts[col].map(normalize_cn_label)
        counts = counts.groupby(["sample", "feature_id"], as_index=False, observed=True)["count"].sum()
        name_map = dict(zip(niche_map["CN"], niche_map["functional_niche"]))
        family_map = dict(zip(niche_map["CN"], niche_map["niche_family"]))
        cn_longs.append(
            _complete_feature_counts(
                counts=counts,
                samples=sample_order,
                features=features,
                denominator=denominator,
                feature_domain="CN",
                k=k,
                feature_name_map=name_map,
                family_map=family_map,
                v1_like_samples=v1_like_samples,
            )
        )

    long_df = pd.concat([cell_long, *cn_longs], ignore_index=True)
    long_df = long_df.merge(sample_meta, on="sample", how="left", validate="many_to_one")
    long_df["denominator_matches_qc"] = long_df["denominator"].eq(long_df["cn_participant_cells"])

    columns = [
        "sample",
        "Time_short",
        "Outcome_short",
        "Group",
        "panel_layout",
        "is_v1_like",
        "b0_defined_rows",
        "defined_fraction",
        "cn_participant_cells",
        "cn_filtered_cells",
        "cn_participant_fraction",
        "low_cell_lt30000",
        "passes_lowcell_30000",
        "passes_lowcell_20000",
        "denominator",
        "denominator_matches_qc",
        "feature_domain",
        "k",
        "feature_id",
        "feature_name",
        "niche_family",
        "count",
        "proportion",
        "panel_confounded_direct",
    ]
    return long_df[columns].sort_values(["sample", "feature_domain", "k", "feature_id"]).reset_index(drop=True)


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
    denom_df = len(cr) + len(nr) - 2
    pooled_var = ((len(cr) - 1) * var_cr + (len(nr) - 1) * var_nr) / denom_df
    pooled_sd = math.sqrt(pooled_var) if pooled_var > 0 else 0.0
    if pooled_sd == 0:
        return 0.0 if diff == 0 else np.nan
    return diff / pooled_sd


def bootstrap_cohens_d_ci(
    cr_values: Iterable[float],
    nr_values: Iterable[float],
    *,
    n_boot: int = 5000,
    seed: int = 0,
) -> tuple[float, float]:
    cr = np.asarray(list(cr_values), dtype=float)
    nr = np.asarray(list(nr_values), dtype=float)
    cr = cr[np.isfinite(cr)]
    nr = nr[np.isfinite(nr)]
    if len(cr) == 0 or len(nr) == 0:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        cr_b = rng.choice(cr, size=len(cr), replace=True)
        nr_b = rng.choice(nr, size=len(nr), replace=True)
        boot[i] = cohens_d_nr_minus_cr(cr_b, nr_b)
    boot = boot[np.isfinite(boot)]
    if len(boot) == 0:
        return (np.nan, np.nan)
    return (float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975)))


def mannwhitney_p(cr_values: Iterable[float], nr_values: Iterable[float]) -> float:
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
    mask = np.isfinite(p)
    if not np.any(mask):
        return q
    finite_p = p[mask]
    order = np.argsort(finite_p)
    ranked = finite_p[order]
    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty_like(finite_p)
    out[order] = adjusted
    q[mask] = out
    return q


def _feature_filter_for_config(long_df: pd.DataFrame, config: AnalysisConfig) -> pd.Series:
    mask = long_df["Time_short"].eq("Pre") & long_df["Outcome_short"].isin(["CR", "NR"])
    if config.low_cell_threshold is not None:
        mask &= long_df["cn_participant_cells"].ge(config.low_cell_threshold)

    feature_mask = pd.Series(False, index=long_df.index)
    if config.include_celltype:
        feature_mask |= long_df["feature_domain"].eq("celltype")
    if config.cn_k is not None:
        feature_mask |= long_df["feature_domain"].eq("CN") & long_df["k"].eq(config.cn_k)
    return mask & feature_mask


def compute_crnr_stats(
    long_df: pd.DataFrame,
    *,
    bootstrap_n: int = 5000,
    random_seed: int = 0,
) -> pd.DataFrame:
    records = []
    group_cols = [
        "feature_domain",
        "k",
        "feature_id",
        "feature_name",
        "niche_family",
        "panel_confounded_direct",
    ]
    for config in ANALYSIS_CONFIGS:
        analysis_df = long_df.loc[_feature_filter_for_config(long_df, config)].copy()
        for key, sub in analysis_df.groupby(group_cols, dropna=False, sort=False):
            (
                feature_domain,
                k,
                feature_id,
                feature_name,
                niche_family,
                panel_confounded_direct,
            ) = key
            feature_df = sub
            if config.exclude_v1_all_features:
                feature_df = feature_df.loc[feature_df["panel_layout"].ne(V1_LAYOUT)]
                panel_handling = "excluded_v1_all_features"
            elif bool(panel_confounded_direct):
                feature_df = feature_df.loc[feature_df["panel_layout"].ne(V1_LAYOUT)]
                panel_handling = "excluded_v1_for_direct_panel_confounded_feature"
            else:
                panel_handling = "included_all_panel_layouts"

            cr = feature_df.loc[feature_df["Outcome_short"].eq("CR"), "proportion"].astype(float).to_numpy()
            nr = feature_df.loc[feature_df["Outcome_short"].eq("NR"), "proportion"].astype(float).to_numpy()
            d = cohens_d_nr_minus_cr(cr, nr)
            seed = (
                random_seed
                + zlib.crc32(f"{config.name}|{feature_domain}|{k}|{feature_id}".encode("utf-8"))
            ) & 0xFFFFFFFF
            ci_low, ci_high = bootstrap_cohens_d_ci(cr, nr, n_boot=bootstrap_n, seed=seed)
            p_value = mannwhitney_p(cr, nr)
            if np.isfinite(d):
                direction = "higher_in_NR" if d > 0 else "higher_in_CR" if d < 0 else "no_difference"
            else:
                diff = np.nanmean(nr) - np.nanmean(cr) if len(nr) and len(cr) else np.nan
                direction = "higher_in_NR" if diff > 0 else "higher_in_CR" if diff < 0 else "no_difference"
            records.append(
                {
                    "analysis_set": config.name,
                    "low_cell_threshold": config.low_cell_threshold if config.low_cell_threshold is not None else "none",
                    "exclude_v1_all_features": config.exclude_v1_all_features,
                    "panel_handling": panel_handling,
                    "feature_domain": feature_domain,
                    "k": k,
                    "feature_id": feature_id,
                    "feature_name": feature_name,
                    "niche_family": niche_family,
                    "panel_confounded_direct": bool(panel_confounded_direct),
                    "n_CR": int(len(cr)),
                    "n_NR": int(len(nr)),
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
        stats_df.loc[idx, "FDR_BH_within_analysis_set"] = bh_fdr(
            stats_df.loc[idx, "wilcoxon_mannwhitney_p"].to_numpy()
        )
    stats_df["passes_FDR_0.10"] = stats_df["FDR_BH_within_analysis_set"].lt(0.10)
    return stats_df.sort_values(
        ["analysis_set", "FDR_BH_within_analysis_set", "feature_domain", "feature_id"],
        na_position="last",
    ).reset_index(drop=True)


def _fisher_test_record(meta: pd.DataFrame, scope: str, variable: str, levels: list[str]) -> dict:
    rows = [V1_LAYOUT, V2_LAYOUT]
    table = []
    for layout in rows:
        table.append([int(((meta["panel_layout"] == layout) & (meta[variable] == level)).sum()) for level in levels])
    if len(levels) == 2:
        oddsratio, p_value = scipy_stats.fisher_exact(table)
    else:
        oddsratio, p_value = (np.nan, np.nan)
    return {
        "record_type": "test",
        "scope": scope,
        "panel_layout": "",
        "Time_short": "",
        "Outcome_short": "",
        "n_samples": int(len(meta)),
        "test": f"Fisher exact panel_layout x {variable}",
        "variable": variable,
        "levels": "|".join(levels),
        "contingency_table_json": json.dumps({"rows": rows, "columns": levels, "counts": table}),
        "oddsratio": oddsratio,
        "p_value": p_value,
    }


def make_panel_qc(sample_meta: pd.DataFrame) -> pd.DataFrame:
    records = []
    counts = (
        sample_meta.groupby(["panel_layout", "Time_short", "Outcome_short"], dropna=False)
        .size()
        .rename("n_samples")
        .reset_index()
    )
    for _, row in counts.iterrows():
        records.append(
            {
                "record_type": "count",
                "scope": "all_samples_panel_by_time_outcome",
                "panel_layout": row["panel_layout"],
                "Time_short": row["Time_short"],
                "Outcome_short": row["Outcome_short"],
                "n_samples": int(row["n_samples"]),
                "test": "",
                "variable": "",
                "levels": "",
                "contingency_table_json": "",
                "oddsratio": np.nan,
                "p_value": np.nan,
            }
        )

    records.append(_fisher_test_record(sample_meta, "all_samples", "Outcome_short", ["CR", "NR"]))
    records.append(_fisher_test_record(sample_meta, "all_samples", "Time_short", ["Pre", "Post"]))
    pre = sample_meta.loc[sample_meta["Time_short"].eq("Pre")]
    post = sample_meta.loc[sample_meta["Time_short"].eq("Post")]
    records.append(_fisher_test_record(pre, "baseline_pre_only", "Outcome_short", ["CR", "NR"]))
    records.append(_fisher_test_record(post, "post_only", "Outcome_short", ["CR", "NR"]))
    cr = sample_meta.loc[sample_meta["Outcome_short"].eq("CR")]
    nr = sample_meta.loc[sample_meta["Outcome_short"].eq("NR")]
    records.append(_fisher_test_record(cr, "CR_only", "Time_short", ["Pre", "Post"]))
    records.append(_fisher_test_record(nr, "NR_only", "Time_short", ["Pre", "Post"]))
    return pd.DataFrame(records)


def make_cd21_zero_inflation(long_df: pd.DataFrame) -> pd.DataFrame:
    target_mask = (
        (long_df["feature_domain"].eq("celltype") & long_df["feature_id"].eq("CD21+ B cell"))
        | (long_df["feature_domain"].eq("CN") & long_df["niche_family"].eq("B_cell_lymphoid"))
    )
    base = long_df.loc[
        target_mask
        & long_df["Time_short"].eq("Pre")
        & long_df["Outcome_short"].isin(["CR", "NR"])
        & long_df["cn_participant_cells"].ge(PRIMARY_LOW_CELL_THRESHOLD)
    ].copy()
    records = []
    thresholds = [("gt0", 0.0, "proportion > 0"), ("ge1pct", 0.01, "proportion >= 0.01")]
    group_cols = ["feature_domain", "k", "feature_id", "feature_name", "niche_family"]
    for key, sub in base.groupby(group_cols, sort=False, dropna=False):
        feature_domain, k, feature_id, feature_name, niche_family = key
        for threshold_name, threshold, threshold_label in thresholds:
            if threshold_name == "gt0":
                present = sub["proportion"].gt(0)
            else:
                present = sub["proportion"].ge(threshold)
            sub = sub.assign(_present=present.to_numpy())
            nr_present = int(((sub["Outcome_short"] == "NR") & sub["_present"]).sum())
            nr_absent = int(((sub["Outcome_short"] == "NR") & ~sub["_present"]).sum())
            cr_present = int(((sub["Outcome_short"] == "CR") & sub["_present"]).sum())
            cr_absent = int(((sub["Outcome_short"] == "CR") & ~sub["_present"]).sum())
            oddsratio, p_value = scipy_stats.fisher_exact([[nr_present, nr_absent], [cr_present, cr_absent]])
            cr_pos = sub.loc[sub["Outcome_short"].eq("CR") & sub["_present"], "proportion"].astype(float).to_numpy()
            nr_pos = sub.loc[sub["Outcome_short"].eq("NR") & sub["_present"], "proportion"].astype(float).to_numpy()
            records.append(
                {
                    "feature_domain": feature_domain,
                    "k": k,
                    "feature_id": feature_id,
                    "feature_name": feature_name,
                    "niche_family": niche_family,
                    "presence_definition": threshold_label,
                    "presence_threshold_name": threshold_name,
                    "n_CR": int((sub["Outcome_short"] == "CR").sum()),
                    "n_NR": int((sub["Outcome_short"] == "NR").sum()),
                    "CR_present": cr_present,
                    "CR_absent": cr_absent,
                    "NR_present": nr_present,
                    "NR_absent": nr_absent,
                    "presence_oddsratio_NR_vs_CR": oddsratio,
                    "presence_fisher_p": p_value,
                    "CR_positive_n": int(len(cr_pos)),
                    "NR_positive_n": int(len(nr_pos)),
                    "positive_mean_CR": float(np.nanmean(cr_pos)) if len(cr_pos) else np.nan,
                    "positive_mean_NR": float(np.nanmean(nr_pos)) if len(nr_pos) else np.nan,
                    "positive_median_CR": float(np.nanmedian(cr_pos)) if len(cr_pos) else np.nan,
                    "positive_median_NR": float(np.nanmedian(nr_pos)) if len(nr_pos) else np.nan,
                    "positive_cohens_d_NR_minus_CR": cohens_d_nr_minus_cr(cr_pos, nr_pos),
                    "positive_mannwhitney_p": mannwhitney_p(cr_pos, nr_pos),
                }
            )
    return pd.DataFrame(records)


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


def _markdown_table(df: pd.DataFrame, columns: list[str], max_rows: int = 10) -> str:
    if df.empty:
        return "None."
    view = df.loc[:, columns].head(max_rows).copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            if "p" in col.lower() or "fdr" in col.lower():
                view[col] = view[col].map(fmt_p)
            else:
                view[col] = view[col].map(lambda x: fmt_float(x, 3))
    view = view.fillna("NA").astype(str)
    header = "| " + " | ".join(view.columns) + " |"
    separator = "| " + " | ".join(["---"] * len(view.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in view.to_numpy()]
    return "\n".join([header, separator, *rows])


def build_report(
    sample_meta: pd.DataFrame,
    stats_df: pd.DataFrame,
    panel_qc: pd.DataFrame,
    cd21_df: pd.DataFrame,
) -> str:
    baseline = sample_meta.loc[sample_meta["Time_short"].eq("Pre") & sample_meta["Outcome_short"].isin(["CR", "NR"])]
    excluded30 = baseline.loc[~baseline["cn_participant_cells"].ge(PRIMARY_LOW_CELL_THRESHOLD), "sample"].tolist()
    excluded20 = baseline.loc[~baseline["cn_participant_cells"].ge(SENSITIVITY_LOW_CELL_THRESHOLD), "sample"].tolist()
    counts_by_outcome = baseline.groupby("Outcome_short").size().to_dict()
    counts30 = baseline.loc[baseline["cn_participant_cells"].ge(PRIMARY_LOW_CELL_THRESHOLD)].groupby("Outcome_short").size().to_dict()

    primary = stats_df.loc[stats_df["analysis_set"].eq("primary_k10_threshold30k_panelaware")].copy()
    primary["abs_d"] = primary["cohens_d_NR_minus_CR"].abs()
    top_by_abs_d = primary.sort_values("abs_d", ascending=False)
    sig = primary.loc[primary["FDR_BH_within_analysis_set"].lt(0.10)].sort_values("FDR_BH_within_analysis_set")

    panel_tests = panel_qc.loc[panel_qc["record_type"].eq("test")].copy()
    baseline_panel = panel_tests.loc[panel_tests["scope"].eq("baseline_pre_only")]
    baseline_panel_p = baseline_panel["p_value"].iloc[0] if not baseline_panel.empty else np.nan
    all_panel_outcome = panel_tests.loc[panel_tests["scope"].eq("all_samples") & panel_tests["variable"].eq("Outcome_short")]
    all_panel_outcome_p = all_panel_outcome["p_value"].iloc[0] if not all_panel_outcome.empty else np.nan

    cd21_main = cd21_df.loc[
        (cd21_df["feature_domain"].eq("CN"))
        & (cd21_df["k"].eq("10"))
        & (cd21_df["presence_threshold_name"].eq("ge1pct"))
    ]

    lines = [
        "# WP1.3 Baseline CR-vs-NR Proportion Association",
        "",
        "## Scope and clustering provenance",
        "",
        "WP1.3 did not rerun clustering. It directly read the WP1.1 B0 CN assignment table "
        f"(`{CN_TABLE}`) and the WP1.2 functional niche maps. The clustering provenance is: "
        "WP0.3/B0 cleaned cell phenotypes -> WP1.1 radius-200 neighborhood cell-type "
        "composition vectors for valid cells -> KMeans for k=10 (primary) and k=11 "
        "(sensitivity, random_state=0, n_init=50, no extra scaling of composition vectors) -> "
        "WP1.2 functional names from cluster centers/enrichments -> WP1.3 sample-level "
        "proportion aggregation and CR-vs-NR tests.",
        "",
        "## Cohort and QC gates",
        "",
        f"Baseline Pre CR/NR samples before the low-cell gate: CR={counts_by_outcome.get('CR', 0)}, "
        f"NR={counts_by_outcome.get('NR', 0)}.",
        f"Primary low-cell gate: cn_participant_cells >= {PRIMARY_LOW_CELL_THRESHOLD:,}. "
        f"After the gate: CR={counts30.get('CR', 0)}, NR={counts30.get('NR', 0)}. "
        f"Excluded baseline samples: {', '.join(excluded30) if excluded30 else 'none'}.",
        f"Low-cell sensitivity also used threshold {SENSITIVITY_LOW_CELL_THRESHOLD:,}; excluded: "
        f"{', '.join(excluded20) if excluded20 else 'none'}.",
        f"Panel-layout QC: all-sample panel x outcome Fisher p={fmt_p(all_panel_outcome_p)}; "
        f"baseline-only panel x outcome Fisher p={fmt_p(baseline_panel_p)}. Direct panel-dependent "
        "features (Fibroblast / SMA+ stromal celltype and fibroblast/stromal CNs) exclude the "
        "v1-like samples in the primary panel-aware analysis. A stricter sensitivity excludes "
        "v1-like samples for every feature.",
        "",
        "## Primary k=10 results",
        "",
    ]

    if sig.empty:
        lines.append("No primary k=10 celltype/CN feature passed BH-FDR < 0.10.")
    else:
        lines.append(f"{len(sig)} primary features passed BH-FDR < 0.10:")
        lines.append(
            _markdown_table(
                sig,
                [
                    "feature_domain",
                    "k",
                    "feature_id",
                    "feature_name",
                    "cohens_d_NR_minus_CR",
                    "FDR_BH_within_analysis_set",
                    "direction",
                ],
                max_rows=20,
            )
        )

    lines.extend(
        [
            "",
            "Top primary features by absolute Cohen's d (NR minus CR):",
            "",
            _markdown_table(
                top_by_abs_d,
                [
                    "feature_domain",
                    "k",
                    "feature_id",
                    "feature_name",
                    "n_CR",
                    "n_NR",
                    "mean_CR",
                    "mean_NR",
                    "cohens_d_NR_minus_CR",
                    "cohens_d_bootstrap95_low",
                    "cohens_d_bootstrap95_high",
                    "wilcoxon_mannwhitney_p",
                    "FDR_BH_within_analysis_set",
                    "direction",
                    "panel_handling",
                ],
                max_rows=12,
            ),
            "",
            "## CD21+ B-cell lymphoid zero-inflation",
            "",
        ]
    )

    if cd21_main.empty:
        lines.append("CD21+ B-cell lymphoid niche row was not found.")
    else:
        row = cd21_main.iloc[0]
        lines.append(
            "Primary k=10 CD21+ B-cell lymphoid niche presence at >=1%: "
            f"CR {int(row['CR_present'])}/{int(row['n_CR'])}, "
            f"NR {int(row['NR_present'])}/{int(row['n_NR'])}; "
            f"Fisher p={fmt_p(row['presence_fisher_p'])}, "
            f"oddsratio_NR_vs_CR={fmt_float(row['presence_oddsratio_NR_vs_CR'])}. "
            f"Among positive samples, magnitude Mann-Whitney p={fmt_p(row['positive_mannwhitney_p'])}, "
            f"Cohen d={fmt_float(row['positive_cohens_d_NR_minus_CR'])}."
        )
    lines.extend(
        [
            "",
            "Full CD21 two-stage results are in `WP1.3_cd21_zero_inflation.csv`.",
            "",
            "## Sensitivity analyses",
            "",
        ]
    )
    sensitivity_summary = (
        stats_df.groupby("analysis_set")
        .agg(
            n_features=("feature_id", "size"),
            n_FDR_lt_0_10=("passes_FDR_0.10", "sum"),
            min_FDR=("FDR_BH_within_analysis_set", "min"),
        )
        .reset_index()
        .sort_values("analysis_set")
    )
    lines.append(_markdown_table(sensitivity_summary, ["analysis_set", "n_features", "n_FDR_lt_0_10", "min_FDR"], 20))
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- `{OUT_LONG}`",
            f"- `{OUT_STATS}`",
            f"- `{OUT_QC_PANEL}`",
            f"- `{OUT_CD21}`",
            f"- `{OUT_REPORT}`",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    DELIVERABLE_DIR.mkdir(parents=True, exist_ok=True)
    sample_meta = read_sample_meta()
    niche_map_k10 = prepare_niche_map(NICHE_MAP_K10)
    niche_map_k11 = prepare_niche_map(NICHE_MAP_K11)
    labels = read_cn_labels()

    long_df = build_sample_proportion_long(labels, sample_meta, niche_map_k10, niche_map_k11, V1_LIKE_SAMPLES)
    long_df.to_csv(OUT_LONG, index=False)

    stats_df = compute_crnr_stats(long_df, bootstrap_n=5000, random_seed=13)
    stats_df.to_csv(OUT_STATS, index=False)

    panel_qc = make_panel_qc(sample_meta)
    panel_qc.to_csv(OUT_QC_PANEL, index=False)

    cd21_df = make_cd21_zero_inflation(long_df)
    cd21_df.to_csv(OUT_CD21, index=False)

    report = build_report(sample_meta, stats_df, panel_qc, cd21_df)
    OUT_REPORT.write_text(report, encoding="utf-8")

    print(f"Wrote {OUT_LONG}")
    print(f"Wrote {OUT_STATS}")
    print(f"Wrote {OUT_QC_PANEL}")
    print(f"Wrote {OUT_CD21}")
    print(f"Wrote {OUT_REPORT}")


if __name__ == "__main__":
    main()
