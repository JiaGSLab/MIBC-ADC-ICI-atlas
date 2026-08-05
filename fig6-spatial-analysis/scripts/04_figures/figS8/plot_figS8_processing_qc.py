#!/usr/bin/env python3
"""WP3.3 supplementary Fig. 6 candidate panels.

Builds S8-S12 from existing deliverables only. No upstream analysis is rerun.
Python/matplotlib is the sole plotting backend for all exports.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle
from PIL import Image


PROJECT_DIR = Path("/data/home/chenwenhui/mIHC/fig6_pm")
DELIVERABLE_DIR = PROJECT_DIR / "deliverables"
CELLINFO_DIR = Path("/data/home/chenwenhui/mIHC/CellInfo")
B0_DIR = CELLINFO_DIR / "BLCA_fig6_B0_refresh"
PANEL_DIR = B0_DIR / "fig6_panels_supp"
SCRATCH_DIR = PROJECT_DIR / "scratch_wp3.3"
SAMPLE_INFO = CELLINFO_DIR / "sample_info _IT.csv"

OUT_SCRIPT = DELIVERABLE_DIR / "WP3_3_supplementary_figures.py"
OUT_INDEX = DELIVERABLE_DIR / "WP3.3_supp_index.csv"
OUT_REPORT = DELIVERABLE_DIR / "WP3.3_supp_report.md"
OUT_METADATA = DELIVERABLE_DIR / "WP3.3_supp_metadata.json"

PANEL_IDS = ["S8", "S9", "S10", "S11", "S12"]
CR_COLOR = "#E4212B"
NR_COLOR = "#2C5F9E"
PRE_COLOR = "#E08A3C"
POST_COLOR = "#4A78B5"
NEUTRAL = "#777777"
PALE = "#EDEDED"
COLORS_15 = [
    "#A52A2A",
    "#053E7A",
    "#E46B4F",
    "#8B8B00",
    "#B8860B",
    "#CD661D",
    "#8B5A2B",
    "#CDAD00",
    "#2F4F4F",
    "#483D8B",
    "#7A378B",
    "#4A7088",
    "#556B2F",
    "#008B8B",
    "#CD7054",
]
GROUP_COLORS = {"CR": CR_COLOR, "NR": NR_COLOR}
TIME_COLORS = {"Pre": PRE_COLOR, "Post": POST_COLOR}

EXCLUDED_ACTIVE_SECTIONS = [
    "pt16-post",
    "pt18-post",
    "pt18-pre",
    "pt33-pre",
    "pt50-post",
    "pt51-post",
    "pt55-post",
    "pt56-post",
    "pt56-pre",
    "pt57-post",
]


@dataclass(frozen=True)
class PanelExport:
    panel_id: str
    title: str
    basename: str
    paths: dict[str, Path]
    source_rows: int
    conclusion: str
    data_sources: str


def required_panel_ids() -> list[str]:
    return PANEL_IDS.copy()


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "DejaVu Sans",
            "font.size": 7,
            "axes.edgecolor": "#333333",
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.dpi": 150,
        }
    )


def fmt(value: float, digits: int = 2) -> str:
    if value is None or not np.isfinite(float(value)):
        return "NA"
    value = float(value)
    if 0 < abs(value) < 0.001:
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def short_label(text: str, max_len: int = 35) -> str:
    text = str(text)
    if len(text) <= max_len:
        return text
    return text[: max_len - 1] + "."


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")


def read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(DELIVERABLE_DIR / name)


def read_json(name: str) -> dict:
    return json.loads((DELIVERABLE_DIR / name).read_text(encoding="utf-8"))


def normalize_sample_info(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["patient_id"] = out["Sample_ID"].astype(str)
    out["Outcome_short"] = out["Outcome"].map({"response": "CR", "non-response": "NR"}).fillna(out["Outcome"].astype(str))
    out["Time_short"] = out["Time_point"].map({"pre-treatment": "Pre", "post-treatment": "Post"}).fillna(out["Time_point"].astype(str))
    out = out.sort_values(["patient_id", "Time_short", "Date", "orig.ident"]).reset_index(drop=True)
    sample_names = []
    for _, group in out.groupby(["patient_id", "Time_short"], sort=False):
        if group["Time_short"].iloc[0] == "Post" and len(group) > 1:
            for rank, (_, row) in enumerate(group.iterrows(), start=1):
                sample_names.append((row.name, f"{row['patient_id']}-post{rank}"))
        else:
            suffix = "pre" if group["Time_short"].iloc[0] == "Pre" else "post"
            for _, row in group.iterrows():
                sample_names.append((row.name, f"{row['patient_id']}-{suffix}"))
    name_lookup = dict(sample_names)
    out["sample"] = [name_lookup[i] for i in out.index]
    return out


def sample_inclusion_summary(sample_info: pd.DataFrame, active: pd.DataFrame, baseline: pd.DataFrame) -> pd.DataFrame:
    sample_info = sample_info.copy()
    active = active.copy()
    baseline = baseline.copy()
    rows = [
        {
            "stage": "experimental_mIHC_sections",
            "n_samples": int(sample_info["sample"].nunique()),
            "n_patients": int(sample_info["patient_id"].nunique()) if "patient_id" in sample_info.columns else np.nan,
            "n_CR": int(sample_info.loc[sample_info["Outcome_short"].eq("CR"), "sample"].nunique()),
            "n_NR": int(sample_info.loc[sample_info["Outcome_short"].eq("NR"), "sample"].nunique()),
            "description": "all IT mIHC sections from sample_info",
        },
        {
            "stage": "downstream_active_sections",
            "n_samples": int(active["sample"].nunique()),
            "n_patients": int(active["sample"].str.extract(r"^(pt\\d+)")[0].nunique()) if "sample" in active.columns else np.nan,
            "n_CR": int(active.loc[active["Outcome_short"].eq("CR"), "sample"].nunique()) if "Outcome_short" in active.columns else np.nan,
            "n_NR": int(active.loc[active["Outcome_short"].eq("NR"), "sample"].nunique()) if "Outcome_short" in active.columns else np.nan,
            "description": "active B0 downstream samples",
        },
        {
            "stage": "baseline_pre_sections",
            "n_samples": int(baseline.loc[baseline["Time_short"].eq("Pre"), "sample"].nunique()),
            "n_patients": int(baseline.loc[baseline["Time_short"].eq("Pre"), "sample"].str.extract(r"^(pt\\d+)")[0].nunique()),
            "n_CR": int(baseline.loc[baseline["Time_short"].eq("Pre") & baseline["Outcome_short"].eq("CR"), "sample"].nunique()),
            "n_NR": int(baseline.loc[baseline["Time_short"].eq("Pre") & baseline["Outcome_short"].eq("NR"), "sample"].nunique()),
            "description": "baseline Pre subset",
        },
        {
            "stage": "baseline_pre_ge30k_sections",
            "n_samples": int(
                baseline.loc[
                    baseline["Time_short"].eq("Pre") & baseline["passes_lowcell_30000"].astype(bool),
                    "sample",
                ].nunique()
            ),
            "n_patients": int(
                baseline.loc[
                    baseline["Time_short"].eq("Pre") & baseline["passes_lowcell_30000"].astype(bool),
                    "sample",
                ]
                .str.extract(r"^(pt\\d+)")[0]
                .nunique()
            ),
            "n_CR": int(
                baseline.loc[
                    baseline["Time_short"].eq("Pre")
                    & baseline["Outcome_short"].eq("CR")
                    & baseline["passes_lowcell_30000"].astype(bool),
                    "sample",
                ].nunique()
            ),
            "n_NR": int(
                baseline.loc[
                    baseline["Time_short"].eq("Pre")
                    & baseline["Outcome_short"].eq("NR")
                    & baseline["passes_lowcell_30000"].astype(bool),
                    "sample",
                ].nunique()
            ),
            "description": "baseline Pre subset used for primary CR/NR statistics",
        },
    ]
    return pd.DataFrame(rows)


def missing_marker_summary(coverage: pd.DataFrame, markers: set[str]) -> pd.DataFrame:
    cols = ["marker", "n_measured_samples", "n_missing_samples", "missing_samples"]
    return coverage.loc[coverage["marker"].isin(markers), cols].copy().sort_values("marker").reset_index(drop=True)


def paths_for_basename(basename: str) -> dict[str, Path]:
    return {
        "svg": PANEL_DIR / f"{basename}.svg",
        "pdf": PANEL_DIR / f"{basename}.pdf",
        "png": PANEL_DIR / f"{basename}.png",
        "tiff": PANEL_DIR / f"{basename}.tiff",
        "source_data": PANEL_DIR / f"{basename}_source_data.csv",
    }


def save_panel(
    fig: plt.Figure,
    *,
    panel_id: str,
    title: str,
    basename: str,
    source: pd.DataFrame,
    conclusion: str,
    data_sources: str,
    dpi: int = 600,
) -> PanelExport:
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    paths = paths_for_basename(basename)
    source.to_csv(paths["source_data"], index=False)
    fig.savefig(paths["svg"], bbox_inches="tight", dpi=dpi)
    fig.savefig(paths["pdf"], bbox_inches="tight", dpi=dpi)
    fig.savefig(paths["png"], bbox_inches="tight", dpi=300)
    fig.savefig(paths["tiff"], bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return PanelExport(panel_id, title, basename, paths, int(len(source)), conclusion, data_sources)


def tagged(df: pd.DataFrame, source_table: str, panel_part: str) -> pd.DataFrame:
    out = df.copy()
    out.insert(0, "panel_part", panel_part)
    out.insert(0, "source_table", source_table)
    return out


def combine_source(frames: Iterable[pd.DataFrame]) -> pd.DataFrame:
    parts = [frame for frame in frames if frame is not None and not frame.empty]
    if not parts:
        return pd.DataFrame({"source_table": [], "panel_part": []})
    return pd.concat(parts, ignore_index=True, sort=False)


def plot_errorbarh(ax: plt.Axes, rows: pd.DataFrame, *, label_col: str, effect_col: str, low_col: str, high_col: str, color_col: str | None = None) -> None:
    rows = rows.copy().reset_index(drop=True)
    y = np.arange(len(rows))
    colors = [rows[color_col].iloc[i] if color_col else NEUTRAL for i in range(len(rows))]
    ax.axvline(0, color="#444444", lw=0.8, ls="--", zorder=0)
    for i, row in rows.iterrows():
        effect = float(row[effect_col])
        lo = float(row[low_col])
        hi = float(row[high_col])
        ax.plot([lo, hi], [i, i], color=colors[i], lw=1.4)
        ax.scatter([effect], [i], color=colors[i], s=28, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([short_label(x, 45) for x in rows[label_col]])
    ax.invert_yaxis()
    ax.set_xlabel("Cohen's d or delta (NR - CR)")


def panel_s8_qc() -> PanelExport:
    sample_info = normalize_sample_info(pd.read_csv(SAMPLE_INFO))
    active = read_csv("WP0.3_sample_validation.csv")
    baseline = read_csv("WP1.3_sample_proportion_long.csv").drop_duplicates("sample")
    numeric = read_json("WP0.1_numeric_summary.json")
    summary = sample_inclusion_summary(sample_info, active, baseline)
    excluded = pd.DataFrame({"sample": EXCLUDED_ACTIVE_SECTIONS, "reason": "excluded before downstream B0/CN analysis"})
    flow = pd.DataFrame(
        {
            "quantity": ["union rows", "B0 defined rows", "undefined removed"],
            "n_cells": [
                int(numeric["parquet"]["num_rows_metadata"]),
                int(active["b0_defined_rows"].sum()),
                int(active["dropped_undefined_rows"].sum()),
            ],
        }
    )

    fig = plt.figure(figsize=(8.2, 5.0))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.25, 1.15], height_ratios=[1, 1], wspace=0.55, hspace=0.55)
    ax = fig.add_subplot(gs[0, 0])
    sns.barplot(data=summary, x="stage", y="n_samples", ax=ax, color=COLORS_15[1])
    ax.set_xticks(np.arange(len(summary)))
    ax.set_xticklabels(["82 all", "72 active", "40 Pre", "36 Pre>=30k"], rotation=25, ha="right")
    ax.set_ylabel("sections")
    ax.set_xlabel("")
    for patch, value in zip(ax.patches, summary["n_samples"]):
        ax.text(patch.get_x() + patch.get_width() / 2, patch.get_height() + 1, str(int(value)), ha="center", va="bottom", fontsize=7)
    ax.set_title("sample inclusion")

    ax = fig.add_subplot(gs[0, 1])
    flow_palette = {"union rows": COLORS_15[2], "B0 defined rows": COLORS_15[1], "undefined removed": "#B8B8B8"}
    sns.barplot(data=flow, x="quantity", y="n_cells", hue="quantity", ax=ax, palette=flow_palette, legend=False)
    ax.set_ylabel("cells / rows")
    ax.set_xlabel("")
    ax.set_xticks(np.arange(len(flow)))
    ax.set_xticklabels(["union\n57.3M", "B0\n23.1M", "undefined\nremoved"], rotation=0)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(6, 6))
    ax.set_title("B0 reconciliation")

    ax = fig.add_subplot(gs[0, 2])
    sns.boxplot(
        data=active,
        x="Outcome_short",
        y="defined_fraction",
        hue="Outcome_short",
        ax=ax,
        palette=GROUP_COLORS,
        width=0.55,
        fliersize=2,
        legend=False,
    )
    sns.stripplot(data=active, x="Outcome_short", y="defined_fraction", ax=ax, color="#333333", size=2.3, alpha=0.55)
    ax.set_xlabel("")
    ax.set_ylabel("defined fraction")
    ax.set_title("QC not outcome-confounded")

    ax = fig.add_subplot(gs[1, :2])
    counts = active.sort_values("b0_defined_rows")
    colors = counts["Time_short"].map(TIME_COLORS).fillna(NEUTRAL)
    ax.bar(np.arange(len(counts)), counts["b0_defined_rows"], color=colors, width=0.8)
    ax.axhline(30000, color="#333333", lw=1.0, ls="--")
    ax.set_ylabel("B0 cells per active section")
    ax.set_xlabel("72 active sections, sorted")
    ax.set_title("low-cell sections flagged, not silently dropped")
    ax.text(0.01, 0.92, "30k primary-stat threshold", transform=ax.transAxes, fontsize=7)

    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    text = "\n".join(excluded["sample"].tolist())
    ax.text(0.0, 1.0, "excluded before downstream\nactive set (10 sections)", fontweight="bold", va="top")
    ax.text(0.0, 0.78, text, va="top", family="monospace", fontsize=7)

    source = combine_source(
        [
            tagged(summary, "sample_info + WP0.3 + WP1.3", "sample_inclusion"),
            tagged(flow, "WP0.1_numeric_summary + WP0.3_sample_validation", "reconciliation_counts"),
            tagged(active, "WP0.3_sample_validation", "sample_qc"),
            tagged(excluded, "manager SAMPLE COUNTS FINAL + active sample list", "excluded_sections"),
        ]
    )
    fig.suptitle("S8. mIHC processing and QC", x=0.01, y=0.995, ha="left", fontsize=9, fontweight="bold")
    return save_panel(
        fig,
        panel_id="S8",
        title="mIHC processing and QC",
        basename="S8_mihc_processing_qc",
        source=source,
        conclusion="The downstream Fig. 6 mIHC analysis uses a reconciled 72-section/23.1M-cell B0 set after explicit undefined-cell removal and QC annotation.",
        data_sources="sample_info _IT.csv; WP0.1_numeric_summary.json; WP0.3_sample_validation.csv; WP1.3_sample_proportion_long.csv",
    )


def panel_s9_celltype_biology() -> PanelExport:
    stats13 = read_csv("WP1.3_CRvsNR_stats.csv")
    cell_primary = stats13[
        (stats13["analysis_set"] == "primary_k10_threshold30k_panelaware")
        & (stats13["feature_domain"] == "celltype")
        & (~stats13["panel_confounded_direct"].astype(bool))
    ].copy()
    cell_primary["abs_d"] = cell_primary["cohens_d_NR_minus_CR"].abs()
    top = cell_primary.sort_values("abs_d", ascending=False).head(10).sort_values("cohens_d_NR_minus_CR")
    top["plot_color"] = np.where(top["cohens_d_NR_minus_CR"] < 0, CR_COLOR, NR_COLOR)

    cd21 = read_csv("WP1.3_cd21_zero_inflation.csv")
    coherence_stats = read_csv("WP1.5_CRvsNR_stats.csv")
    spatial_keys = [
        ("coherence", "same_type_knn_fraction", "CD21+ B cell"),
        ("coherence", "binary_moran_i_knn", "CD21+ B cell"),
        ("proximity", "median_nn_distance", "PanCK+ Tumor -> CD8 T cell"),
        ("proximity", "median_nn_distance", "HER2+ Tumor -> Vascular"),
    ]
    spatial_rows = []
    for feature_class, metric, name in spatial_keys:
        cand = coherence_stats[
            (coherence_stats["analysis_set"] == "primary_threshold30k_panelaware")
            & (coherence_stats["feature_class"] == feature_class)
            & (coherence_stats["metric"] == metric)
            & (coherence_stats["feature_name"] == name)
        ].copy()
        if not cand.empty:
            spatial_rows.append(cand.iloc[0])
    spatial = pd.DataFrame(spatial_rows)
    if not spatial.empty:
        spatial["plot_color"] = np.where(spatial["cohens_d_NR_minus_CR"] < 0, CR_COLOR, NR_COLOR)

    sample_props = read_csv("WP1.3_sample_proportion_long.csv")
    cd21_samples = sample_props[
        (sample_props["analysis_set"].eq("primary_k10_threshold30k_panelaware") if "analysis_set" in sample_props.columns else sample_props["Time_short"].eq("Pre"))
        & (sample_props["Time_short"].eq("Pre"))
        & (sample_props["feature_domain"].eq("celltype"))
        & (sample_props["feature_name"].eq("CD21+ B cell"))
        & (sample_props["passes_lowcell_30000"].astype(bool))
    ].copy()

    fig = plt.figure(figsize=(8.2, 5.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1.0], height_ratios=[1, 0.95], hspace=0.55, wspace=0.45)
    ax = fig.add_subplot(gs[:, 0])
    plot_errorbarh(
        ax,
        top,
        label_col="feature_name",
        effect_col="cohens_d_NR_minus_CR",
        low_col="cohens_d_bootstrap95_low",
        high_col="cohens_d_bootstrap95_high",
        color_col="plot_color",
    )
    ax.set_title("celltype proportion associations")
    ax.set_xlabel("Cohen's d (NR - CR); negative = higher in CR")

    ax = fig.add_subplot(gs[0, 1])
    sns.stripplot(
        data=cd21_samples,
        x="Outcome_short",
        y="proportion",
        hue="Outcome_short",
        ax=ax,
        palette=GROUP_COLORS,
        size=4,
        jitter=0.18,
        legend=False,
    )
    sns.boxplot(data=cd21_samples, x="Outcome_short", y="proportion", ax=ax, showcaps=False, boxprops={"facecolor": "none"}, showfliers=False, whiskerprops={"linewidth": 0})
    ax.set_xlabel("")
    ax.set_ylabel("CD21+ B fraction")
    ax.set_title("CD21+ B/TLS is presence driven")
    row = cd21[(cd21["feature_domain"] == "celltype") & (cd21["presence_threshold_name"] == "ge1pct")]
    if not row.empty:
        r = row.iloc[0]
        ax.text(
            0.02,
            0.95,
            f">=1%: CR {int(r.CR_present)}/{int(r.n_CR)} vs NR {int(r.NR_present)}/{int(r.n_NR)}",
            transform=ax.transAxes,
            va="top",
            fontsize=7,
        )

    ax = fig.add_subplot(gs[1, 1])
    if not spatial.empty:
        spatial_plot = spatial.sort_values("cohens_d_NR_minus_CR")
        plot_errorbarh(
            ax,
            spatial_plot,
            label_col="feature_name",
            effect_col="cohens_d_NR_minus_CR",
            low_col="cohens_d_bootstrap95_low",
            high_col="cohens_d_bootstrap95_high",
            color_col="plot_color",
        )
    ax.set_title("celltype spatial support")
    ax.set_xlabel("Cohen's d (NR - CR)")
    fig.suptitle("S9. Celltype biology remains supportive but candidate-level", x=0.01, y=0.995, ha="left", fontsize=9, fontweight="bold")
    source = combine_source(
        [
            tagged(top, "WP1.3_CRvsNR_stats.csv", "celltype_effects"),
            tagged(cd21, "WP1.3_cd21_zero_inflation.csv", "cd21_presence"),
            tagged(spatial, "WP1.5_CRvsNR_stats.csv", "celltype_spatial"),
            tagged(cd21_samples, "WP1.3_sample_proportion_long.csv", "cd21_sample_values"),
        ]
    )
    return save_panel(
        fig,
        panel_id="S9",
        title="celltype biology",
        basename="S9_celltype_biology",
        source=source,
        conclusion="Celltype-level signals support CD21+ B/TLS biology, but remain candidate-level and do not form a robust standalone predictor.",
        data_sources="WP1.3_CRvsNR_stats.csv; WP1.3_cd21_zero_inflation.csv; WP1.3_sample_proportion_long.csv; WP1.5_CRvsNR_stats.csv",
    )


def load_overlay_thumbnail(path: Path, max_size: tuple[int, int] = (760, 520)) -> np.ndarray:
    image = Image.open(path).convert("RGB")
    image.thumbnail(max_size)
    return np.asarray(image)


def panel_s10_correction_validation() -> PanelExport:
    decision = read_csv("WPB4_marker_decision_table.csv")
    decision = decision.copy()
    decision["plot_color"] = np.where(decision["needs_manual"].astype(bool), "#999999", COLORS_15[1])
    decision = decision.sort_values("mean_delta_crosstalk_minus_residual")
    coverage = read_csv("WPB4_marker_coverage.csv")
    metadata = read_json("WPB4_run_metadata.json")
    overlays = read_csv("WPB4_overlay_index.csv")
    thumb_targets = []
    for marker in ["FOXP3", "CD8", "SMA"]:
        subset = overlays[overlays["marker"].eq(marker)].sort_values("overlay_rank")
        if not subset.empty:
            thumb_targets.append(subset.iloc[min(1, len(subset) - 1)])

    fig = plt.figure(figsize=(8.2, 5.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.08, 0.92], wspace=0.35, hspace=0.42)
    ax = fig.add_subplot(gs[0, :2])
    y = np.arange(len(decision))
    ax.axvline(0, color="#444444", lw=0.8, ls="--")
    for i, row in decision.reset_index(drop=True).iterrows():
        ax.plot([row["delta_ci_low"], row["delta_ci_high"]], [i, i], color=row["plot_color"], lw=1.1)
        ax.scatter(row["mean_delta_crosstalk_minus_residual"], i, color=row["plot_color"], s=20)
    ax.set_yticks(y)
    ax.set_yticklabels(decision["marker"])
    ax.set_xlabel("paired FOV composite delta (crosstalk - residual)")
    ax.set_title("full-cohort image calibration favors global crosstalk")

    ax = fig.add_subplot(gs[0, 2])
    conf = decision["confidence"].value_counts().reindex(["high", "medium", "low", "needs-manual"]).dropna().reset_index()
    conf.columns = ["confidence", "n_markers"]
    sns.barplot(data=conf, x="confidence", y="n_markers", ax=ax, color=COLORS_15[1])
    ax.set_xticks(np.arange(len(conf)))
    ax.set_xticklabels(conf["confidence"], rotation=30, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("markers")
    ax.set_title("decision confidence")
    ax.text(
        0.0,
        0.98,
        f"{metadata['n_requested_samples']}/72 samples\n{metadata['n_fovs']} FOVs\n{metadata['n_markers']} markers",
        transform=ax.transAxes,
        va="top",
        fontsize=7,
    )

    for idx, target in enumerate(thumb_targets):
        ax = fig.add_subplot(gs[1, idx])
        image_path = PROJECT_DIR / str(target["overlay_png"])
        if image_path.exists():
            ax.imshow(load_overlay_thumbnail(image_path))
        ax.set_title(f"{target['marker']} overlay\n{target['sample']} {target['panel']}", fontsize=7)
        ax.set_axis_off()

    fig.suptitle("S10. Correction-layer image validation", x=0.01, y=0.995, ha="left", fontsize=9, fontweight="bold")
    source = combine_source(
        [
            tagged(decision, "WPB4_marker_decision_table.csv", "marker_decisions"),
            tagged(coverage, "WPB4_marker_coverage.csv", "marker_coverage"),
            tagged(pd.DataFrame([metadata]), "WPB4_run_metadata.json", "run_metadata"),
            tagged(pd.DataFrame(thumb_targets), "WPB4_overlay_index.csv", "overlay_thumbnails"),
        ]
    )
    return save_panel(
        fig,
        panel_id="S10",
        title="correction-layer image validation",
        basename="S10_correction_layer_validation",
        source=source,
        conclusion="Full-cohort FOV-level image validation supports global crosstalk plus production thresholds; per-marker recalibration/layer switching is not adopted.",
        data_sources="WPB4_marker_decision_table.csv; WPB4_marker_coverage.csv; WPB4_overlay_index.csv; WPB4_run_metadata.json",
    )


def panel_s11_cn_robustness() -> PanelExport:
    match = read_csv("WP1.2_k10_k11_matching.csv")
    no_sma = read_csv("WP1.2_exclude_noSMA_matching_to_current.csv")
    coverage = read_csv("WP1.2_cn_coverage.csv")
    defined_group = read_csv("WP1.1_defined_fraction_group_stats.csv")
    defined_comp = read_csv("WP1.1_defined_fraction_group_comparison.csv")
    missing = read_csv("WPB4_missing_SMA_PDL1_samples.csv")
    excluded_check = read_csv("WP1.2_noSMA_excluded_sample_celltype_check.csv")

    stable_counts = match["stability_call"].value_counts().reset_index()
    stable_counts.columns = ["stability_call", "n_k11_clusters"]
    nosma10 = no_sma[no_sma["k"].eq(10)].copy()

    fig = plt.figure(figsize=(8.2, 5.0))
    gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.38)
    ax = fig.add_subplot(gs[0, 0])
    sns.histplot(match["best_cosine_similarity"], bins=np.linspace(0.9, 1.0, 16), color=COLORS_15[1], ax=ax)
    ax.set_xlabel("best k11-to-k10 center cosine")
    ax.set_ylabel("clusters")
    ax.set_title("k10/k11 centers are stable")
    ax.text(0.02, 0.92, f"stable/split: {int((match['stability_call'].isin(['stable','split'])).sum())}/{len(match)}", transform=ax.transAxes, va="top")

    ax = fig.add_subplot(gs[0, 1])
    sns.barplot(data=stable_counts, x="stability_call", y="n_k11_clusters", ax=ax, color=COLORS_15[3])
    ax.set_xlabel("")
    ax.set_ylabel("k11 clusters")
    ax.set_title("k selection changes split HER2 niche")

    ax = fig.add_subplot(gs[1, 0])
    sns.stripplot(data=nosma10, x="k", y="cosine_similarity", ax=ax, color=COLORS_15[1], size=5, jitter=0.12)
    ax.axhline(0.98, color="#333333", ls="--", lw=0.9)
    ax.set_ylim(0.95, 1.003)
    ax.set_xlabel("exclude 7 no-SMA samples, k=10")
    ax.set_ylabel("center cosine to current CN")
    ax.set_title("noSMA means excluding samples")
    missing_samples = missing[missing["marker"].eq("SMA")]["sample"].tolist()
    missing_wrapped = "\n".join(
        [
            "excluded:",
            ", ".join(missing_samples[:3]),
            ", ".join(missing_samples[3:6]),
            ", ".join(missing_samples[6:]),
        ]
    )
    ax.text(0.02, 0.04, missing_wrapped, transform=ax.transAxes, fontsize=6, va="bottom")

    ax = fig.add_subplot(gs[1, 1])
    cov10 = coverage[coverage["k"].eq(10)].sort_values("CN")
    ax.bar(cov10["CN"], cov10["sample_coverage_fraction_ge_1pct"], color=[COLORS_15[i % len(COLORS_15)] for i in range(len(cov10))])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("sample coverage >=1%")
    ax.set_xlabel("CN")
    ax.set_title("CN coverage across samples")
    comp_row = defined_comp.iloc[0]
    ax.text(
        0.02,
        0.96,
        f"defined fraction CR/NR p={fmt(comp_row['welch_t_p'], 3)}",
        transform=ax.transAxes,
        va="top",
        fontsize=7,
    )

    fig.suptitle("S11. CN robustness and panel-layout sensitivity", x=0.01, y=0.995, ha="left", fontsize=9, fontweight="bold")
    source = combine_source(
        [
            tagged(match, "WP1.2_k10_k11_matching.csv", "k10_k11_center_matching"),
            tagged(no_sma, "WP1.2_exclude_noSMA_matching_to_current.csv", "exclude_noSMA_center_matching"),
            tagged(coverage, "WP1.2_cn_coverage.csv", "cn_coverage"),
            tagged(defined_group, "WP1.1_defined_fraction_group_stats.csv", "defined_fraction_group_stats"),
            tagged(defined_comp, "WP1.1_defined_fraction_group_comparison.csv", "defined_fraction_group_comparison"),
            tagged(missing, "WPB4_missing_SMA_PDL1_samples.csv", "missing_SMA_PDL1_samples"),
            tagged(excluded_check, "WP1.2_noSMA_excluded_sample_celltype_check.csv", "exclude_noSMA_sample_check"),
        ]
    )
    return save_panel(
        fig,
        panel_id="S11",
        title="CN robustness",
        basename="S11_cn_robustness",
        source=source,
        conclusion="CN interpretation is stable to k selection and to excluding the seven no-SMA/v1-like samples; noSMA is a sample-exclusion robustness analysis.",
        data_sources="WP1.2_k10_k11_matching.csv; WP1.2_exclude_noSMA_matching_to_current.csv; WP1.2_cn_coverage.csv; WP1.1_defined_fraction*; WPB4_missing_SMA_PDL1_samples.csv",
    )


def select_prediction_rows() -> pd.DataFrame:
    summary = read_csv("WP2.2_nested_prediction_summary.csv")
    rows = []
    specs = [
        ("celltype-only", (summary["feature_pool"] == "celltype_proportion_only") & (summary["model"] == "logistic_l2")),
        ("CN prop", (summary["feature_pool"] == "cn_proportion_only") & (summary["model"] == "logistic_l2")),
        ("external prior", (summary["feature_pool"] == "external_prior") & (summary["model"] == "external_prior_composite")),
        ("CN coupling top3", (summary["feature_pool"] == "full_candidate_pool") & (summary["top_n"].astype(str) == "3") & (summary["model"] == "fold_composite")),
    ]
    for label, mask in specs:
        sub = summary[mask & (summary["analysis_set"] == "primary_threshold30k")]
        if sub.empty:
            sub = summary[mask]
        if not sub.empty:
            row = sub.sort_values("auc", ascending=False).iloc[0].copy()
            row["ladder_label"] = label
            rows.append(row)
    return pd.DataFrame(rows)


def panel_s12_prediction_integrity() -> PanelExport:
    adj = read_csv("WP1.4_CRvsNR_stats.csv")
    keys = ["CN1__CN7", "CN3__CN7", "CN4__CN7"]
    coupling = adj[
        (adj["analysis_set"] == "primary_k10_threshold30k_panelaware")
        & (adj["metric"] == "z_score")
        & (adj["feature_id"].isin(keys))
    ].copy()
    coupling["plot_color"] = np.where(coupling["direction"].eq("higher_in_CR"), CR_COLOR, NR_COLOR)
    frag = adj[
        (adj["analysis_set"] == "primary_k10_threshold30k_panelaware")
        & (adj["feature_class"].isin(["fragmentation", "coherence"]))
    ].copy()
    frag["abs_d"] = frag["cohens_d_NR_minus_CR"].abs()
    frag_top = frag.sort_values("abs_d", ascending=False).head(5).copy()
    frag_top["plot_color"] = np.where(frag_top["cohens_d_NR_minus_CR"] < 0, CR_COLOR, NR_COLOR)

    leak = read_csv("WP2.2_leakage_comparison.csv")
    ladder = select_prediction_rows()
    null = read_csv("WP2.2_permutation_null.csv")
    null_top3 = null[
        (null["analysis_set"] == "primary_threshold30k")
        & (null["feature_pool"] == "full_candidate_pool")
        & (null["top_n"].astype(str) == "3")
        & (null["model"] == "fold_composite")
    ].copy()

    fig = plt.figure(figsize=(8.2, 5.6))
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.38)
    ax = fig.add_subplot(gs[0, 0])
    if not coupling.empty:
        plot_errorbarh(
            ax,
            coupling.sort_values("cohens_d_NR_minus_CR"),
            label_col="feature_id",
            effect_col="cohens_d_NR_minus_CR",
            low_col="cohens_d_bootstrap95_low",
            high_col="cohens_d_bootstrap95_high",
            color_col="plot_color",
        )
    ax.set_title("vascular niche-coupling axis")
    ax.set_xlabel("Cohen's d (NR - CR)")

    ax = fig.add_subplot(gs[0, 1])
    if not frag_top.empty:
        plot_errorbarh(
            ax,
            frag_top.sort_values("cohens_d_NR_minus_CR"),
            label_col="feature_name",
            effect_col="cohens_d_NR_minus_CR",
            low_col="cohens_d_bootstrap95_low",
            high_col="cohens_d_bootstrap95_high",
            color_col="plot_color",
        )
    ax.set_title("fragmentation/coherence candidates")
    ax.set_xlabel("Cohen's d (NR - CR)")

    ax = fig.add_subplot(gs[1, 0])
    leak_plot = leak.melt(
        id_vars=["comparison"],
        value_vars=["wp2_1_auc", "wp2_2_auc"],
        var_name="analysis",
        value_name="auc",
    )
    leak_plot["analysis"] = leak_plot["analysis"].map({"wp2_1_auc": "leaky WP2.1", "wp2_2_auc": "nested WP2.2"})
    sns.barplot(data=leak_plot, x="auc", y="comparison", hue="analysis", ax=ax, palette=["#AAAAAA", COLORS_15[1]])
    ax.set_xlim(0, 1)
    ax.set_xlabel("AUC")
    ax.set_ylabel("")
    ax.set_title("leakage comparison")
    ax.legend(loc="lower right", fontsize=6)

    ax = fig.add_subplot(gs[1, 1])
    if not ladder.empty:
        x = np.arange(len(ladder))
        ax.bar(x, ladder["auc"], color=[COLORS_15[2], COLORS_15[3], COLORS_15[7], COLORS_15[1]][: len(ladder)])
        for i, row in ladder.reset_index(drop=True).iterrows():
            ax.plot([i, i], [row["auc_bootstrap95_low"], row["auc_bootstrap95_high"]], color="#222222", lw=1)
            ax.text(i, min(0.98, row["auc"] + 0.035), fmt(row["auc"], 2), ha="center", fontsize=7)
        ax.set_xticks(x)
        ax.set_xticklabels(ladder["ladder_label"], rotation=30, ha="right")
    if not null_top3.empty:
        obs = float(null_top3["observed_auc"].iloc[0])
        ax.axhline(obs, color=CR_COLOR, lw=1.0, ls="--")
        ax.text(0.02, 0.93, f"top3 coupling observed AUC={fmt(obs, 3)}", transform=ax.transAxes, va="top", fontsize=7)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("AUC")
    ax.set_title("AUC ladder: celltype to coupling")

    fig.suptitle("S12. Spatial coupling and prediction integrity", x=0.01, y=0.995, ha="left", fontsize=9, fontweight="bold")
    source = combine_source(
        [
            tagged(coupling, "WP1.4_CRvsNR_stats.csv", "vascular_coupling_axis"),
            tagged(frag_top, "WP1.4_CRvsNR_stats.csv", "fragmentation_coherence_top"),
            tagged(leak, "WP2.2_leakage_comparison.csv", "leakage_comparison"),
            tagged(ladder, "WP2.2_nested_prediction_summary.csv", "auc_ladder"),
            tagged(null_top3.head(1000), "WP2.2_permutation_null.csv", "top3_permutation_null_first1000"),
        ]
    )
    return save_panel(
        fig,
        panel_id="S12",
        title="adjacency and prediction integrity",
        basename="S12_adjacency_prediction_integrity",
        source=source,
        conclusion="The only defensible predictive signal is the potential vascular niche-coupling axis; leakage comparisons show why old all-data composites are not valid.",
        data_sources="WP1.4_CRvsNR_stats.csv; WP2.2_nested_prediction_summary.csv; WP2.2_leakage_comparison.csv; WP2.2_permutation_null.csv",
    )


def build_supp_index(exports: list[PanelExport]) -> pd.DataFrame:
    expected = set(required_panel_ids())
    seen = {export.panel_id for export in exports}
    if seen != expected:
        raise ValueError(f"Expected panels {sorted(expected)}, got {sorted(seen)}")
    rows = []
    for export in exports:
        for ext in ["svg", "pdf", "png", "tiff", "source_data"]:
            path = export.paths[ext]
            if not path.exists() or path.stat().st_size <= 0:
                raise FileNotFoundError(f"Missing or empty {ext} for {export.panel_id}: {path}")
        rows.append(
            {
                "panel_id": export.panel_id,
                "title": export.title,
                "basename": export.basename,
                "conclusion": export.conclusion,
                "data_sources": export.data_sources,
                "source_rows": export.source_rows,
                "svg": str(export.paths["svg"]),
                "pdf": str(export.paths["pdf"]),
                "png": str(export.paths["png"]),
                "tiff": str(export.paths["tiff"]),
                "source_data": str(export.paths["source_data"]),
                "svg_bytes": export.paths["svg"].stat().st_size,
                "png_bytes": export.paths["png"].stat().st_size,
                "tiff_bytes": export.paths["tiff"].stat().st_size,
            }
        )
    order = {panel_id: idx for idx, panel_id in enumerate(required_panel_ids())}
    return pd.DataFrame(rows).sort_values("panel_id", key=lambda s: s.map(order)).reset_index(drop=True)


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    render = df.copy()
    for col in render.columns:
        render[col] = render[col].map(lambda x: fmt(x, 4) if isinstance(x, (float, np.floating)) else str(x))
    header = "| " + " | ".join(render.columns) + " |"
    divider = "| " + " | ".join("---" for _ in render.columns) + " |"
    rows = ["| " + " | ".join(str(row[col]) for col in render.columns) + " |" for _, row in render.iterrows()]
    return "\n".join([header, divider] + rows)


def write_report(index: pd.DataFrame) -> None:
    compact = index[["panel_id", "title", "conclusion", "data_sources", "source_rows"]].copy()
    report = f"""# WP3.3 supplementary Fig. 6 candidate panels

## Figure Contract

- Core conclusion: supplementary Fig. 6 carries the celltype-layer biology and method-integrity evidence behind the CN-centric main figure.
- Evidence chain: S8 QC/inclusion; S9 celltype biology; S10 correction-layer image validation; S11 CN robustness and noSMA sample-exclusion robustness; S12 adjacency, leakage comparison, and AUC ladder.
- Archetype: quantitative grid plus image-validation plate.
- Backend/export: Python-only matplotlib; SVG/PDF/PNG/TIFF exports plus per-panel source data CSV.
- Integrity language: candidate-level associations, potential prediction, noSMA means excluding seven samples, not removing an SMA feature.

## Candidate Index

{markdown_table(compact)}

## Output Directory

`{PANEL_DIR}`
"""
    OUT_REPORT.write_text(report, encoding="utf-8")


def run() -> None:
    configure_matplotlib()
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    exports = [
        panel_s8_qc(),
        panel_s9_celltype_biology(),
        panel_s10_correction_validation(),
        panel_s11_cn_robustness(),
        panel_s12_prediction_integrity(),
    ]
    index = build_supp_index(exports)
    index.to_csv(OUT_INDEX, index=False)
    write_report(index)
    metadata = {
        "panel_dir": str(PANEL_DIR),
        "index": str(OUT_INDEX),
        "report": str(OUT_REPORT),
        "n_panels": int(len(index)),
        "panels": index[["panel_id", "basename", "source_rows", "svg_bytes", "png_bytes"]].to_dict(orient="records"),
        "backend": "python/matplotlib",
        "exports": ["svg", "pdf", "png", "tiff", "source_data"],
    }
    OUT_METADATA.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


if __name__ == "__main__":
    run()
