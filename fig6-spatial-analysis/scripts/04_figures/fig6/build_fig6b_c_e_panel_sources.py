#!/usr/bin/env python3
"""WP3.1b Fig. 6 revision panels.

Revision scope:
- manuscript colors_15 palette,
- crop-based cell-level spatial panels with phenotype/CN rows,
- rasterized scatter artists for lightweight SVG/PDF,
- revised 6d sample selection by headline vascular-coupling metrics.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from PIL import Image
from sklearn.metrics import roc_auc_score, roc_curve


PROJECT_DIR = Path("/data/home/chenwenhui/mIHC/fig6_pm")
DELIVERABLE_DIR = PROJECT_DIR / "deliverables"
CELLINFO_DIR = Path("/data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh")
NEIGHBORHOOD_DIR = CELLINFO_DIR / "neighborhood_radius200_B0"
PANEL_DIR = CELLINFO_DIR / "fig6_panels_rev"
OLD_PANEL_DIR = CELLINFO_DIR / "fig6_panels"
SCRATCH_DIR = PROJECT_DIR / "scratch_wp3.1b"

PARQUET_PATH = NEIGHBORHOOD_DIR / "neighborhood_composition.parquet"
OUT_INDEX = DELIVERABLE_DIR / "WP3.1b_fig6_candidate_index.csv"
OUT_REPORT = DELIVERABLE_DIR / "WP3.1b_fig6_revision_report.md"
OUT_SCRIPT = DELIVERABLE_DIR / "WP3_1b_fig6_revision.py"
OUT_METADATA = DELIVERABLE_DIR / "WP3.1b_fig6_revision_metadata.json"

PRIMARY_PROP_ANALYSIS = "primary_k10_threshold30k_panelaware"
PIXEL_SIZE_UM = 0.4995
SCALE_BAR_UM = 100
SCALE_BAR_PX = SCALE_BAR_UM / PIXEL_SIZE_UM

NO_SMA_SAMPLES = {
    "pt16-pre",
    "pt17-pre",
    "pt19-pre",
    "pt20-pre",
    "pt21-pre",
    "pt23-post",
    "pt23-pre",
}

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
CN_ORDER = [f"CN{i}" for i in range(10)]
CN_COLORS = {cn: COLORS_15[i] for i, cn in enumerate(CN_ORDER)}
CELLTYPE_COLORS = {
    "B cell": "#A52A2A",
    "CD163+ Macrophage": "#053E7A",
    "CD163- Macrophage": "#E46B4F",
    "CD21+ B cell": "#8B8B00",
    "CD4 T cell": "#B8860B",
    "CD8 T cell": "#CD661D",
    "CD8 exhausted T cell": "#8B5A2B",
    "DC": "#CDAD00",
    "Vascular": "#2F4F4F",
    "Fibroblast / SMA+ stromal": "#483D8B",
    "HER2+ Tumor": "#7A378B",
    "Treg": "#4A7088",
    "PanCK+ Tumor": "#CD7054",
}
GROUP_COLORS = {"CR": "#E4212B", "NR": "#2C5F9E"}
TIME_COLORS = {"Pre": "#E08A3C", "Post": "#4A78B5"}
NEUTRAL_CELL = "#D9D9D9"
NEUTRAL_FAINT = "#ECECEC"

# 空间 crop 的点样式:重点 CN/celltype 放大加实,背景细胞缩小压淡。
# 与 WP3_2_fig6_assembly.CONFIG['crop_point_style'] 保持同一套默认值。
CROP_POINT_STYLE = {
    "highlight_size_scale": 2.0,
    "other_keep_color": True,    # 背景点保留类别原色(只缩小压淡);False = 压成中性灰
    "other_size_scale": 0.8,
    "highlight_alpha": 0.95,
    "other_alpha": 0.28,
}
MANUSCRIPT_BWR = LinearSegmentedColormap.from_list("manuscript_bwr", ["#053E7A", "#FFFFFF", "#A52A2A"])

CN_SHORT = {
    "CN0": "Mac-rich",
    "CN1": "HER2 tumor-Mac",
    "CN2": "CD4/Treg",
    "CN3": "PanCK tumor",
    "CN4": "CD8 cytotoxic",
    "CN5": "CD21+B/TLS",
    "CN6": "HER2 dense",
    "CN7": "Vascular",
    "CN8": "B tumor-immune",
    "CN9": "Fibro-stromal",
}


@dataclass(frozen=True)
class FigureExport:
    basename: str
    files: dict[str, Path]
    source_rows: int
    svg_bytes: int


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
            "font.size": 7,
            "axes.edgecolor": "#333333",
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.dpi": 140,
        }
    )


def required_export_paths(panel_dir: Path | str, basename: str) -> dict[str, Path]:
    panel_dir = Path(panel_dir)
    return {
        "svg": panel_dir / f"{basename}.svg",
        "pdf": panel_dir / f"{basename}.pdf",
        "png": panel_dir / f"{basename}.png",
        "tiff": panel_dir / f"{basename}.tiff",
        "source_data": panel_dir / f"{basename}_source_data.csv",
    }


def save_figure(fig: plt.Figure, basename: str, source: pd.DataFrame, *, dpi: int = 600) -> FigureExport:
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    paths = required_export_paths(PANEL_DIR, basename)
    source.to_csv(paths["source_data"], index=False)
    fig.savefig(paths["svg"], bbox_inches="tight", dpi=dpi)
    fig.savefig(paths["pdf"], bbox_inches="tight", dpi=dpi)
    fig.savefig(paths["png"], bbox_inches="tight", dpi=300)
    fig.savefig(paths["tiff"], bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    return FigureExport(basename, paths, int(len(source)), int(paths["svg"].stat().st_size))


def robust_scale(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if len(values) < 2:
        return 1.0
    iqr = float(values.quantile(0.75) - values.quantile(0.25))
    if np.isfinite(iqr) and iqr > 0:
        return iqr / 1.349
    std = float(values.std(ddof=1))
    return std if np.isfinite(std) and std > 0 else 1.0


def select_representative_samples(
    metrics: pd.DataFrame,
    *,
    metric_cols: list[str],
    per_group: int = 3,
    group_col: str = "Outcome_short",
) -> pd.DataFrame:
    scales = {col: robust_scale(metrics[col]) for col in metric_cols}
    parts = []
    for group in ["CR", "NR"]:
        sub = metrics.loc[metrics[group_col].eq(group)].copy()
        med = sub[metric_cols].median(numeric_only=True)
        distances = []
        for _, row in sub.iterrows():
            chunks = []
            for col in metric_cols:
                value = pd.to_numeric(pd.Series([row[col]]), errors="coerce").iloc[0]
                if np.isfinite(value) and np.isfinite(med[col]):
                    chunks.append(((float(value) - float(med[col])) / scales[col]) ** 2)
            distances.append(float(np.sqrt(np.mean(chunks))) if chunks else np.inf)
        sub["selection_distance"] = distances
        sub["selection_target_median"] = ";".join(f"{col}={med[col]:.4g}" for col in metric_cols)
        sub["selection_basis"] = "robust distance to group median: " + ";".join(metric_cols)
        sub = sub.sort_values(["selection_distance", "sample"]).head(per_group)
        sub["selection_rank"] = np.arange(1, len(sub) + 1)
        parts.append(sub)
    return pd.concat(parts, ignore_index=True)


def select_headline_effect_side_samples(
    metrics: pd.DataFrame,
    *,
    group_specs: dict[str, dict[str, object]],
    per_group: int = 3,
) -> pd.DataFrame:
    """Select samples on the expected effect side closest to the group median."""
    parts = []
    for group in ["CR", "NR"]:
        spec = group_specs[group]
        metric = str(spec["metric"])
        direction = int(spec.get("direction", 1))
        sub = metrics.loc[metrics["Outcome_short"].eq(group)].dropna(subset=[metric]).copy()
        median = float(sub[metric].median())
        if direction >= 0:
            effect_side = sub[sub[metric] >= median].copy()
        else:
            effect_side = sub[sub[metric] <= median].copy()
        if len(effect_side) < per_group:
            effect_side = sub.copy()
        scale = robust_scale(sub[metric])
        effect_side["selection_distance"] = (effect_side[metric] - median).abs() / scale
        effect_side["selection_target_median"] = f"{metric} group median={median:.4g}; effect-side filter direction={direction:+d}"
        effect_side["selection_basis"] = f"headline metric effect-side distance to group median: {metric}"
        effect_side = effect_side.sort_values(["selection_distance", "sample"]).head(per_group)
        effect_side["selection_rank"] = np.arange(1, len(effect_side) + 1)
        parts.append(effect_side)
    return pd.concat(parts, ignore_index=True)


def load_niche_map() -> pd.DataFrame:
    return pd.read_csv(DELIVERABLE_DIR / "WP1.2_functional_niche_map_k10.csv").sort_values("CN")


def load_primary_meta(prop: pd.DataFrame) -> pd.DataFrame:
    meta_cols = [
        "sample",
        "Time_short",
        "Outcome_short",
        "Group",
        "passes_lowcell_30000",
        "passes_lowcell_20000",
        "cn_participant_cells",
        "panel_layout",
        "is_v1_like",
    ]
    meta = prop[meta_cols].drop_duplicates("sample")
    mask = meta["Time_short"].eq("Pre") & meta["Outcome_short"].isin(["CR", "NR"]) & meta["passes_lowcell_30000"].astype(bool)
    return meta.loc[mask].sort_values("sample").reset_index(drop=True)


def series_by_sample(df: pd.DataFrame, key_col: str, key: str, value_col: str, samples: list[str]) -> pd.Series:
    rows = df.loc[df[key_col].astype(str).eq(str(key)), ["sample", value_col]].copy()
    return rows.drop_duplicates("sample").set_index("sample")[value_col].astype(float).reindex(samples)


def build_sample_metrics() -> pd.DataFrame:
    prop = pd.read_csv(DELIVERABLE_DIR / "WP1.3_sample_proportion_long.csv")
    meta = load_primary_meta(prop)
    samples = meta["sample"].tolist()
    out = meta.set_index("sample").copy()
    cn_rows = prop[prop["feature_domain"].eq("CN") & prop["k"].astype(str).eq("10")].copy()
    cell_rows = prop[prop["feature_domain"].eq("celltype") & prop["k"].astype(str).eq("celltype")].copy()
    for cn in CN_ORDER:
        out[f"{cn.lower()}_prop"] = series_by_sample(cn_rows, "feature_id", cn, "proportion", samples)
    out["cd21b_prop"] = series_by_sample(cell_rows, "feature_id", "CD21+ B cell", "proportion", samples)
    out["cd8_prop"] = series_by_sample(cell_rows, "feature_id", "CD8 T cell", "proportion", samples)
    out["her2_prop"] = series_by_sample(cell_rows, "feature_id", "HER2+ Tumor", "proportion", samples)
    out["panck_prop"] = series_by_sample(cell_rows, "feature_id", "PanCK+ Tumor", "proportion", samples)

    frag = pd.read_csv(DELIVERABLE_DIR / "WP1.4_fragmentation_coherence.csv")
    frag = frag[frag["k"].eq(10)].copy()
    cn5 = frag[frag["CN"].eq("CN5")].drop_duplicates("sample").set_index("sample")
    for metric in ["mean_patch_size", "frag_components_per_cell", "cn_coherence", "lcc_fraction"]:
        out[f"cn5_{metric}"] = cn5[metric].astype(float).reindex(samples)

    adj = pd.read_csv(DELIVERABLE_DIR / "WP1.4_adjacency_sample_metrics.csv")
    adj = adj[adj["k"].eq(10)].copy()
    for pair in ["CN1__CN7", "CN3__CN7", "CN4__CN7", "CN5__CN7"]:
        pair_df = adj[adj["CN_pair"].eq(pair)].drop_duplicates("sample").set_index("sample")
        key = pair.replace("__", "_")
        out[f"coupling_z_{key}"] = pair_df["z_score"].astype(float).reindex(samples)
        out[f"coupling_roe_{key}"] = pair_df["RoE"].astype(float).reindex(samples)

    coh = pd.read_csv(DELIVERABLE_DIR / "WP1.5_coherence_sample_metrics.csv")
    cd21 = coh[coh["celltype"].eq("CD21+ B cell")].drop_duplicates("sample").set_index("sample")
    out["cd21b_same_type_knn_fraction"] = cd21["same_type_knn_fraction"].astype(float).reindex(samples)
    out["cd21b_binary_moran_i_knn"] = cd21["binary_moran_i_knn"].astype(float).reindex(samples)
    return out.reset_index()


def load_spatial_for_samples(samples: list[str]) -> pd.DataFrame:
    wanted = set(samples)
    cols = ["sample", "spatial_x", "spatial_y", "celltype", "CN_k10"]
    parts = []
    parquet = pq.ParquetFile(PARQUET_PATH)
    for row_group in range(parquet.num_row_groups):
        df = parquet.read_row_group(row_group, columns=cols).to_pandas()
        df = df[df["sample"].isin(wanted)]
        if not df.empty:
            parts.append(df)
    if not parts:
        raise RuntimeError(f"no spatial data for {samples}")
    spatial = pd.concat(parts, ignore_index=True).dropna(subset=["spatial_x", "spatial_y", "CN_k10"])
    spatial["CN_k10"] = spatial["CN_k10"].astype(str)
    spatial["celltype"] = spatial["celltype"].astype(str)
    return spatial


def assign_plot_colors(
    df: pd.DataFrame,
    column: str,
    highlights: set[str],
    palette: dict[str, str],
) -> pd.Series:
    return df[column].astype(str).map(lambda value: palette.get(value, NEUTRAL_CELL) if value in highlights else NEUTRAL_CELL)


def select_crop_window(
    df: pd.DataFrame,
    *,
    crop_size_px: int,
    cn_targets: list[str],
    cell_targets: list[str],
    cn_weights: dict[str, float] | None = None,
    cell_weights: dict[str, float] | None = None,
    require_first_two: bool = False,
    interface_weight: float = 3.0,
    grid_step_px: int | None = None,
) -> dict[str, float]:
    """Find a deterministic crop with target niche/celltype co-localization."""
    grid_step_px = grid_step_px or max(120, crop_size_px // 4)
    key = df[df["CN_k10"].isin(cn_targets) | df["celltype"].isin(cell_targets)].copy()
    if key.empty:
        key = df.copy()
    x_min, x_max = float(df["spatial_x"].min()), float(df["spatial_x"].max())
    y_min, y_max = float(df["spatial_y"].min()), float(df["spatial_y"].max())
    half = crop_size_px / 2
    xs = np.arange(max(x_min + half, key["spatial_x"].quantile(0.05)), min(x_max - half, key["spatial_x"].quantile(0.95)) + 1, grid_step_px)
    ys = np.arange(max(y_min + half, key["spatial_y"].quantile(0.05)), min(y_max - half, key["spatial_y"].quantile(0.95)) + 1, grid_step_px)
    if len(xs) == 0:
        xs = np.array([float(key["spatial_x"].median())])
    if len(ys) == 0:
        ys = np.array([float(key["spatial_y"].median())])
    cn_weights = cn_weights or {cn: 1.0 for cn in cn_targets}
    cell_weights = cell_weights or {cell: 1.0 for cell in cell_targets}
    best = None
    for cx in xs:
        in_x = (df["spatial_x"] >= cx - half) & (df["spatial_x"] <= cx + half)
        for cy in ys:
            window = df.loc[in_x & (df["spatial_y"] >= cy - half) & (df["spatial_y"] <= cy + half)]
            if window.empty:
                continue
            cn_counts = {cn: int(window["CN_k10"].eq(cn).sum()) for cn in cn_targets}
            cell_counts = {ct: int(window["celltype"].eq(ct).sum()) for ct in cell_targets}
            target_total = sum(float(cn_weights.get(cn, 1.0)) * min(count, 1500) for cn, count in cn_counts.items())
            target_total += 0.55 * sum(float(cell_weights.get(ct, 1.0)) * min(count, 1500) for ct, count in cell_counts.items())
            interface = 0.0
            if len(cn_targets) >= 2:
                first_two_min = min(cn_counts.get(cn_targets[0], 0), cn_counts.get(cn_targets[1], 0), 1500)
                interface += first_two_min * float(interface_weight)
            if len(cn_targets) >= 3:
                interface += min(cn_counts.get(cn_targets[0], 0), cn_counts.get(cn_targets[2], 0), 1500) * 1.0
            density_bonus = min(len(window), 18000) * 0.015
            score = float(target_total + interface + density_bonus)
            if require_first_two and len(cn_targets) >= 2 and first_two_min == 0:
                score *= 0.05
            candidate = {
                "center_x": float(cx),
                "center_y": float(cy),
                "xmin": float(cx - half),
                "xmax": float(cx + half),
                "ymin": float(cy - half),
                "ymax": float(cy + half),
                "score": score,
                "n_cells_window": int(len(window)),
                **{f"n_{cn}": cn_counts.get(cn, 0) for cn in cn_targets},
                **{f"n_{ct}": cell_counts.get(ct, 0) for ct in cell_targets},
            }
            if best is None or (candidate["score"], -abs(cx), -abs(cy)) > (best["score"], -abs(best["center_x"]), -abs(best["center_y"])):
                best = candidate
    if best is None:
        cx, cy = float(key["spatial_x"].median()), float(key["spatial_y"].median())
        best = {
            "center_x": cx,
            "center_y": cy,
            "xmin": cx - half,
            "xmax": cx + half,
            "ymin": cy - half,
            "ymax": cy + half,
            "score": 0.0,
            "n_cells_window": 0,
        }
    return best


def crop_dataframe(df: pd.DataFrame, crop: dict[str, float]) -> pd.DataFrame:
    return df[
        (df["spatial_x"] >= crop["xmin"])
        & (df["spatial_x"] <= crop["xmax"])
        & (df["spatial_y"] >= crop["ymin"])
        & (df["spatial_y"] <= crop["ymax"])
    ].copy()


def draw_scale_bar(ax: plt.Axes, crop: dict[str, float]) -> None:
    x0 = crop["xmax"] - SCALE_BAR_PX - 35
    y0 = crop["ymax"] - 35
    ax.plot([x0, x0 + SCALE_BAR_PX], [y0, y0], color="black", lw=1.4, solid_capstyle="butt", zorder=20)
    ax.text(x0 + SCALE_BAR_PX / 2, y0 - 18, f"{SCALE_BAR_UM} µm", ha="center", va="top", fontsize=6.4, color="black")


def add_locator(ax: plt.Axes, sample_df: pd.DataFrame, crop: dict[str, float], color: str) -> None:
    inset = inset_axes(ax, width="24%", height="24%", loc="upper right", borderpad=0.45)
    if len(sample_df) > 12000:
        sample = sample_df.sample(n=12000, random_state=17)
    else:
        sample = sample_df
    inset.scatter(sample["spatial_x"], sample["spatial_y"], s=0.08, color="#CFCFCF", linewidths=0, rasterized=True)
    inset.add_patch(
        Rectangle(
            (crop["xmin"], crop["ymin"]),
            crop["xmax"] - crop["xmin"],
            crop["ymax"] - crop["ymin"],
            fill=False,
            edgecolor=color,
            linewidth=0.8,
        )
    )
    inset.set_aspect("equal")
    inset.invert_yaxis()
    inset.set_xticks([])
    inset.set_yticks([])
    for spine in inset.spines.values():
        spine.set_color(color)
        spine.set_linewidth(0.5)


def scatter_crop(
    ax: plt.Axes,
    crop_df: pd.DataFrame,
    *,
    color_col: str,
    highlights: Iterable[str] | None = None,
    palette: dict[str, str],
    crop: dict[str, float],
    point_size: float = 1.6,
    style: dict | None = None,
    emphasis_col: str | None = None,
    emphasis_values: Iterable[str] | None = None,
    other_keep_color: bool | None = None,
) -> None:
    """重点细胞画大画实,其余画小画淡(与 WP3.2 assembly 同一套口径)。

    emphasis_col / emphasis_values
        决定"谁被放大"的依据列。默认沿用 color_col + highlights,即各图按各自的列判定。
        传 emphasis_col="CN_k10" 则 phenotype 图也按 CN 判定放大 —— 两张图强调的是
        同一批细胞,可以上下对读(scimap_IT.ipynb 的 export_manuscript_roi_v4 就是这个口径)。
        此时颜色仍按 color_col 上,所以 phenotype 图显示的是该 niche 的细胞构成。

    other_keep_color
        背景点是否保留自己的类别颜色(只缩小压淡),False 则统一压成中性灰。
    """
    st = {**CROP_POINT_STYLE, **(style or {})}
    keep_color = st["other_keep_color"] if other_keep_color is None else other_keep_color
    hl = None if highlights is None else set(highlights)

    if emphasis_col is not None:
        mask = crop_df[emphasis_col].astype(str).isin(set(emphasis_values or hl or []))
    else:
        mask = crop_df[color_col].astype(str).isin(hl or set())
    fore, back = crop_df[mask], crop_df[~mask]

    def _color(label: str) -> str:
        # highlights=None -> 所有类别都上自己的色(参考图口径);否则只有 highlights 里的上色
        if hl is not None and label not in hl:
            return NEUTRAL_CELL
        return palette.get(label, NEUTRAL_CELL)

    def _draw(part, *, size_scale: str, alpha_key: str, zorder: int, flat: bool):
        if part.empty:
            return
        if flat:
            ax.scatter(
                part["spatial_x"], part["spatial_y"],
                s=point_size * float(st[size_scale]), color=NEUTRAL_FAINT,
                alpha=float(st[alpha_key]), linewidths=0, rasterized=True, zorder=zorder,
            )
            return
        for label, grp in part.groupby(part[color_col].astype(str), sort=True):
            ax.scatter(
                grp["spatial_x"], grp["spatial_y"],
                s=point_size * float(st[size_scale]), color=_color(label),
                alpha=float(st[alpha_key]), linewidths=0, rasterized=True, zorder=zorder,
            )

    _draw(back, size_scale="other_size_scale", alpha_key="other_alpha",
          zorder=2, flat=not keep_color)
    _draw(fore, size_scale="highlight_size_scale", alpha_key="highlight_alpha",
          zorder=3, flat=False)
    ax.set_xlim(crop["xmin"], crop["xmax"])
    ax.set_ylim(crop["ymax"], crop["ymin"])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#CCCCCC")
        spine.set_linewidth(0.45)


def plot_workflow_placeholder() -> FigureExport:
    basename = "fig6a_workflow_schematic_placeholder"
    source = pd.DataFrame(
        [
            {"panel": "6a", "role": "schematic placeholder", "note": "PI will redraw with BioRender; WP3.1b keeps manuscript-color placeholder only."},
            {"panel": "6b", "role": "CN annotation", "note": "k10 composition heatmap using manuscript blue-white-red."},
            {"panel": "6c/6d", "role": "spatial crops", "note": "cell-level phenotype and CN rows with locator and scale bar."},
            {"panel": "6e/6f", "role": "quantification/prediction", "note": "sample-level effect sizes and unchanged WP2.2 ROC."},
        ]
    )
    fig, ax = plt.subplots(figsize=(9.4, 1.8))
    ax.set_axis_off()
    ax.text(0.01, 0.96, "a", fontsize=13, weight="bold", transform=ax.transAxes, va="top")
    ax.text(0.5, 0.62, "Fig. 6 schematic placeholder", ha="center", va="center", fontsize=13, weight="bold", transform=ax.transAxes)
    ax.text(
        0.5,
        0.38,
        "CN territories -> immune/TLS organization -> vascular niche coupling -> potential baseline prediction",
        ha="center",
        va="center",
        fontsize=8.4,
        color="#333333",
        transform=ax.transAxes,
    )
    ax.text(0.5, 0.16, "schematic placeholder; BioRender redraw planned", ha="center", va="center", fontsize=7.2, color="#666666", transform=ax.transAxes)
    return save_figure(fig, basename, source, dpi=600)


def plot_composition_heatmap(niche_map: pd.DataFrame) -> FigureExport:
    basename = "fig6b_cn_composition_heatmap"
    z = pd.read_csv(DELIVERABLE_DIR / "WP1.2_cn_composition_current_k10_heatmap_zscore_matrix.csv")
    comp = pd.read_csv(DELIVERABLE_DIR / "WP1.2_cn_composition_current_k10_stacked_bar_composition_matrix.csv")
    z = z.rename(columns={z.columns[0]: "cell_phenotype"}).set_index("cell_phenotype").loc[:, CN_ORDER]
    comp = comp.rename(columns={comp.columns[0]: "cell_phenotype"}).set_index("cell_phenotype").loc[:, CN_ORDER]
    source = (
        z.reset_index()
        .melt(id_vars="cell_phenotype", var_name="CN", value_name="z_score")
        .merge(comp.reset_index().melt(id_vars="cell_phenotype", var_name="CN", value_name="composition_fraction"), on=["cell_phenotype", "CN"], how="left")
        .merge(niche_map[["CN", "functional_niche", "niche_family", "dominant_celltype"]], on="CN", how="left")
    )

    fig = plt.figure(figsize=(9.1, 4.7))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.25, 1.25], wspace=0.32)
    ax = fig.add_subplot(gs[0, 0])
    sns.heatmap(
        z,
        ax=ax,
        cmap=MANUSCRIPT_BWR,
        norm=TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3),
        cbar_kws={"label": "Z-score", "shrink": 0.62, "ticks": [-2, 0, 3]},
        linewidths=0.35,
        linecolor="white",
        xticklabels=[f"{cn}\n{CN_SHORT[cn]}" for cn in CN_ORDER],
    )
    ax.set_xlabel("Cellular neighborhoods (k=10)")
    ax.set_ylabel("Cell phenotypes")
    ax.tick_params(axis="x", rotation=45)
    for label in ax.get_xticklabels():
        label.set_ha("right")

    ax2 = fig.add_subplot(gs[0, 1])
    bottom = np.zeros(len(CN_ORDER))
    for celltype in comp.index:
        values = comp.loc[celltype, CN_ORDER].to_numpy(dtype=float)
        ax2.bar(np.arange(len(CN_ORDER)), values, bottom=bottom, width=0.76, color=CELLTYPE_COLORS.get(celltype, "#AAAAAA"), edgecolor="none", label=celltype)
        bottom += values
    ax2.set_ylim(0, 1)
    ax2.set_xticks(np.arange(len(CN_ORDER)))
    ax2.set_xticklabels(CN_ORDER, rotation=45, ha="right")
    ax2.set_ylabel("Composition fraction")
    ax2.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=5.9, title="Phenotype", title_fontsize=6.3)
    fig.text(0.01, 0.98, "b", fontsize=13, weight="bold", va="top")
    fig.subplots_adjust(left=0.08, right=0.84, top=0.96, bottom=0.20)
    return save_figure(fig, basename, source, dpi=600)


def spatial_panel_specs(panel: str, outcome: str) -> dict[str, object]:
    if panel == "6c":
        if outcome == "CR":
            return {
                "cn_targets": ["CN5", "CN4", "CN6"],
                "cell_targets": ["CD21+ B cell", "CD8 T cell", "HER2+ Tumor", "PanCK+ Tumor"],
                "cn_weights": {"CN5": 8.0, "CN4": 3.0, "CN6": 0.8},
                "cell_weights": {"CD21+ B cell": 5.0, "CD8 T cell": 2.0, "HER2+ Tumor": 0.35, "PanCK+ Tumor": 0.35},
                "require_first_two": True,
                "interface_weight": 8.0,
                "cn_highlights": {"CN5", "CN4", "CN6", "CN3", "CN1"},
                "cell_highlights": {"CD21+ B cell", "CD8 T cell", "HER2+ Tumor", "PanCK+ Tumor"},
                "crop_size_px": 1250,
            }
        return {
            "cn_targets": ["CN6", "CN3", "CN4", "CN5"],
            "cell_targets": ["HER2+ Tumor", "PanCK+ Tumor", "CD8 T cell", "CD21+ B cell"],
            "cn_weights": {"CN6": 3.0, "CN3": 3.0, "CN4": 1.0, "CN5": 2.0},
            "cell_weights": {"HER2+ Tumor": 1.0, "PanCK+ Tumor": 1.0, "CD8 T cell": 1.0, "CD21+ B cell": 2.0},
            "require_first_two": False,
            "interface_weight": 4.0,
            "cn_highlights": {"CN5", "CN4", "CN6", "CN3", "CN1"},
            "cell_highlights": {"CD21+ B cell", "CD8 T cell", "HER2+ Tumor", "PanCK+ Tumor"},
            "crop_size_px": 1250,
        }
    if outcome == "CR":
        return {
            "cn_targets": ["CN7", "CN1"],
            "cell_targets": ["Vascular", "HER2+ Tumor"],
            "cn_weights": {"CN7": 4.0, "CN1": 5.0},
            "cell_weights": {"Vascular": 4.0, "HER2+ Tumor": 2.0},
            "require_first_two": True,
            "interface_weight": 12.0,
            "cn_highlights": {"CN7", "CN1"},
            "cell_highlights": {"Vascular", "HER2+ Tumor", "CD8 T cell"},
            "crop_size_px": 1150,
        }
    return {
        "cn_targets": ["CN7", "CN3"],
        "cell_targets": ["Vascular", "PanCK+ Tumor"],
        "cn_weights": {"CN7": 4.0, "CN3": 5.0},
        "cell_weights": {"Vascular": 4.0, "PanCK+ Tumor": 2.0},
        "require_first_two": True,
        "interface_weight": 12.0,
        "cn_highlights": {"CN7", "CN3"},
        "cell_highlights": {"Vascular", "PanCK+ Tumor", "HER2+ Tumor"},
        "crop_size_px": 1150,
    }


def plot_spatial_crop_panel(
    *,
    basename: str,
    letter: str,
    title: str,
    panel_id: str,
    selected: pd.DataFrame,
    spatial: pd.DataFrame,
    metrics: pd.DataFrame,
) -> tuple[FigureExport, pd.DataFrame]:
    selected = selected.copy().reset_index(drop=True)
    samples = selected["sample"].tolist()
    ncols = len(samples)
    fig = plt.figure(figsize=(max(9.8, ncols * 1.85), 4.3))
    gs = fig.add_gridspec(2, ncols, height_ratios=[1, 1], hspace=0.04, wspace=0.045)
    crop_records = []
    source_parts = []
    metric_lookup = metrics.set_index("sample")
    for col, sample in enumerate(samples):
        row = selected[selected["sample"].eq(sample)].iloc[0]
        sample_df = spatial[spatial["sample"].eq(sample)].copy()
        spec = spatial_panel_specs(panel_id, str(row["Outcome_short"]))
        crop = select_crop_window(
            sample_df,
            crop_size_px=int(spec["crop_size_px"]),
            cn_targets=list(spec["cn_targets"]),
            cell_targets=list(spec["cell_targets"]),
            cn_weights=dict(spec.get("cn_weights", {})),
            cell_weights=dict(spec.get("cell_weights", {})),
            require_first_two=bool(spec.get("require_first_two", False)),
            interface_weight=float(spec.get("interface_weight", 3.0)),
        )
        crop_df = crop_dataframe(sample_df, crop)
        crop_df = crop_df.sort_values(["CN_k10", "celltype"]).copy()
        if len(crop_df) > 30000:
            crop_df = crop_df.sample(n=30000, random_state=20260708 + col).sort_values(["CN_k10", "celltype"])
        crop_row = {
            "panel": panel_id,
            "sample": sample,
            "Outcome_short": row["Outcome_short"],
            "selection_rank": int(row["selection_rank"]),
            "selection_basis": row["selection_basis"],
            **{k: v for k, v in crop.items() if isinstance(v, (int, float, np.integer, np.floating))},
        }
        for metric in ["cn5_prop", "cd21b_same_type_knn_fraction", "coupling_z_CN1_CN7", "coupling_z_CN3_CN7", "coupling_z_CN4_CN7"]:
            if metric in metric_lookup.columns:
                crop_row[metric] = metric_lookup.loc[sample, metric]
                crop_df[metric] = metric_lookup.loc[sample, metric]
        crop_records.append(crop_row)
        crop_df["panel"] = panel_id
        crop_df["crop_sample"] = sample
        crop_df["crop_xmin"] = crop["xmin"]
        crop_df["crop_xmax"] = crop["xmax"]
        crop_df["crop_ymin"] = crop["ymin"]
        crop_df["crop_ymax"] = crop["ymax"]
        crop_df["selection_rank"] = int(row["selection_rank"])
        source_parts.append(crop_df)

        ax_ph = fig.add_subplot(gs[0, col])
        ax_cn = fig.add_subplot(gs[1, col])
        scatter_crop(
            ax_ph,
            crop_df,
            color_col="celltype",
            highlights=set(spec["cell_highlights"]),
            palette=CELLTYPE_COLORS,
            crop=crop,
            point_size=1.65,
        )
        scatter_crop(
            ax_cn,
            crop_df,
            color_col="CN_k10",
            highlights=set(spec["cn_highlights"]),
            palette=CN_COLORS,
            crop=crop,
            point_size=1.65,
        )
        add_locator(ax_ph, sample_df, crop, GROUP_COLORS[str(row["Outcome_short"])])
        draw_scale_bar(ax_cn, crop)
        ax_ph.set_title(f"{sample} {row['Outcome_short']} #{int(row['selection_rank'])}", color=GROUP_COLORS[str(row["Outcome_short"])], fontsize=8.2, pad=3)
        if col == 0:
            ax_ph.text(-0.05, 0.5, "phenotype", rotation=0, ha="right", va="center", transform=ax_ph.transAxes, fontsize=7.5)
            ax_cn.text(-0.05, 0.5, "CN", rotation=0, ha="right", va="center", transform=ax_cn.transAxes, fontsize=7.5)
        annotation = ""
        if panel_id == "6c":
            annotation = f"CN5 {metric_lookup.loc[sample, 'cn5_prop']:.3f}"
        elif row["Outcome_short"] == "CR":
            annotation = f"z CN1-CN7 {metric_lookup.loc[sample, 'coupling_z_CN1_CN7']:.2f}"
        else:
            annotation = f"z CN3-CN7 {metric_lookup.loc[sample, 'coupling_z_CN3_CN7']:.2f}"
        ax_cn.text(0.02, 0.05, annotation, transform=ax_cn.transAxes, ha="left", va="bottom", fontsize=6.6, bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=1.5))

    fig.text(0.01, 0.99, letter, fontsize=13, weight="bold", va="top")
    fig.suptitle(title, y=0.995, fontsize=10.5, weight="bold")
    legend_cn = []
    if panel_id == "6c":
        cn_items = ["CN5", "CN4", "CN6", "CN3", "CN1"]
        cell_items = ["CD21+ B cell", "CD8 T cell", "HER2+ Tumor", "PanCK+ Tumor"]
    else:
        cn_items = ["CN7", "CN1", "CN3"]
        cell_items = ["Vascular", "HER2+ Tumor", "PanCK+ Tumor", "CD8 T cell"]
    for item in cn_items:
        legend_cn.append(Line2D([0], [0], marker="o", color="w", label=f"{item} {CN_SHORT[item]}", markerfacecolor=CN_COLORS[item], markersize=5.0))
    for item in cell_items:
        legend_cn.append(Line2D([0], [0], marker="o", color="w", label=item, markerfacecolor=CELLTYPE_COLORS[item], markersize=5.0))
    legend_cn.append(Line2D([0], [0], marker="o", color="w", label="other cells", markerfacecolor=NEUTRAL_CELL, markersize=5.0))
    fig.legend(handles=legend_cn, loc="lower center", bbox_to_anchor=(0.5, -0.005), ncol=min(6, len(legend_cn)), fontsize=6.3)
    fig.subplots_adjust(left=0.055, right=0.995, top=0.90, bottom=0.15)
    source = pd.concat(source_parts, ignore_index=True) if source_parts else pd.DataFrame()
    export = save_figure(fig, basename, source, dpi=600)
    crop_table = pd.DataFrame(crop_records)
    crop_table.to_csv(PANEL_DIR / f"{basename}_crop_metadata.csv", index=False)
    return export, crop_table


def stat_row_from_table(table: pd.DataFrame, filters: dict[str, str]) -> pd.Series:
    mask = pd.Series(True, index=table.index)
    for col, value in filters.items():
        mask &= table[col].astype(str).eq(str(value))
    rows = table[mask]
    if rows.empty:
        raise ValueError(f"missing stat row: {filters}")
    return rows.iloc[0]


def plot_quantification(metrics: pd.DataFrame) -> FigureExport:
    basename = "fig6e_niche_quantification"
    wp13 = pd.read_csv(DELIVERABLE_DIR / "WP1.3_CRvsNR_stats.csv")
    wp14 = pd.read_csv(DELIVERABLE_DIR / "WP1.4_CRvsNR_stats.csv")
    wp15 = pd.read_csv(DELIVERABLE_DIR / "WP1.5_CRvsNR_stats.csv")
    specs = [
        ("cn5_prop", "CD21+B/TLS CN5 proportion", stat_row_from_table(wp13, {"analysis_set": PRIMARY_PROP_ANALYSIS, "feature_domain": "CN", "feature_id": "CN5"})),
        (
            "cd21b_same_type_knn_fraction",
            "CD21+B coherence",
            stat_row_from_table(
                wp15,
                {"analysis_set": "primary_threshold30k_panelaware", "feature_class": "coherence", "metric": "same_type_knn_fraction", "feature_id": "CD21+ B cell"},
            ),
        ),
        ("coupling_z_CN1_CN7", "CN1-CN7 coupling z", stat_row_from_table(wp14, {"analysis_set": PRIMARY_PROP_ANALYSIS, "feature_class": "adjacency", "metric": "z_score", "feature_id": "CN1__CN7"})),
        ("coupling_z_CN3_CN7", "CN3-CN7 coupling z", stat_row_from_table(wp14, {"analysis_set": PRIMARY_PROP_ANALYSIS, "feature_class": "adjacency", "metric": "z_score", "feature_id": "CN3__CN7"})),
    ]
    rows = []
    for key, label, stat in specs:
        for _, row in metrics.iterrows():
            rows.append(
                {
                    "sample": row["sample"],
                    "Outcome_short": row["Outcome_short"],
                    "metric_key": key,
                    "metric_label": label,
                    "value": row[key],
                    "cohens_d_NR_minus_CR": stat["cohens_d_NR_minus_CR"],
                    "cohens_d_bootstrap95_low": stat["cohens_d_bootstrap95_low"],
                    "cohens_d_bootstrap95_high": stat["cohens_d_bootstrap95_high"],
                    "wilcoxon_mannwhitney_p": stat["wilcoxon_mannwhitney_p"],
                    "FDR_BH_within_analysis_set": stat["FDR_BH_within_analysis_set"],
                    "direction": stat["direction"],
                }
            )
    source = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, len(specs), figsize=(10.4, 3.0))
    rng = np.random.default_rng(31)
    for ax, (key, label, stat) in zip(axes, specs):
        data = source[source["metric_key"].eq(key)]
        for xpos, group in enumerate(["CR", "NR"]):
            vals = data.loc[data["Outcome_short"].eq(group), "value"].dropna().astype(float).to_numpy()
            ax.scatter(np.full(len(vals), xpos) + rng.normal(0, 0.035, len(vals)), vals, color=GROUP_COLORS[group], s=15, edgecolor="white", linewidth=0.25, alpha=0.82)
            if len(vals):
                mean = float(vals.mean())
                boot = [float(np.mean(rng.choice(vals, len(vals), replace=True))) for _ in range(1000)]
                low, high = np.quantile(boot, [0.025, 0.975])
                ax.errorbar([xpos], [mean], yerr=[[mean - low], [high - mean]], fmt="_", color="black", ms=14, lw=1.0, capsize=2.5)
        ax.set_title(label, fontsize=7.8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["CR", "NR"])
        ax.grid(axis="y", color="#E7E7E7", linewidth=0.6)
        ax.text(
            0.5,
            0.98,
            f"d={float(stat['cohens_d_NR_minus_CR']):.2f}\n95%CI {float(stat['cohens_d_bootstrap95_low']):.2f},{float(stat['cohens_d_bootstrap95_high']):.2f}\np={float(stat['wilcoxon_mannwhitney_p']):.3g}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=6.2,
            bbox=dict(facecolor="white", edgecolor="#DDDDDD", alpha=0.88, pad=2.0),
        )
    axes[0].set_ylabel("Per-sample value")
    fig.text(0.01, 0.98, "e", fontsize=13, weight="bold", va="top")
    fig.subplots_adjust(left=0.06, right=0.995, top=0.80, bottom=0.18, wspace=0.34)
    return save_figure(fig, basename, source, dpi=600)


def reconstruct_wp22_folds() -> tuple[pd.DataFrame, pd.DataFrame]:
    script = DELIVERABLE_DIR / "WP2.2_fully_nested_prediction.py"
    spec = importlib.util.spec_from_file_location("wp22_prediction", script)
    if spec is None or spec.loader is None:
        raise ImportError(f"could not import {script}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    matrices, y, meta, _catalog = module.build_feature_matrices("primary_threshold30k")
    result = module.nested_loocv_topn_predictions(
        matrices["full_candidate_pool"],
        y,
        model_name="fold_composite",
        top_n=3,
        random_seed=20260708,
        min_per_class=5,
        inner_splits=3,
        selection_mode="data_driven_topn",
    )
    folds = result.fold_predictions.copy().merge(meta, on="sample", how="left", validate="one_to_one")
    return folds, result.selected_features


def roc_curve_with_bootstrap_ci(y_true: Iterable[int], scores: Iterable[float], *, n_bootstrap: int = 2000, random_seed: int = 7) -> dict[str, np.ndarray | float]:
    y = np.asarray(list(y_true), dtype=int)
    s = np.asarray(list(scores), dtype=float)
    fpr, tpr, _ = roc_curve(y, s)
    grid = np.linspace(0, 1, 151)
    obs = np.interp(grid, fpr, tpr)
    obs[0], obs[-1] = 0, 1
    rng = np.random.default_rng(random_seed)
    boots = []
    for _ in range(n_bootstrap):
        idx = rng.choice(len(y), len(y), replace=True)
        if len(np.unique(y[idx])) < 2:
            continue
        bfpr, btpr, _ = roc_curve(y[idx], s[idx])
        b = np.interp(grid, bfpr, btpr)
        b[0], b[-1] = 0, 1
        boots.append(b)
    arr = np.vstack(boots)
    return {"fpr": fpr, "tpr": tpr, "grid": grid, "mean": obs, "low": np.quantile(arr, 0.025, axis=0), "high": np.quantile(arr, 0.975, axis=0), "auc": float(roc_auc_score(y, s))}


def plot_prediction_roc() -> FigureExport:
    basename = "fig6f_prediction_roc"
    folds, selected = reconstruct_wp22_folds()
    summary = pd.read_csv(DELIVERABLE_DIR / "WP2.2_nested_prediction_summary.csv")
    row = summary[
        summary["analysis_set"].eq("primary_threshold30k")
        & summary["feature_pool"].eq("full_candidate_pool")
        & summary["selection_mode"].eq("data_driven_topn")
        & summary["top_n"].astype(str).eq("3")
        & summary["model"].eq("fold_composite")
    ].iloc[0]
    roc = roc_curve_with_bootstrap_ci(folds["truth"], folds["score"], n_bootstrap=2000, random_seed=20260708)
    curve_source = pd.DataFrame({"source_type": "roc_ci_grid", "fpr": roc["grid"], "tpr": roc["mean"], "tpr_low": roc["low"], "tpr_high": roc["high"]})
    fold_source = folds.copy()
    fold_source.insert(0, "source_type", "fold_prediction")
    source = pd.concat([fold_source, curve_source], ignore_index=True, sort=False)
    source["auc_wp22_summary"] = float(row["auc"])
    source["auc_bootstrap95_low_wp22_summary"] = float(row["auc_bootstrap95_low"])
    source["auc_bootstrap95_high_wp22_summary"] = float(row["auc_bootstrap95_high"])
    source["permutation_p_wp22_summary"] = float(row["permutation_p"])

    fig, ax = plt.subplots(figsize=(3.2, 3.05))
    ax.fill_between(roc["grid"], roc["low"], roc["high"], color="#A6BDD9", alpha=0.32, linewidth=0, label="Bootstrap 95% CI band")
    ax.plot(roc["grid"], roc["mean"], color="#2C5F9E", lw=1.7, label="Top-3 fold composite")
    ax.plot([0, 1], [0, 1], color="#999999", lw=0.8, linestyle="--", label="Chance")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.text(
        0.06,
        0.94,
        f"AUC {float(row['auc']):.3f}\n95% CI {float(row['auc_bootstrap95_low']):.2f}-{float(row['auc_bootstrap95_high']):.2f}\npermutation p={float(row['permutation_p']):.3f}\nwarrants validation",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.8,
        bbox=dict(facecolor="white", edgecolor="#DDDDDD", alpha=0.92, pad=2.2),
    )
    ax.legend(loc="lower right", fontsize=6.2)
    fig.text(0.02, 0.98, "f", fontsize=13, weight="bold", va="top")
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.16, top=0.97)
    export = save_figure(fig, basename, source, dpi=600)
    selected.to_csv(PANEL_DIR / f"{basename}_selected_features_source_data.csv", index=False)
    return export


def build_candidate_index(immune: pd.DataFrame, vascular: pd.DataFrame, crops: dict[str, pd.DataFrame], exports: dict[str, FigureExport]) -> pd.DataFrame:
    rows = []
    for panel_label, selected, export_key, crop_key in [
        ("6c immune/TLS crop", immune, "fig6c_spatial_immune_tls_crop_candidates", "6c"),
        ("6d vascular-axis crop", vascular, "fig6d_spatial_vascular_axis_crop_candidates", "6d"),
    ]:
        crop_table = crops[crop_key].set_index("sample")
        for _, row in selected.iterrows():
            sample = row["sample"]
            crop_row = crop_table.loc[sample].to_dict()
            rows.append(
                {
                    "panel": panel_label,
                    "sample": sample,
                    "Outcome_short": row["Outcome_short"],
                    "selection_rank": int(row["selection_rank"]),
                    "selection_distance": float(row["selection_distance"]),
                    "selection_basis": row["selection_basis"],
                    "selection_target_median": row["selection_target_median"],
                    "crop_xmin": crop_row.get("xmin"),
                    "crop_xmax": crop_row.get("xmax"),
                    "crop_ymin": crop_row.get("ymin"),
                    "crop_ymax": crop_row.get("ymax"),
                    "crop_size_um": (crop_row.get("xmax") - crop_row.get("xmin")) * PIXEL_SIZE_UM if pd.notna(crop_row.get("xmax")) else np.nan,
                    "cn5_prop": row.get("cn5_prop", np.nan),
                    "cd21b_same_type_knn_fraction": row.get("cd21b_same_type_knn_fraction", np.nan),
                    "coupling_z_CN1_CN7": row.get("coupling_z_CN1_CN7", np.nan),
                    "coupling_z_CN3_CN7": row.get("coupling_z_CN3_CN7", np.nan),
                    "coupling_z_CN4_CN7": row.get("coupling_z_CN4_CN7", np.nan),
                    "figure_svg": str(exports[export_key].files["svg"]),
                    "source_data": str(exports[export_key].files["source_data"]),
                }
            )
    return pd.DataFrame(rows)


def audit_exports(exports: dict[str, FigureExport]) -> pd.DataFrame:
    rows = []
    for export in exports.values():
        row = {"figure": export.basename, "source_rows": export.source_rows, "svg_bytes": export.svg_bytes}
        for fmt, path in export.files.items():
            row[f"{fmt}_path"] = str(path)
            row[f"{fmt}_bytes"] = int(path.stat().st_size) if path.exists() else 0
        if export.files["tiff"].exists():
            with Image.open(export.files["tiff"]) as img:
                row["tiff_width_px"] = int(img.size[0])
                row["tiff_height_px"] = int(img.size[1])
                dpi = img.info.get("dpi", (np.nan, np.nan))
                row["tiff_dpi_x"] = float(dpi[0]) if dpi else np.nan
                row["tiff_dpi_y"] = float(dpi[1]) if dpi else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def markdown_table(df: pd.DataFrame, columns: list[str], max_rows: int = 20) -> str:
    view = df.loc[:, columns].head(max_rows).copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "NA" if pd.isna(x) else f"{x:.3f}")
    view = view.fillna("NA").astype(str).map(lambda x: x.replace("|", "\\|"))
    lines = ["| " + " | ".join(view.columns) + " |", "| " + " | ".join(["---"] * len(view.columns)) + " |"]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(row.tolist()) + " |")
    return "\n".join(lines)


def build_report(index: pd.DataFrame, audit: pd.DataFrame, crops: dict[str, pd.DataFrame], exports: dict[str, FigureExport]) -> str:
    vascular = index[index["panel"].eq("6d vascular-axis crop")].copy()
    immune = index[index["panel"].eq("6c immune/TLS crop")].copy()
    spatial_audit = audit[audit["figure"].isin(["fig6c_spatial_immune_tls_crop_candidates", "fig6d_spatial_vascular_axis_crop_candidates"])].copy()
    report = f"""# WP3.1b Fig. 6 Revision Report

Status: completed.

## Revision Summary

- Output directory: `{PANEL_DIR}`.
- Manuscript palette applied: CN0-CN9 use locked `colors_15`; CR uses `#E4212B`; NR uses `#2C5F9E`; Pre/Post constants are recorded as `#E08A3C`/`#4A78B5`.
- 6c/6d are now cell-level crop panels with phenotype and CN rows, locator insets, 100 µm scale bars, and rasterized scatter layers with vector labels/boxes.
- 6a is a schematic placeholder for BioRender redraw. 6f retains the WP2.2 ROC logic and "warrants validation" annotation.

## Figure Contract

Core conclusion: Fig. 6 should show that CR is associated with immune/TLS organization and HER2-tumor/vascular niche coupling, whereas NR shows PanCK/HER2-low vascular coupling, with this vascular-coupling axis acting as a potential baseline predictor requiring validation.

Figure archetype: image plate + quant / asymmetric mixed-modality figure.

Backend: Python only for plotting, export, and visual QA.

## Spatial Candidate Selection

6c retains the WP3.1 objective immune/TLS representative rule: robust distance to group medians for CN5 proportion, CD21+B proportion, CD21+B KNN coherence, and CN5 patch size among baseline v2-like samples.

6d uses the revised headline rule: CR candidates are selected on the effect side closest to the CR median of CN1-CN7 coupling z; NR candidates are selected on the effect side closest to the NR median of CN3-CN7 coupling z. This avoids the previous composite-distance selection problem.

### 6c Candidates

{markdown_table(immune, ["sample", "Outcome_short", "selection_rank", "selection_distance", "cn5_prop", "cd21b_same_type_knn_fraction", "crop_size_um"], 10)}

### 6d Candidates

{markdown_table(vascular, ["sample", "Outcome_short", "selection_rank", "selection_distance", "coupling_z_CN1_CN7", "coupling_z_CN3_CN7", "coupling_z_CN4_CN7", "crop_size_um"], 10)}

## SVG Size Audit

Target for spatial SVGs was approximately <1.5 MB. Rasterized scatter artists reduced spatial SVGs to the sizes below while keeping text, scale bars, legends, and locator boxes vector.

{markdown_table(spatial_audit, ["figure", "svg_bytes", "pdf_bytes", "png_bytes", "tiff_dpi_x", "source_rows"], 10)}

## Output Manifest

{markdown_table(audit, ["figure", "svg_path", "pdf_path", "png_path", "tiff_path", "source_rows", "svg_bytes"], 20)}

## Source Data

Every panel has a source-data CSV next to the figure files. Spatial source data include crop bounds and the cells plotted in each crop. Additional crop metadata tables are written as:

- `{PANEL_DIR / "fig6c_spatial_immune_tls_crop_candidates_crop_metadata.csv"}`
- `{PANEL_DIR / "fig6d_spatial_vascular_axis_crop_candidates_crop_metadata.csv"}`

## Notes

- Spatial wording should remain "niche/territory coupling", not direct single-cell vascular contact.
- Raw HER2/CD31 crops remain outside this WP because the revision requested registered B0 phenotype/CN crops and did not specify objective raw-image crop coordinates.
- Candidate index: `{OUT_INDEX}`.
"""
    return report


def main() -> None:
    configure_matplotlib()
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    source_script = Path(__file__).resolve()
    if source_script != OUT_SCRIPT:
        OUT_SCRIPT.write_text(source_script.read_text(encoding="utf-8"), encoding="utf-8")

    niche_map = load_niche_map()
    metrics = build_sample_metrics()
    spatial_pool = metrics[
        metrics["Time_short"].eq("Pre")
        & metrics["Outcome_short"].isin(["CR", "NR"])
        & metrics["passes_lowcell_30000"].astype(bool)
        & (~metrics["is_v1_like"].astype(bool))
        & (~metrics["sample"].isin(NO_SMA_SAMPLES))
    ].copy()
    immune_candidates = select_representative_samples(
        spatial_pool,
        metric_cols=["cn5_prop", "cd21b_prop", "cd21b_same_type_knn_fraction", "cn5_mean_patch_size"],
        per_group=3,
    )
    vascular_candidates = select_headline_effect_side_samples(
        spatial_pool,
        group_specs={
            "CR": {"metric": "coupling_z_CN1_CN7", "direction": 1},
            "NR": {"metric": "coupling_z_CN3_CN7", "direction": 1},
        },
        per_group=3,
    )
    all_samples = sorted(set(immune_candidates["sample"]).union(vascular_candidates["sample"]))
    spatial = load_spatial_for_samples(all_samples)

    exports: dict[str, FigureExport] = {}
    crop_tables: dict[str, pd.DataFrame] = {}
    for export in [
        plot_workflow_placeholder(),
        plot_composition_heatmap(niche_map),
        plot_quantification(metrics),
        plot_prediction_roc(),
    ]:
        exports[export.basename] = export

    export_6c, crops_6c = plot_spatial_crop_panel(
        basename="fig6c_spatial_immune_tls_crop_candidates",
        letter="c",
        title="Immune/TLS organization candidates (cell-level crop)",
        panel_id="6c",
        selected=immune_candidates,
        spatial=spatial,
        metrics=metrics,
    )
    exports[export_6c.basename] = export_6c
    crop_tables["6c"] = crops_6c

    export_6d, crops_6d = plot_spatial_crop_panel(
        basename="fig6d_spatial_vascular_axis_crop_candidates",
        letter="d",
        title="Vascular niche-coupling candidates (cell-level crop)",
        panel_id="6d",
        selected=vascular_candidates,
        spatial=spatial,
        metrics=metrics,
    )
    exports[export_6d.basename] = export_6d
    crop_tables["6d"] = crops_6d

    candidate_index = build_candidate_index(immune_candidates, vascular_candidates, crop_tables, exports)
    candidate_index.to_csv(OUT_INDEX, index=False)
    audit = audit_exports(exports)
    audit.to_csv(PANEL_DIR / "WP3.1b_export_audit.csv", index=False)
    report = build_report(candidate_index, audit, crop_tables, exports)
    OUT_REPORT.write_text(report, encoding="utf-8")
    metadata = {
        "status": "completed",
        "panel_dir": str(PANEL_DIR),
        "candidate_index": str(OUT_INDEX),
        "report": str(OUT_REPORT),
        "palette": {"colors_15": COLORS_15, "CR": GROUP_COLORS["CR"], "NR": GROUP_COLORS["NR"], "Pre": TIME_COLORS["Pre"], "Post": TIME_COLORS["Post"]},
        "exports": {
            name: {fmt: str(path) for fmt, path in export.files.items()} | {"source_rows": export.source_rows, "svg_bytes": export.svg_bytes}
            for name, export in exports.items()
        },
        "immune_candidates": immune_candidates["sample"].tolist(),
        "vascular_candidates": vascular_candidates["sample"].tolist(),
        "scale_bar_um": SCALE_BAR_UM,
        "pixel_size_um": PIXEL_SIZE_UM,
        "spatial_svg_target_bytes": 1_500_000,
    }
    OUT_METADATA.write_text(json.dumps(metadata, indent=2, ensure_ascii=False, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
