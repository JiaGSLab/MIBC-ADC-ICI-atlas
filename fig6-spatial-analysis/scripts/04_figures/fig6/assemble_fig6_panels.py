from __future__ import annotations

import importlib.util
import json
import math
import os
import shutil
import sys
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve


PROJECT_DIR = Path(os.environ.get("FIG6_PROJECT_DIR", "/data/home/chenwenhui/mIHC/fig6_pm"))
DELIVERABLE_DIR = PROJECT_DIR / "deliverables"
CELLINFO_DIR = Path(os.environ.get("FIG6_CELLINFO_DIR", "/data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh"))
PANEL_DIR = CELLINFO_DIR / "fig6_panels_rev"
SPATIAL_PARQUET = CELLINFO_DIR / "neighborhood_radius200_B0" / "neighborhood_composition.parquet"
SCRIPT_NAME = "assemble_fig6_panels.py"

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
CELLTYPE_ORDER = [
    "B cell",
    "CD163+ Macrophage",
    "CD163- Macrophage",
    "CD21+ B cell",
    "CD4 T cell",
    "CD8 T cell",
    "CD8 exhausted T cell",
    "DC",
    "Vascular",
    "Fibroblast / SMA+ stromal",
    "HER2+ Tumor",
    "Treg",
    "PanCK+ Tumor",
]
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
PIXEL_SIZE_UM = 0.4995            # µm / px,与 WP3_1b 一致
PIXELS_PER_UM = 1.0 / PIXEL_SIZE_UM
SCALE_BAR_UM = 100
SCALE_BAR_PX = SCALE_BAR_UM / PIXEL_SIZE_UM

CONFIG = {
    "project_dir": str(PROJECT_DIR),
    "script_dir": str(Path(__file__).resolve().parent),
    "cellinfo_dir": str(CELLINFO_DIR),
    "figure_size": [16.0, 12.8],
    "preview_dpi": 220,
    "export_dpi": 600,
    "background": "white",
    "font_family": "Arial",
    "gridspec": {
        "height_ratios": [1.25, 2.25, 2.25, 1.85],
        "width_ratios": [1, 1, 1, 1, 1, 1],
        "hspace": 0.34,
        "wspace": 0.25,
        "margins": {"left": 0.045, "right": 0.985, "top": 0.975, "bottom": 0.055},
    },
    "panel_6c_samples": ["pt42-pre", "pt63-pre", "pt50-pre", "pt48-pre"],
    "panel_6d_samples": ["pt60-pre", "pt58-pre", "pt55-pre", "pt22-pre"],
    "crop_size_px": {"6c": 1250, "6d": 3000},
    "crop_min_cells": {"6c": 350, "6d": 250},
    "crop_candidate_n": 3,
    "crop_overrides": {},
    # phenotype 行是否也按 CN 判定放大(联动高亮);False = 各图按各自的列
    "link_emphasis_by_cn": True,
    "crop_point_style": {
        "point_size": 1.3,          # 基准点径 (s)
        "other_keep_color": True,   # 背景点保留类别原色;False = 压成中性灰
        "highlight_size_scale": 2.0,  # 重点 CN/celltype 相对基准的放大倍数
        "other_size_scale": 0.8,      # 非重点细胞相对基准的缩小倍数
        "highlight_alpha": 0.95,
        "other_alpha": 0.28,
    },
    "sample_counts": {
        "experimental_patients": 48,
        "experimental_sections": 82,
        "experimental_pre": 43,
        "experimental_post": 39,
        "downstream_patients": 45,
        "downstream_sections": 72,
        "downstream_cells_million": 23.1,
        "baseline_pre_sections": 40,
        "baseline_pre_CR": 22,
        "baseline_pre_NR": 18,
        "prediction_threshold30k_sections": 36,
        "prediction_threshold30k_CR": 19,
        "prediction_threshold30k_NR": 17,
    },
}


def configure_matplotlib(config: dict | None = None) -> None:
    cfg = config or CONFIG
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [cfg.get("font_family", "Arial"), "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
            "font.size": 7.2,
            "axes.linewidth": 0.55,
            "axes.edgecolor": "#333333",
            "axes.labelsize": 7.4,
            "axes.titlesize": 8.2,
            "xtick.labelsize": 6.6,
            "ytick.labelsize": 6.6,
            "legend.fontsize": 6.3,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": cfg.get("background", "white"),
            "axes.facecolor": cfg.get("background", "white"),
        }
    )


def _wp31b_module():
    candidates = [
        Path(__file__).with_name("build_fig6b_c_e_panel_sources.py"),
        DELIVERABLE_DIR / "WP3_1b_fig6_revision.py",
        PROJECT_DIR / "scratch_wp3.1b" / "WP3_1b_fig6_revision.py",
    ]
    for path in candidates:
        if path.exists():
            spec = importlib.util.spec_from_file_location("wp31b_revision", path)
            if spec is None or spec.loader is None:
                continue
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            return module
    raise FileNotFoundError("Could not locate build_fig6b_c_e_panel_sources.py or the original WP3_1b helper")


def ensure_dirs() -> None:
    DELIVERABLE_DIR.mkdir(parents=True, exist_ok=True)
    PANEL_DIR.mkdir(parents=True, exist_ok=True)


def save_figure(fig: plt.Figure, basename: str, source: pd.DataFrame, *, dpi: int = 600, panel_dir: Path = PANEL_DIR) -> dict[str, str]:
    panel_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "svg": panel_dir / f"{basename}.svg",
        "pdf": panel_dir / f"{basename}.pdf",
        "png": panel_dir / f"{basename}.png",
        "tiff": panel_dir / f"{basename}.tiff",
        "source_data": panel_dir / f"{basename}_source_data.csv",
    }
    source.to_csv(paths["source_data"], index=False)
    for fmt in ["svg", "pdf", "png", "tiff"]:
        fig.savefig(paths[fmt], dpi=dpi, transparent=False, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return {key: str(value) for key, value in paths.items()}


def load_primary_metrics() -> pd.DataFrame:
    return _wp31b_module().build_sample_metrics()


def load_niche_map() -> pd.DataFrame:
    return pd.read_csv(DELIVERABLE_DIR / "WP1.2_functional_niche_map_k10.csv")


def load_spatial_for_samples(samples: list[str]) -> pd.DataFrame:
    return _wp31b_module().load_spatial_for_samples(samples)


def legend_payloads() -> dict[str, list[str]]:
    return {
        "phenotype": ["CD21+ B cell", "CD8 T cell", "HER2+ Tumor", "PanCK+ Tumor", "Vascular"],
        "CN": ["CN1", "CN3", "CN4", "CN5", "CN6", "CN7"],
    }


def extract_auc_ladder(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []

    def best_row(feature_pool: str, label: str, tier: str) -> pd.Series:
        subset = summary[
            summary["analysis_set"].astype(str).eq("primary_threshold30k")
            & summary["feature_pool"].astype(str).eq(feature_pool)
            & summary["selection_mode"].astype(str).eq("fixed_all")
            & summary["top_n"].astype(str).eq("all")
        ].copy()
        if subset.empty:
            raise ValueError(f"Missing WP2.2 summary rows for {feature_pool}")
        row = subset.sort_values("auc", ascending=False).iloc[0].copy()
        row["display_label"] = label
        row["tier"] = tier
        return row

    rows.append(best_row("celltype_proportion_only", "Celltype proportions", "celltype"))
    rows.append(best_row("cn_proportion_only", "CN proportions", "cn_prop"))
    coupling = summary[
        summary["analysis_set"].astype(str).eq("primary_threshold30k")
        & summary["feature_pool"].astype(str).eq("full_candidate_pool")
        & summary["selection_mode"].astype(str).eq("data_driven_topn")
        & summary["top_n"].astype(str).eq("3")
        & summary["model"].astype(str).eq("fold_composite")
    ].copy()
    if coupling.empty:
        raise ValueError("Missing WP2.2 top-3 fold_composite row")
    row = coupling.iloc[0].copy()
    row["display_label"] = "CN vascular coupling"
    row["tier"] = "cn_coupling"
    rows.append(row)

    out = pd.DataFrame(rows).reset_index(drop=True)
    out.insert(0, "rank", [1, 2, 3])
    keep = [
        "rank",
        "display_label",
        "tier",
        "auc",
        "auc_bootstrap95_low",
        "auc_bootstrap95_high",
        "permutation_p",
        "model",
        "n_samples",
        "n_CR",
        "n_NR",
        "decision_tier",
    ]
    return out[[c for c in keep if c in out.columns]]


def compute_tls_presence(
    prop: pd.DataFrame,
    zero_inflation: pd.DataFrame,
    *,
    threshold: float = 0.01,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = prop[
        prop["Time_short"].astype(str).eq("Pre")
        & prop["Outcome_short"].astype(str).isin(["CR", "NR"])
        & prop["passes_lowcell_30000"].astype(bool)
        & prop["feature_domain"].astype(str).eq("CN")
        & prop["k"].astype(str).eq("10")
        & prop["feature_id"].astype(str).eq("CN5")
    ].copy()
    rows = rows.reset_index(drop=True)
    rows["tls_present"] = rows["proportion"].astype(float).ge(threshold)
    rows["presence_threshold"] = threshold
    summary = (
        rows.groupby("Outcome_short")["tls_present"]
        .agg(n="size", present="sum")
        .reindex(["CR", "NR"])
        .fillna(0)
        .astype(int)
    )
    summary["absent"] = summary["n"] - summary["present"]
    summary["fraction_present"] = np.where(summary["n"] > 0, summary["present"] / summary["n"], np.nan)

    zrow = zero_inflation[
        zero_inflation["feature_domain"].astype(str).eq("CN")
        & zero_inflation["k"].astype(str).eq("10")
        & zero_inflation["feature_id"].astype(str).eq("CN5")
        & zero_inflation["presence_threshold_name"].astype(str).eq("ge1pct")
    ]
    if not zrow.empty:
        z = zrow.iloc[0]
        summary.attrs["fisher_p"] = float(z["presence_fisher_p"])
        summary.attrs["source_CR_present"] = int(z["CR_present"])
        summary.attrs["source_NR_present"] = int(z["NR_present"])
    else:
        table = [[summary.loc["CR", "present"], summary.loc["CR", "absent"]], [summary.loc["NR", "present"], summary.loc["NR", "absent"]]]
        summary.attrs["fisher_p"] = float(stats.fisher_exact(table)[1])
    return rows, summary


def panel_crop_targets(panel: str, outcome: str) -> dict[str, object]:
    if panel == "6c":
        if outcome == "CR":
            return {
                "cn_targets": ["CN5", "CN4", "CN6"],
                "cell_targets": ["CD21+ B cell", "CD8 T cell", "HER2+ Tumor", "PanCK+ Tumor"],
                "weights": {"CN5": 8.0, "CN4": 2.2, "CN6": 0.9, "CD21+ B cell": 5.0, "CD8 T cell": 1.8},
                "interface_weight": 8.0,
            }
        return {
            "cn_targets": ["CN6", "CN3", "CN5", "CN4"],
            "cell_targets": ["HER2+ Tumor", "PanCK+ Tumor", "CD21+ B cell", "CD8 T cell"],
            "weights": {"CN6": 3.0, "CN3": 3.0, "CN5": 2.0, "CD21+ B cell": 2.5},
            "interface_weight": 4.0,
        }
    if outcome == "CR":
        return {
            "cn_targets": ["CN7", "CN1"],
            "cell_targets": ["Vascular", "HER2+ Tumor", "CD8 T cell"],
            "weights": {"CN7": 5.0, "CN1": 5.0, "Vascular": 4.0, "HER2+ Tumor": 2.0},
            "interface_weight": 50.0,
        }
    return {
        "cn_targets": ["CN7", "CN3"],
        "cell_targets": ["Vascular", "PanCK+ Tumor", "HER2+ Tumor"],
        "weights": {"CN7": 5.0, "CN3": 5.0, "Vascular": 4.0, "PanCK+ Tumor": 2.0},
        "interface_weight": 50.0,
    }


def _candidate_centers(df: pd.DataFrame, key: pd.DataFrame, half: float, grid_n: int = 10) -> list[tuple[float, float]]:
    x_min, x_max = float(df["spatial_x"].min()), float(df["spatial_x"].max())
    y_min, y_max = float(df["spatial_y"].min()), float(df["spatial_y"].max())
    lo_x, hi_x = x_min + half, x_max - half
    lo_y, hi_y = y_min + half, y_max - half
    if hi_x < lo_x:
        lo_x = hi_x = float(df["spatial_x"].median())
    if hi_y < lo_y:
        lo_y = hi_y = float(df["spatial_y"].median())
    xs = np.linspace(lo_x, hi_x, min(grid_n, max(2, grid_n)))
    ys = np.linspace(lo_y, hi_y, min(grid_n, max(2, grid_n)))
    centers = [(float(x), float(y)) for x in xs for y in ys]
    if not key.empty:
        for qx in [0.15, 0.35, 0.50, 0.65, 0.85]:
            for qy in [0.15, 0.35, 0.50, 0.65, 0.85]:
                centers.append((float(key["spatial_x"].quantile(qx)), float(key["spatial_y"].quantile(qy))))
        centers.append((float(key["spatial_x"].median()), float(key["spatial_y"].median())))
    clipped = []
    seen = set()
    for cx, cy in centers:
        cx = min(max(cx, lo_x), hi_x)
        cy = min(max(cy, lo_y), hi_y)
        token = (round(cx, 1), round(cy, 1))
        if token not in seen:
            clipped.append((cx, cy))
            seen.add(token)
    return clipped


def smart_crop_candidates(
    df: pd.DataFrame,
    *,
    panel: str,
    outcome: str,
    sample: str,
    crop_size: float,
    n_candidates: int = 3,
    min_cells: int = 300,
    manual_override: dict | None = None,
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    spec = panel_crop_targets(panel, outcome)
    cn_targets = list(spec["cn_targets"])
    cell_targets = list(spec["cell_targets"])
    weights = dict(spec.get("weights", {}))
    interface_weight = float(spec.get("interface_weight", 4.0))
    half = crop_size / 2.0
    key = df[df["CN_k10"].astype(str).isin(cn_targets) | df["celltype"].astype(str).isin(cell_targets)].copy()
    centers = _candidate_centers(df, key, half)
    if manual_override:
        if {"center_x", "center_y"}.issubset(manual_override):
            centers.insert(0, (float(manual_override["center_x"]), float(manual_override["center_y"])))
        elif {"xmin", "xmax", "ymin", "ymax"}.issubset(manual_override):
            centers.insert(
                0,
                (
                    (float(manual_override["xmin"]) + float(manual_override["xmax"])) / 2.0,
                    (float(manual_override["ymin"]) + float(manual_override["ymax"])) / 2.0,
                ),
            )
    records = []
    x = df["spatial_x"].to_numpy(dtype=float)
    y = df["spatial_y"].to_numpy(dtype=float)
    cn = df["CN_k10"].astype(str).to_numpy()
    cell = df["celltype"].astype(str).to_numpy()
    for idx, (cx, cy) in enumerate(centers):
        mask = (x >= cx - half) & (x <= cx + half) & (y >= cy - half) & (y <= cy + half)
        n = int(mask.sum())
        if n == 0:
            continue
        cn_counts = {target: int(np.sum(mask & (cn == target))) for target in cn_targets}
        cell_counts = {target: int(np.sum(mask & (cell == target))) for target in cell_targets}
        target_score = 0.0
        for target, count in {**cn_counts, **cell_counts}.items():
            target_score += weights.get(target, 1.0) * min(count, 2500)
        interface = 0
        if len(cn_targets) >= 2:
            interface = min(cn_counts.get(cn_targets[0], 0), cn_counts.get(cn_targets[1], 0))
        density_score = min(n / max(min_cells, 1), 4.0) * 150.0
        score = target_score + interface_weight * min(interface, 2500) + density_score
        if panel == "6d" and interface == 0:
            score *= 0.02
        if panel == "6c" and outcome == "CR":
            primary_count = cn_counts.get(cn_targets[0], 0)
            if primary_count == 0:
                score *= 0.12
            else:
                score += 18.0 * min(primary_count, 1500)
        if n < min_cells:
            score *= max(n / max(min_cells, 1), 0.05)
        records.append(
            {
                "panel": panel,
                "sample": sample,
                "Outcome_short": outcome,
                "center_x": float(cx),
                "center_y": float(cy),
                "xmin": float(cx - half),
                "xmax": float(cx + half),
                "ymin": float(cy - half),
                "ymax": float(cy + half),
                "n_cells": n,
                "density_pass": bool(n >= min_cells),
                "target_score": float(target_score),
                "interface_count": int(interface),
                "score": float(score),
                "manual_override": bool(manual_override and idx == 0),
                **{f"n_{k}": v for k, v in cn_counts.items()},
                **{f"n_{k}": v for k, v in cell_counts.items()},
            }
        )
    if not records:
        cx, cy = float(df["spatial_x"].median()), float(df["spatial_y"].median())
        records.append(
            {
                "panel": panel,
                "sample": sample,
                "Outcome_short": outcome,
                "center_x": cx,
                "center_y": cy,
                "xmin": cx - half,
                "xmax": cx + half,
                "ymin": cy - half,
                "ymax": cy + half,
                "n_cells": 0,
                "density_pass": False,
                "target_score": 0.0,
                "interface_count": 0,
                "score": 0.0,
                "manual_override": False,
            }
        )
    ranked = pd.DataFrame(records).sort_values(["manual_override", "density_pass", "score", "n_cells"], ascending=[False, False, False, False])
    selected = []
    for _, row in ranked.iterrows():
        too_close = any(math.hypot(row["center_x"] - r["center_x"], row["center_y"] - r["center_y"]) < crop_size * 0.35 for r in selected)
        if not too_close or len(selected) == 0:
            selected.append(row.to_dict())
        if len(selected) >= n_candidates:
            break
    if len(selected) < min(n_candidates, len(ranked)):
        for _, row in ranked.iterrows():
            if len(selected) >= n_candidates:
                break
            if not any(abs(row["center_x"] - r["center_x"]) < 1 and abs(row["center_y"] - r["center_y"]) < 1 for r in selected):
                selected.append(row.to_dict())
    out = pd.DataFrame(selected).reset_index(drop=True)
    out.insert(0, "candidate_rank", np.arange(1, len(out) + 1))
    return out


def crop_dataframe(df: pd.DataFrame, crop: pd.Series | dict) -> pd.DataFrame:
    return df[
        (df["spatial_x"] >= float(crop["xmin"]))
        & (df["spatial_x"] <= float(crop["xmax"]))
        & (df["spatial_y"] >= float(crop["ymin"]))
        & (df["spatial_y"] <= float(crop["ymax"]))
    ].copy()


def scatter_crop(
    ax: plt.Axes,
    crop_df: pd.DataFrame,
    *,
    color_col: str,
    highlights: Iterable[str] | None = None,
    palette: dict[str, str],
    crop: pd.Series | dict,
    point_size: float,
    style: dict | None = None,
    emphasis_col: str | None = None,
    emphasis_values: Iterable[str] | None = None,
    other_keep_color: bool | None = None,
) -> None:
    """重点细胞画大画实,其余画小画淡。

    emphasis_col / emphasis_values 决定"谁被放大"的依据列(默认按 color_col + highlights)。
    传 emphasis_col="CN_k10" 则 phenotype 图也按 CN 判定,两张图强调同一批细胞。
    other_keep_color 控制背景点是保留类别原色还是压成中性灰。
    """
    st = {**CONFIG["crop_point_style"], **(style or {})}
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
    ax.set_xlim(float(crop["xmin"]), float(crop["xmax"]))
    ax.set_ylim(float(crop["ymax"]), float(crop["ymin"]))
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#C9C9C9")
        spine.set_linewidth(0.42)


def draw_scale_bar(ax: plt.Axes, crop: pd.Series | dict) -> None:
    x0 = float(crop["xmax"]) - SCALE_BAR_PX - 32
    y0 = float(crop["ymax"]) - 34
    ax.plot([x0, x0 + SCALE_BAR_PX], [y0, y0], color="black", lw=1.2, solid_capstyle="butt", zorder=20)
    ax.text(x0 + SCALE_BAR_PX / 2, y0 - 16, f"{SCALE_BAR_UM} µm", ha="center", va="top", fontsize=5.9, color="black")


def add_locator(ax: plt.Axes, sample_df: pd.DataFrame, crop: pd.Series | dict, color: str) -> None:
    inset = inset_axes(ax, width="23%", height="23%", loc="upper right", borderpad=0.36)
    small = sample_df.sample(n=min(len(sample_df), 9000), random_state=17) if len(sample_df) > 9000 else sample_df
    inset.scatter(small["spatial_x"], small["spatial_y"], s=0.06, color="#CFCFCF", linewidths=0, rasterized=True)
    inset.add_patch(
        Rectangle(
            (float(crop["xmin"]), float(crop["ymin"])),
            float(crop["xmax"]) - float(crop["xmin"]),
            float(crop["ymax"]) - float(crop["ymin"]),
            fill=False,
            edgecolor=color,
            linewidth=0.75,
        )
    )
    inset.set_aspect("equal")
    inset.invert_yaxis()
    inset.set_xticks([])
    inset.set_yticks([])
    for spine in inset.spines.values():
        spine.set_color(color)
        spine.set_linewidth(0.45)


def draw_workflow_placeholder(ax: plt.Axes, config: dict) -> pd.DataFrame:
    counts = config["sample_counts"]
    ax.set_axis_off()
    ax.text(0.0, 1.0, "a", fontsize=13, weight="bold", transform=ax.transAxes, va="top")
    ax.text(0.50, 0.73, "Fig. 6 workflow placeholder", ha="center", va="center", fontsize=10.6, weight="bold", transform=ax.transAxes)
    ax.text(
        0.50,
        0.47,
        "Phenotypes -> k10 CN territories\n-> vascular coupling -> baseline prediction",
        ha="center",
        va="center",
        fontsize=7.0,
        color="#333333",
        transform=ax.transAxes,
        linespacing=1.35,
        clip_on=True,
    )
    ax.text(
        0.50,
        0.18,
        f"82 mIHC sections; downstream 72 sections / {counts['downstream_cells_million']:.1f}M cells; prediction n=36",
        ha="center",
        va="center",
        fontsize=6.8,
        color="#555555",
        transform=ax.transAxes,
    )
    return pd.DataFrame(
        [
            {
                "panel": "6a",
                "note": "BioRender placeholder",
                "experimental_sections": counts["experimental_sections"],
                "downstream_sections": counts["downstream_sections"],
                "downstream_cells_million": counts["downstream_cells_million"],
                "prediction_threshold30k_sections": counts["prediction_threshold30k_sections"],
            }
        ]
    )


def draw_composition_heatmap(ax: plt.Axes, cax: plt.Axes | None = None) -> pd.DataFrame:
    z = pd.read_csv(DELIVERABLE_DIR / "WP1.2_cn_composition_current_k10_heatmap_zscore_matrix.csv")
    z = z.rename(columns={z.columns[0]: "cell_phenotype"}).set_index("cell_phenotype").loc[:, CN_ORDER]
    niche = load_niche_map()[["CN", "functional_niche", "niche_family", "dominant_celltype"]]
    sns.heatmap(
        z,
        ax=ax,
        cmap=MANUSCRIPT_BWR,
        norm=TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3),
        cbar=cax is not None,
        cbar_ax=cax,
        cbar_kws={"label": "Z-score", "ticks": [-2, 0, 3]} if cax is not None else None,
        linewidths=0.35,
        linecolor="white",
        xticklabels=[f"{cn}\n{CN_SHORT[cn]}" for cn in CN_ORDER],
    )
    ax.set_xlabel("Cellular neighborhoods (k=10)")
    ax.set_ylabel("Cell phenotypes")
    ax.tick_params(axis="x", rotation=45)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    ax.text(-0.13, 1.06, "b", transform=ax.transAxes, fontsize=13, weight="bold", va="top")
    return z.reset_index().melt(id_vars="cell_phenotype", var_name="CN", value_name="z_score").merge(niche, on="CN", how="left")


def draw_tls_presence_panel(
    fig: plt.Figure,
    subspec,
    prop_rows: pd.DataFrame,
    summary: pd.DataFrame,
    stats13: pd.DataFrame,
    stats14: pd.DataFrame,
    *,
    compact: bool = False,
) -> pd.DataFrame:
    sub = subspec.subgridspec(1, 3, width_ratios=[1.05, 0.80, 1.25], wspace=0.42)
    ax0 = fig.add_subplot(sub[0, 0])
    ax1 = fig.add_subplot(sub[0, 1])
    ax2 = fig.add_subplot(sub[0, 2])
    ax0.text(-0.24, 1.06, "e", transform=ax0.transAxes, fontsize=13, weight="bold", va="top")
    rng = np.random.default_rng(20260708)
    for i, outcome in enumerate(["CR", "NR"]):
        vals = prop_rows[prop_rows["Outcome_short"].eq(outcome)]["proportion"].astype(float).to_numpy()
        jitter = rng.normal(0, 0.045, len(vals))
        ax0.scatter(np.full(len(vals), i) + jitter, vals, s=14, color=GROUP_COLORS[outcome], alpha=0.82, linewidths=0, zorder=3)
        if len(vals):
            mean = float(np.mean(vals))
            ci = stats.bootstrap((vals,), np.mean, confidence_level=0.95, n_resamples=2000, random_state=7).confidence_interval if len(vals) > 1 else None
            low = float(ci.low) if ci is not None else mean
            high = float(ci.high) if ci is not None else mean
            ax0.plot([i - 0.22, i + 0.22], [mean, mean], color="black", lw=1.0)
            ax0.vlines(i, low, high, color="black", lw=0.8)
    ax0.axhline(0.01, color="#777777", lw=0.8, ls="--")
    ax0.set_xticks([0, 1])
    ax0.set_xticklabels(["CR", "NR"])
    ax0.set_ylabel("CN5 proportion")
    ax0.set_title("CD21+B/TLS niche abundance")
    row = stats13[
        stats13["analysis_set"].astype(str).eq("primary_k10_threshold30k_panelaware")
        & stats13["feature_domain"].astype(str).eq("CN")
        & stats13["k"].astype(str).eq("10")
        & stats13["feature_id"].astype(str).eq("CN5")
    ].iloc[0]
    ax0.text(
        0.02,
        0.98,
        f"d={float(row['cohens_d_NR_minus_CR']):.2f}\n95% CI {float(row['cohens_d_bootstrap95_low']):.2f} to {float(row['cohens_d_bootstrap95_high']):.2f}",
        transform=ax0.transAxes,
        ha="left",
        va="top",
        fontsize=6.1,
        bbox=dict(facecolor="white", edgecolor="#DDDDDD", pad=1.6),
    )
    max_y = max(0.02, prop_rows["proportion"].max() * 1.12)
    ax0.set_ylim(0, max_y)

    p = summary.attrs.get("fisher_p", np.nan)
    for i, outcome in enumerate(["CR", "NR"]):
        frac = float(summary.loc[outcome, "fraction_present"])
        ax1.bar(i, frac, color=GROUP_COLORS[outcome], width=0.62)
        ax1.bar(i, 1 - frac, bottom=frac, color="#E5E5E5", width=0.62)
        ax1.text(i, min(frac + 0.05, 0.96), f"{int(summary.loc[outcome, 'present'])}/{int(summary.loc[outcome, 'n'])}", ha="center", va="bottom", fontsize=7.0)
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(["CR", "NR"])
    ax1.set_ylim(0, 1)
    ax1.set_ylabel("Presence fraction")
    ax1.set_title("TLS presence >=1%")
    ax1.text(0.50, 0.04, f"Fisher p={p:.3f}", transform=ax1.transAxes, ha="center", va="bottom", fontsize=6.3)

    effect_ids = ["CN1__CN7", "CN3__CN7", "CN4__CN7"]
    effect_rows = stats14[
        stats14["analysis_set"].astype(str).eq("primary_k10_threshold30k_panelaware")
        & stats14["feature_class"].astype(str).eq("adjacency")
        & stats14["metric"].astype(str).eq("z_score")
        & stats14["k"].astype(str).eq("10")
        & stats14["feature_id"].astype(str).isin(effect_ids)
    ].copy()
    if compact:
        label_map = {"CN1__CN7": "CN1-CN7", "CN3__CN7": "CN3-CN7", "CN4__CN7": "CN4-CN7"}
    else:
        label_map = {
            "CN1__CN7": "HER2 tumor-Mac | vascular",
            "CN3__CN7": "PanCK tumor | vascular",
            "CN4__CN7": "CD8 cytotoxic | vascular",
        }
    effect_rows["label"] = effect_rows["feature_id"].map(label_map)
    effect_rows["plot_order"] = effect_rows["feature_id"].map({v: i for i, v in enumerate(effect_ids)})
    effect_rows = effect_rows.sort_values("plot_order")
    ypos = np.arange(len(effect_rows))
    d = effect_rows["cohens_d_NR_minus_CR"].astype(float).to_numpy()
    low = effect_rows["cohens_d_bootstrap95_low"].astype(float).to_numpy()
    high = effect_rows["cohens_d_bootstrap95_high"].astype(float).to_numpy()
    ax2.axvline(0, color="#777777", lw=0.8, ls="--")
    ax2.hlines(ypos, low, high, color="#333333", lw=1.0)
    ax2.scatter(d, ypos, s=24, color=[CN_COLORS["CN1"], CN_COLORS["CN3"], CN_COLORS["CN4"]], zorder=3)
    ax2.set_yticks(ypos)
    ax2.set_yticklabels(effect_rows["label"].tolist())
    ax2.tick_params(axis="y", labelsize=5.8 if compact else 6.4, pad=1.0)
    ax2.invert_yaxis()
    ax2.set_xlabel("Cohen d (NR - CR)")
    ax2.set_title("Vascular coupling z")
    ax2.text(0.00, -0.22, "CR higher", transform=ax2.transAxes, ha="left", va="top", fontsize=6.0)
    ax2.text(1.00, -0.22, "NR higher", transform=ax2.transAxes, ha="right", va="top", fontsize=6.0)

    source_prop = prop_rows.copy()
    source_prop.insert(0, "source_type", "CN5_sample_proportion")
    source_summary = summary.reset_index().rename(columns={"Outcome_short": "Outcome_short"})
    source_summary.insert(0, "source_type", "TLS_presence_summary")
    source_summary["presence_fisher_p"] = p
    source_effect = effect_rows.copy()
    source_effect.insert(0, "source_type", "vascular_coupling_effect_size")
    return pd.concat([source_prop, source_summary, source_effect], ignore_index=True, sort=False)


def draw_auc_ladder(ax: plt.Axes, ladder: pd.DataFrame, *, compact: bool = False) -> None:
    colors = {"celltype": "#777777", "cn_prop": "#4A7088", "cn_coupling": GROUP_COLORS["CR"]}
    y = np.arange(len(ladder))
    ax.axvline(0.5, color="#999999", lw=0.8, ls="--")
    for i, row in ladder.iterrows():
        auc = float(row["auc"])
        low = float(row["auc_bootstrap95_low"])
        high = float(row["auc_bootstrap95_high"])
        tier = row["tier"]
        ax.hlines(i, low, high, color=colors[tier], lw=1.8)
        ax.scatter([auc], [i], color=colors[tier], s=32, zorder=3)
        ax.text(min(high + 0.018, 0.99), i, f"{auc:.3f}", va="center", ha="left", fontsize=6.8)
    ax.set_yticks(y)
    if compact:
        labels = ["Celltype", "CN prop.", "CN coupling"]
    else:
        labels = ladder["display_label"].tolist()
    ax.set_yticklabels(labels)
    ax.tick_params(axis="y", labelsize=6.2 if compact else 6.6, pad=1.0)
    ax.invert_yaxis()
    ax.set_xlim(0.40, 1.00)
    ax.set_xlabel("Nested LOOCV AUC")
    ax.set_title("Higher-order niche coupling needed")
    ax.text(-0.18, 1.08, "f", transform=ax.transAxes, fontsize=13, weight="bold", va="top")


def reconstruct_roc_source() -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    rev = _wp31b_module()
    folds, selected = rev.reconstruct_wp22_folds()
    summary = pd.read_csv(DELIVERABLE_DIR / "WP2.2_nested_prediction_summary.csv")
    row = summary[
        summary["analysis_set"].astype(str).eq("primary_threshold30k")
        & summary["feature_pool"].astype(str).eq("full_candidate_pool")
        & summary["selection_mode"].astype(str).eq("data_driven_topn")
        & summary["top_n"].astype(str).eq("3")
        & summary["model"].astype(str).eq("fold_composite")
    ].iloc[0]
    y = folds["truth"].astype(int).to_numpy()
    s = folds["score"].astype(float).to_numpy()
    fpr, tpr, _ = roc_curve(y, s)
    grid = np.linspace(0, 1, 151)
    interp = np.interp(grid, fpr, tpr)
    interp[0], interp[-1] = 0, 1
    rng = np.random.default_rng(20260708)
    boots = []
    for _ in range(2000):
        idx = rng.choice(len(y), len(y), replace=True)
        if len(np.unique(y[idx])) < 2:
            continue
        bfpr, btpr, _ = roc_curve(y[idx], s[idx])
        curve = np.interp(grid, bfpr, btpr)
        curve[0], curve[-1] = 0, 1
        boots.append(curve)
    arr = np.vstack(boots)
    curve_source = pd.DataFrame({"source_type": "roc_ci_grid", "fpr": grid, "tpr": interp, "tpr_low": np.quantile(arr, 0.025, axis=0), "tpr_high": np.quantile(arr, 0.975, axis=0)})
    fold_source = folds.copy()
    fold_source.insert(0, "source_type", "fold_prediction")
    source = pd.concat([fold_source, curve_source], ignore_index=True, sort=False)
    source["auc_wp22_summary"] = float(row["auc"])
    source["auc_bootstrap95_low_wp22_summary"] = float(row["auc_bootstrap95_low"])
    source["auc_bootstrap95_high_wp22_summary"] = float(row["auc_bootstrap95_high"])
    source["permutation_p_wp22_summary"] = float(row["permutation_p"])
    return source, row, selected


def draw_prediction_roc(ax: plt.Axes, roc_source: pd.DataFrame, row: pd.Series) -> None:
    curve = roc_source[roc_source["source_type"].eq("roc_ci_grid")].copy()
    ax.fill_between(curve["fpr"].astype(float), curve["tpr_low"].astype(float), curve["tpr_high"].astype(float), color="#A6BDD9", alpha=0.34, linewidth=0)
    ax.plot(curve["fpr"].astype(float), curve["tpr"].astype(float), color=GROUP_COLORS["NR"], lw=1.6)
    ax.plot([0, 1], [0, 1], color="#999999", lw=0.8, ls="--")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Potential baseline predictor")
    ax.text(-0.18, 1.08, "g", transform=ax.transAxes, fontsize=13, weight="bold", va="top")
    ax.text(
        0.05,
        0.95,
        f"AUC {float(row['auc']):.3f}\n95% CI {float(row['auc_bootstrap95_low']):.2f}-{float(row['auc_bootstrap95_high']):.2f}\nperm p={float(row['permutation_p']):.3f}\nwarrants validation",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.1,
        bbox=dict(facecolor="white", edgecolor="#DDDDDD", alpha=0.93, pad=1.6),
    )


def build_crop_candidate_table(spatial: pd.DataFrame, metrics: pd.DataFrame, config: dict) -> pd.DataFrame:
    records = []
    metric_lookup = metrics.set_index("sample")
    samples = [("6c", s) for s in config["panel_6c_samples"]] + [("6d", s) for s in config["panel_6d_samples"]]
    for panel, sample in samples:
        sample_df = spatial[spatial["sample"].eq(sample)].copy()
        if sample_df.empty:
            continue
        outcome = str(metric_lookup.loc[sample, "Outcome_short"])
        manual = config.get("crop_overrides", {}).get(sample)
        candidates = smart_crop_candidates(
            sample_df,
            panel=panel,
            outcome=outcome,
            sample=sample,
            crop_size=float(config["crop_size_px"][panel]),
            n_candidates=int(config.get("crop_candidate_n", 3)),
            min_cells=int(config["crop_min_cells"][panel]),
            manual_override=manual,
        )
        for metric in ["cn5_prop", "cd21b_same_type_knn_fraction", "coupling_z_CN1_CN7", "coupling_z_CN3_CN7", "coupling_z_CN4_CN7"]:
            if metric in metric_lookup.columns:
                candidates[metric] = metric_lookup.loc[sample, metric]
        records.append(candidates)
    out = pd.concat(records, ignore_index=True) if records else pd.DataFrame()
    out.to_csv(PANEL_DIR / "WP3.2_smart_crop_candidates.csv", index=False)
    return out


def draw_spatial_grid(fig: plt.Figure, subspec, panel: str, samples: list[str], spatial: pd.DataFrame, metrics: pd.DataFrame, crops: pd.DataFrame, *, letter: str, title: str, style: dict | None = None, link_by_cn: bool | None = None) -> pd.DataFrame:
    st = {**CONFIG["crop_point_style"], **(style or {})}
    base_size = float(st["point_size"])
    link_cn = bool(CONFIG.get("link_emphasis_by_cn", True) if link_by_cn is None else link_by_cn)
    n = len(samples)
    sub = subspec.subgridspec(2, n + 1, width_ratios=[1] * n + [0.78], hspace=0.05, wspace=0.045)
    metric_lookup = metrics.set_index("sample")
    source_parts = []
    for i, sample in enumerate(samples):
        sample_df = spatial[spatial["sample"].eq(sample)].copy()
        outcome = str(metric_lookup.loc[sample, "Outcome_short"])
        crop_candidates = crops[crops["panel"].eq(panel) & crops["sample"].eq(sample)].sort_values("candidate_rank")
        crop = crop_candidates.iloc[0]
        crop_df = crop_dataframe(sample_df, crop)
        if len(crop_df) > 22000:
            crop_df = crop_df.sample(n=22000, random_state=20260708 + i)
        spec = panel_crop_targets(panel, outcome)
        cn_highlights = set(spec["cn_targets"]) | ({"CN4", "CN5", "CN6", "CN3", "CN1"} if panel == "6c" else {"CN1", "CN3", "CN7"})
        cell_highlights = set(spec["cell_targets"])
        ax_ph = fig.add_subplot(sub[0, i])
        ax_cn = fig.add_subplot(sub[1, i])
        ph_kwargs = dict(emphasis_col="CN_k10", emphasis_values=cn_highlights) if link_cn else {}
        scatter_crop(ax_ph, crop_df, color_col="celltype", highlights=cell_highlights, palette=CELLTYPE_COLORS, crop=crop, point_size=base_size, style=st, **ph_kwargs)
        scatter_crop(ax_cn, crop_df, color_col="CN_k10", highlights=cn_highlights, palette=CN_COLORS, crop=crop, point_size=base_size, style=st)
        add_locator(ax_ph, sample_df, crop, GROUP_COLORS[outcome])
        draw_scale_bar(ax_cn, crop)
        ax_ph.set_title(f"{sample} {outcome}", color=GROUP_COLORS[outcome], fontsize=7.2, pad=2)
        if i == 0:
            ax_ph.text(-0.08, 0.5, "phenotype", transform=ax_ph.transAxes, ha="right", va="center", fontsize=6.6)
            ax_cn.text(-0.08, 0.5, "CN", transform=ax_cn.transAxes, ha="right", va="center", fontsize=6.6)
            ax_ph.text(-0.26, 1.08, letter, transform=ax_ph.transAxes, fontsize=13, weight="bold", va="top")
        label = f"crop #{int(crop['candidate_rank'])}, n={int(crop['n_cells'])}"
        ax_cn.text(0.02, 0.05, label, transform=ax_cn.transAxes, ha="left", va="bottom", fontsize=5.6, bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=1.1))
        part = crop_df.copy()
        part.insert(0, "source_type", f"{panel}_crop_cells")
        part["crop_candidate_rank"] = int(crop["candidate_rank"])
        for col in ["xmin", "xmax", "ymin", "ymax", "score", "target_score", "interface_count", "n_cells"]:
            part[f"crop_{col}"] = crop[col]
        source_parts.append(part)
    leg_ph = fig.add_subplot(sub[0, n])
    leg_cn = fig.add_subplot(sub[1, n])
    for ax in [leg_ph, leg_cn]:
        ax.set_axis_off()
    payload = legend_payloads()
    if panel == "6d":
        phenotype_items = ["Vascular", "HER2+ Tumor", "PanCK+ Tumor", "CD8 T cell"]
        cn_items = ["CN7", "CN1", "CN3"]
    else:
        phenotype_items = payload["phenotype"][:4]
        cn_items = ["CN5", "CN4", "CN6", "CN3", "CN1"]
    ph_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=CELLTYPE_COLORS[item], markersize=4.5, label=item) for item in phenotype_items]
    ph_handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor=NEUTRAL_CELL, markersize=4.5, label="other cells"))
    cn_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=CN_COLORS[item], markersize=4.5, label=f"{item} {CN_SHORT[item]}") for item in cn_items]
    cn_handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor=NEUTRAL_CELL, markersize=4.5, label="other CN"))
    leg_ph.legend(handles=ph_handles, loc="center left", frameon=False, title="Phenotype", title_fontsize=6.4, fontsize=5.9)
    leg_cn.legend(handles=cn_handles, loc="center left", frameon=False, title="CN", title_fontsize=6.4, fontsize=5.9)
    leg_ph.text(0.00, 1.06, title, transform=leg_ph.transAxes, ha="left", va="bottom", fontsize=7.2, weight="bold")
    return pd.concat(source_parts, ignore_index=True) if source_parts else pd.DataFrame()


def export_new_panels(prop_rows: pd.DataFrame, tls_summary: pd.DataFrame, stats13: pd.DataFrame, stats14: pd.DataFrame, ladder: pd.DataFrame, roc_source: pd.DataFrame, roc_row: pd.Series, selected_features: pd.DataFrame, config: dict) -> dict[str, dict[str, str]]:
    exports = {}
    fig = plt.figure(figsize=(7.4, 2.65))
    source_6e = draw_tls_presence_panel(fig, fig.add_gridspec(1, 1)[0, 0], prop_rows, tls_summary, stats13, stats14)
    exports["6e"] = save_figure(fig, "fig6e_tls_presence_and_vascular_coupling", source_6e, dpi=int(config["export_dpi"]))

    fig, ax = plt.subplots(figsize=(3.7, 2.25))
    draw_auc_ladder(ax, ladder)
    fig.subplots_adjust(left=0.42, right=0.97, top=0.88, bottom=0.23)
    exports["6f"] = save_figure(fig, "fig6f_auc_ladder", ladder, dpi=int(config["export_dpi"]))

    fig, ax = plt.subplots(figsize=(3.15, 2.95))
    draw_prediction_roc(ax, roc_source, roc_row)
    fig.subplots_adjust(left=0.18, right=0.97, top=0.88, bottom=0.17)
    exports["6g"] = save_figure(fig, "fig6g_prediction_roc_nested", roc_source, dpi=int(config["export_dpi"]))
    selected_features.to_csv(PANEL_DIR / "fig6g_prediction_roc_nested_selected_features_source_data.csv", index=False)
    return exports


def build_main_assembly(config: dict, metrics: pd.DataFrame, spatial: pd.DataFrame, crops: pd.DataFrame, prop_rows: pd.DataFrame, tls_summary: pd.DataFrame, stats13: pd.DataFrame, stats14: pd.DataFrame, ladder: pd.DataFrame, roc_source: pd.DataFrame, roc_row: pd.Series) -> dict[str, str]:
    fig = plt.figure(figsize=tuple(config["figure_size"]), facecolor=config.get("background", "white"))
    gs_cfg = config["gridspec"]
    gs = fig.add_gridspec(
        4,
        6,
        height_ratios=gs_cfg["height_ratios"],
        width_ratios=gs_cfg["width_ratios"],
        hspace=gs_cfg["hspace"],
        wspace=gs_cfg["wspace"],
    )
    m = gs_cfg["margins"]
    fig.subplots_adjust(left=m["left"], right=m["right"], top=m["top"], bottom=m["bottom"])

    ax_a = fig.add_subplot(gs[0, :2])
    source_a = draw_workflow_placeholder(ax_a, config)
    sub_b = gs[0, 2:].subgridspec(1, 2, width_ratios=[1, 0.035], wspace=0.08)
    ax_b = fig.add_subplot(sub_b[0, 0])
    cax_b = fig.add_subplot(sub_b[0, 1])
    source_b = draw_composition_heatmap(ax_b, cax_b)
    source_c = draw_spatial_grid(fig, gs[1, :], "6c", config["panel_6c_samples"], spatial, metrics, crops, letter="c", title="Immune/TLS organization vs tumor-dense CN territories", style=config.get("crop_point_style"), link_by_cn=config.get("link_emphasis_by_cn"))
    source_d = draw_spatial_grid(fig, gs[2, :], "6d", config["panel_6d_samples"], spatial, metrics, crops, letter="d", title="CN-CN vascular niche coupling", style=config.get("crop_point_style"), link_by_cn=config.get("link_emphasis_by_cn"))
    source_e = draw_tls_presence_panel(fig, gs[3, :3], prop_rows, tls_summary, stats13, stats14, compact=True)
    ax_f = fig.add_subplot(gs[3, 3:5])
    draw_auc_ladder(ax_f, ladder, compact=True)
    ax_g = fig.add_subplot(gs[3, 5])
    draw_prediction_roc(ax_g, roc_source, roc_row)

    manifest = pd.concat(
        [
            source_a.assign(panel="6a"),
            source_b.assign(panel="6b"),
            source_c.assign(panel="6c"),
            source_d.assign(panel="6d"),
            source_e.assign(panel="6e"),
            ladder.assign(panel="6f", source_type="auc_ladder"),
            roc_source.assign(panel="6g"),
        ],
        ignore_index=True,
        sort=False,
    )
    basename = "fig6_main_assembly"
    paths = save_figure(fig, basename, manifest, dpi=int(config["export_dpi"]))
    preview = DELIVERABLE_DIR / "fig6_assembly_preview.png"
    shutil.copyfile(paths["png"], preview)
    paths["preview"] = str(preview)
    return paths


def write_notebook(config: dict) -> Path:
    nb_path = DELIVERABLE_DIR / "fig6_assembly.ipynb"
    config_text = json.dumps(config, indent=4)
    cells = [
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "# Fig. 6 CN-centric assembly\n",
                "\n",
                "Edit `CONFIG`, then run all cells to rebuild 6a-6g, preview, and multi-format exports.\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                f"CONFIG = {config_text}\n",
                f"COLORS_15 = {json.dumps(COLORS_15)}\n",
                f"CN_COLORS = {{f'CN{{i}}': COLORS_15[i] for i in range(10)}}\n",
                f"CELLTYPE_COLORS = {json.dumps(CELLTYPE_COLORS, indent=4)}\n",
                f"GROUP_COLORS = {json.dumps(GROUP_COLORS)}\n",
                f"TIME_COLORS = {json.dumps(TIME_COLORS)}\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "from pathlib import Path\n",
                "import sys\n",
                "project_dir = Path(CONFIG['project_dir'])\n",
                "# script_dir 决定用哪一份 WP3_2/WP3_1b(默认 deliverables;可指向项目文件夹副本)\n",
                "script_dir = CONFIG.get('script_dir') or str(project_dir / 'deliverables')\n",
                "sys.path.insert(0, str(script_dir))\n",
                "import assemble_fig6_panels as fig6\n",
                "fig6.run(CONFIG)\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "from IPython.display import Image, display\n",
                "display(Image(filename=str(Path(CONFIG['project_dir']) / 'deliverables' / 'fig6_assembly_preview.png')))\n",
            ],
        },
    ]
    nb = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    nb_path.write_text(json.dumps(nb, indent=2), encoding="utf-8")
    return nb_path


def write_readme(config: dict, ladder: pd.DataFrame, tls_summary: pd.DataFrame, crops: pd.DataFrame, exports: dict[str, str]) -> Path:
    readme = DELIVERABLE_DIR / "WP3.2_assembly_readme.md"
    auc_lines = "\n".join(
        f"- {r.display_label}: AUC {r.auc:.3f} (95% CI {r.auc_bootstrap95_low:.2f}-{r.auc_bootstrap95_high:.2f}); model `{r.model}`"
        for r in ladder.itertuples(index=False)
    )
    crop_lines = crops[crops["candidate_rank"].eq(1)][["panel", "sample", "Outcome_short", "center_x", "center_y", "n_cells", "score"]].to_string(index=False)
    text = f"""# WP3.2 Fig. 6 Assembly

## Outputs
- Notebook: `{DELIVERABLE_DIR / 'fig6_assembly.ipynb'}`
- Preview: `{DELIVERABLE_DIR / 'fig6_assembly_preview.png'}`
- Main assembly: `{PANEL_DIR / 'fig6_main_assembly.svg'}` plus PDF/PNG/TIFF/source CSV.
- New panels: `fig6e_tls_presence_and_vascular_coupling.*`, `fig6f_auc_ladder.*`, `fig6g_prediction_roc_nested.*`.

## Rebuild
```bash
"$HOME/miniconda3/bin/conda" run --no-capture-output -n scimap python {DELIVERABLE_DIR / SCRIPT_NAME}
```

The notebook has top-level `CONFIG` for figure size, gridspec spacing, background, font, dpi, 6c/6d sample choices, crop sizes, crop overrides, and sample-count constants. Palette constants are centralized in the notebook and script.

## Sample Counts
- Experimental: 48 patients / 82 mIHC sections.
- Downstream: 45 patients / 72 sections / 23.1M reconciled cells.
- Baseline Pre: 40 sections (CR22/NR18); prediction threshold >=30k: 36 sections (CR19/NR17).

## 6e TLS Presence
- CN5 >=1% presence: CR {int(tls_summary.loc['CR', 'present'])}/{int(tls_summary.loc['CR', 'n'])}, NR {int(tls_summary.loc['NR', 'present'])}/{int(tls_summary.loc['NR', 'n'])}; Fisher p={float(tls_summary.attrs.get('fisher_p', float('nan'))):.3f}.
- 6e also includes vascular coupling z effect sizes for CN1-CN7, CN3-CN7, and CN4-CN7.

## 6f AUC Ladder
{auc_lines}

## Spatial Assembly
- Phenotype and CN legends are split into separate legend blocks to avoid manuscript-palette color collisions.
- Smart crop candidates are written to `{PANEL_DIR / 'WP3.2_smart_crop_candidates.csv'}`; the first-ranked crop is used unless overridden in `CONFIG['crop_overrides']`.

Default selected crops:
```text
{crop_lines}
```

## Main Export
{json.dumps(exports, indent=2)}
"""
    readme.write_text(text, encoding="utf-8")
    return readme


def audit_outputs(paths: list[Path]) -> pd.DataFrame:
    rows = []
    for path in paths:
        rows.append({"path": str(path), "exists": path.exists(), "bytes": path.stat().st_size if path.exists() else 0})
    audit = pd.DataFrame(rows)
    audit.to_csv(DELIVERABLE_DIR / "WP3.2_export_audit.csv", index=False)
    return audit


def run(config: dict | None = None) -> dict[str, object]:
    config = {**CONFIG, **(config or {})}
    configure_matplotlib(config)
    ensure_dirs()

    metrics = load_primary_metrics()
    stats13 = pd.read_csv(DELIVERABLE_DIR / "WP1.3_CRvsNR_stats.csv")
    stats14 = pd.read_csv(DELIVERABLE_DIR / "WP1.4_CRvsNR_stats.csv")
    prop = pd.read_csv(DELIVERABLE_DIR / "WP1.3_sample_proportion_long.csv")
    zero = pd.read_csv(DELIVERABLE_DIR / "WP1.3_cd21_zero_inflation.csv")
    summary22 = pd.read_csv(DELIVERABLE_DIR / "WP2.2_nested_prediction_summary.csv")

    ladder = extract_auc_ladder(summary22)
    prop_rows, tls_summary = compute_tls_presence(prop, zero, threshold=0.01)
    roc_source, roc_row, selected_features = reconstruct_roc_source()
    samples = sorted(set(config["panel_6c_samples"] + config["panel_6d_samples"]))
    spatial = load_spatial_for_samples(samples)
    crops = build_crop_candidate_table(spatial, metrics, config)

    panel_exports = export_new_panels(prop_rows, tls_summary, stats13, stats14, ladder, roc_source, roc_row, selected_features, config)
    assembly_paths = build_main_assembly(config, metrics, spatial, crops, prop_rows, tls_summary, stats13, stats14, ladder, roc_source, roc_row)
    notebook = write_notebook(config)
    readme = write_readme(config, ladder, tls_summary, crops, assembly_paths)
    local_script = Path(__file__)
    if local_script.resolve() != (DELIVERABLE_DIR / SCRIPT_NAME).resolve():
        shutil.copyfile(local_script, DELIVERABLE_DIR / SCRIPT_NAME)
    audit_paths = [Path(p) for exports in panel_exports.values() for p in exports.values()] + [Path(p) for p in assembly_paths.values()] + [notebook, readme, DELIVERABLE_DIR / SCRIPT_NAME, PANEL_DIR / "WP3.2_smart_crop_candidates.csv"]
    audit = audit_outputs(audit_paths)
    return {
        "panel_exports": panel_exports,
        "assembly_paths": assembly_paths,
        "notebook": str(notebook),
        "readme": str(readme),
        "auc_ladder": ladder.to_dict(orient="records"),
        "tls_presence": {
            "CR": f"{int(tls_summary.loc['CR', 'present'])}/{int(tls_summary.loc['CR', 'n'])}",
            "NR": f"{int(tls_summary.loc['NR', 'present'])}/{int(tls_summary.loc['NR', 'n'])}",
            "fisher_p": float(tls_summary.attrs.get("fisher_p", np.nan)),
        },
        "audit_failed": audit[(~audit["exists"]) | (audit["bytes"].le(0))].to_dict(orient="records"),
    }


def main() -> None:
    result = run(CONFIG)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
