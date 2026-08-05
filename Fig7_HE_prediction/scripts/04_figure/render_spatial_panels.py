#!/usr/bin/env python3
"""Render corrected Figure 7e-h panels with raster maps and vector text."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from PIL import Image
import numpy as np
import pandas as pd


PALETTE = {
    "Vascular": "#009E73",
    "PanCK signal": "#E69F00",
    "HER2": "#CC79A7",
    "NR": "#0072B2",
    "CR": "#D55E00",
    "grey": "#9E9E9E",
    "dark": "#2F3742",
}
ORDER = ["Vascular", "PanCK signal", "HER2"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_figure(fig: plt.Figure, stem: Path, *, dpi: int = 600) -> dict[str, str]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = {
        "pdf": stem.with_suffix(".pdf"),
        "svg": stem.with_suffix(".svg"),
        "png": stem.with_suffix(".png"),
        "tiff": stem.with_suffix(".tiff"),
    }
    fig.savefig(paths["pdf"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["svg"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["png"], dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(paths["tiff"], dpi=dpi, bbox_inches="tight", facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    return {key: sha256_file(path) for key, path in paths.items()}


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.7, length=2.5, pad=2)


def draw_panel_e(ax: plt.Axes, source: pd.DataFrame, *, panel_label: bool = True) -> None:
    data = source.set_index("phenotype").loc[ORDER].reset_index()
    x = np.arange(len(data))
    for index, row in data.iterrows():
        color = PALETTE[str(row.phenotype)]
        ax.plot(
            [index, index], [float(row.bootstrap_95ci_lower), float(row.bootstrap_95ci_upper)],
            color="#B7B7B7", lw=1.0, zorder=1,
        )
        ax.scatter(index, float(row.max_shuffled_median_spearman), marker="x", s=28, color="#B6B6B6", lw=1.1, zorder=2)
        ax.scatter(index, float(row.median_patient_spearman), s=62, color=color, edgecolor="#4A4A4A", lw=0.65, zorder=3)
        ax.text(index, float(row.median_patient_spearman) + 0.045, f"{float(row.median_patient_spearman):.3f}",
                ha="center", va="bottom", fontsize=6.8, color="#555555")
    ax.axhline(0.30, color="#B8B8B8", lw=0.75, ls=(0, (3, 2)), zorder=0)
    ax.text(2.38, 0.307, "ρ=0.30 reference", ha="right", va="bottom", color="#A3A3A3", fontsize=6.2)
    ax.set_xticks(x, ["Vascular", "PanCK\nsignal", "HER2"])
    ax.set_ylim(0, max(0.66, float(data.bootstrap_95ci_upper.max()) + 0.08))
    ax.set_ylabel("Patient-held-out median Spearman")
    ax.set_title("H&E captures regional phenotypes", fontsize=8.4, pad=4)
    clean_axis(ax)
    if panel_label:
        ax.text(-0.16, 1.08, "e", transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")


def draw_panel_f(ax: plt.Axes, metrics: pd.DataFrame, summary: pd.DataFrame, *, panel_label: bool = True) -> None:
    name_map = {"Epithelial (PanCK + HER2)": "PanCK signal"}
    metrics = metrics.copy()
    summary = summary.copy()
    metrics["display"] = metrics.phenotype.replace(name_map)
    summary["display"] = summary.phenotype.replace(name_map)
    patient_order = {patient: index for index, patient in enumerate(sorted(metrics.patient_id.astype(str).unique()))}
    for index, phenotype in enumerate(ORDER):
        group = metrics.loc[metrics.display.eq(phenotype)].copy()
        jitter = np.asarray([((patient_order[str(patient)] * 37) % 17 - 8) / 70 for patient in group.patient_id], dtype=float)
        ax.scatter(
            np.full(len(group), index, dtype=float) + jitter,
            group.spearman.to_numpy(float), s=17, color=PALETTE[phenotype], alpha=0.68,
            edgecolor="white", linewidth=0.3, zorder=2,
        )
        row = summary.loc[summary.display.eq(phenotype)].iloc[0]
        median = float(row.median_patient_spearman)
        ax.plot([index, index], [float(row.bootstrap_95ci_lower), float(row.bootstrap_95ci_upper)], color="#4A4A4A", lw=1.2, zorder=3)
        ax.scatter(index, median, s=72, color=PALETTE[phenotype], edgecolor="#414141", lw=0.8, zorder=4)
        positive = int((group.spearman > 0).sum())
        ax.text(index, -0.075, f"{positive}/{len(group)} +", ha="center", va="bottom", fontsize=6.6, color="#8B8B8B")
        ax.text(index, median + 0.055, f"{median:.3f}", ha="center", va="bottom", fontsize=6.8, color="#555555")
    ax.axhline(0, color="#B8B8B8", lw=0.75)
    ax.set_xticks(np.arange(3), ["Vascular", "PanCK\nsignal", "HER2"])
    ax.set_ylim(-0.11, max(0.9, float(metrics.spearman.max()) + 0.10))
    ax.set_ylabel("Within-patient spatial Spearman")
    ax.set_title("M0 contribution co-localizes with regional maps", fontsize=8.4, pad=4)
    clean_axis(ax)
    if panel_label:
        ax.text(-0.16, 1.08, "f", transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")


def raster_tile_map(
    frame: pd.DataFrame,
    values: np.ndarray,
    width: int,
    height: int,
    cmap: mpl.colors.Colormap,
    norm: Normalize,
    *,
    max_dim: int = 1800,
) -> np.ndarray:
    scale = min(float(max_dim) / max(width, height), 1.0)
    out_w = max(1, int(round(width * scale)))
    out_h = max(1, int(round(height * scale)))
    canvas = np.ones((out_h, out_w, 4), dtype=np.float32)
    colors = cmap(norm(np.asarray(values, dtype=float)))
    coords = frame[["x_min", "y_min", "x_max", "y_max"]].to_numpy(float)
    for coordinate, color in zip(coords, colors):
        x0 = int(np.floor(coordinate[0] * out_w / width))
        x1 = max(x0 + 1, int(np.ceil(coordinate[2] * out_w / width)))
        y0 = int(np.floor(coordinate[1] * out_h / height))
        y1 = max(y0 + 1, int(np.ceil(coordinate[3] * out_h / height)))
        canvas[max(0, y0):min(out_h, y1), max(0, x0):min(out_w, x1)] = color
    return canvas


def crop_bounds(frame: pd.DataFrame, width: int, height: int, margin_fraction: float = 0.018) -> tuple[float, float, float, float]:
    coords = frame[["x_min", "y_min", "x_max", "y_max"]].to_numpy(float)
    x0, y0 = coords[:, 0].min(), coords[:, 1].min()
    x1, y1 = coords[:, 2].max(), coords[:, 3].max()
    margin = max(width, height) * margin_fraction
    return max(0, x0 - margin), min(width, x1 + margin), max(0, y0 - margin), min(height, y1 + margin)


def scale_bar(ax: plt.Axes, bounds: tuple[float, float, float, float], mpp: float = 0.503) -> None:
    x0, x1, y0, y1 = bounds
    length = 2000.0 / mpp
    start_x = x0 + 0.05 * (x1 - x0)
    y = y1 - 0.06 * (y1 - y0)
    ax.plot([start_x, start_x + length], [y, y], color="black", lw=1.5, solid_capstyle="butt")
    ax.text(start_x + length / 2, y - 0.025 * (y1 - y0), "2 mm", ha="center", va="bottom", fontsize=5.8)


def draw_map_panel(
    fig: plt.Figure,
    slot,
    frame: pd.DataFrame,
    concordance: pd.DataFrame,
    overview_ledger: pd.DataFrame,
    image_root: Path,
    *,
    phenotype: str,
    panel: str,
    cmap_name: str,
    vmin: float,
    vmax: float,
) -> None:
    grid = slot.subgridspec(3, 3, height_ratios=[0.12, 1, 1], hspace=0.12, wspace=0.08)
    headers = ["H&E", f"{phenotype}\nmeasured", f"{phenotype}\nprediction"]
    for column, title in enumerate(headers):
        ax = fig.add_subplot(grid[0, column])
        ax.axis("off")
        ax.text(0.5, 0.2, title, ha="center", va="center", fontsize=7.4, fontweight="bold")
    cmap = mpl.colormaps[cmap_name]
    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    displayed = frame[["sample_id", "response", "response_label"]].drop_duplicates().sort_values("response")
    for row_index, identity in enumerate(displayed.itertuples(index=False), start=1):
        sample = str(identity.sample_id)
        sample_frame = frame.loc[frame.sample_id.astype(str).eq(sample)].copy()
        ledger = overview_ledger.loc[overview_ledger.sample_id.astype(str).eq(sample)].iloc[0]
        width, height = int(ledger.wsi_width), int(ledger.wsi_height)
        image = np.asarray(Image.open(image_root / f"{sample}_HE_overview.png").convert("RGB"))
        measured_map = raster_tile_map(sample_frame, sample_frame.measured.to_numpy(float), width, height, cmap, norm)
        predicted_map = raster_tile_map(sample_frame, sample_frame.predicted.to_numpy(float), width, height, cmap, norm)
        bounds = crop_bounds(sample_frame, width, height)
        maps = [image, measured_map, predicted_map]
        for column, raster in enumerate(maps):
            ax = fig.add_subplot(grid[row_index, column])
            ax.imshow(raster, extent=[0, width, height, 0], interpolation="nearest", rasterized=True)
            ax.set_xlim(bounds[0], bounds[1])
            ax.set_ylim(bounds[3], bounds[2])
            ax.set_aspect("equal")
            ax.axis("off")
            if column == 0:
                color = PALETTE[str(identity.response_label)]
                ax.text(
                    0.018, 0.975, f"{identity.response_label} · {sample.split('-')[0]}",
                    transform=ax.transAxes, ha="left", va="top", fontsize=6.7, color="white", fontweight="bold",
                    bbox={"boxstyle": "square,pad=0.16", "facecolor": color, "edgecolor": "none", "alpha": 0.92},
                )
                scale_bar(ax, bounds)
                for spine in ax.spines.values():
                    spine.set_visible(True)
                    spine.set_color(color)
                    spine.set_linewidth(1.0)
            if column == 2:
                row = concordance.loc[
                    concordance.panel.eq(panel) & concordance.sample_id.astype(str).eq(sample)
                ].iloc[0]
                ax.text(
                    0.02, 0.98, f"ρ={float(row.spearman):.2f}\nn={int(row.n_tiles):,}",
                    transform=ax.transAxes, ha="left", va="top", fontsize=6.2, color="#30343B",
                    bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
                )
    cax = fig.add_axes([slot.get_position(fig).x0 + 0.055 * slot.get_position(fig).width,
                        slot.get_position(fig).y0 - 0.030,
                        slot.get_position(fig).width * 0.89, 0.012])
    colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    colorbar.set_label("Phenotype signal [log1p(cells per mm²)]", fontsize=6.3, labelpad=1)
    colorbar.ax.tick_params(labelsize=5.7, width=0.5, length=2, pad=1)
    fig.text(slot.get_position(fig).x0 - 0.015, slot.get_position(fig).y1 + 0.01, panel,
             fontsize=10, fontweight="bold", ha="left", va="top")


def render_individual_map(
    output: Path,
    name: str,
    frame: pd.DataFrame,
    concordance: pd.DataFrame,
    overview_ledger: pd.DataFrame,
    image_root: Path,
    phenotype: str,
    panel: str,
    cmap: str,
    vmin: float,
    vmax: float,
) -> dict[str, str]:
    fig = plt.figure(figsize=(7.2, 3.4))
    grid = fig.add_gridspec(1, 1, left=0.045, right=0.985, top=0.96, bottom=0.10)
    draw_map_panel(fig, grid[0, 0], frame, concordance, overview_ledger, image_root,
                   phenotype=phenotype, panel=panel, cmap_name=cmap, vmin=vmin, vmax=vmax)
    hashes = save_figure(fig, output / name)
    plt.close(fig)
    return hashes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--images", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.linewidth": 0.75, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "savefig.transparent": False,
    })
    status = json.loads((args.source / "complete.json").read_text(encoding="utf-8"))
    if status.get("status") != "PASS":
        raise ValueError("Corrected source-data completion status is not PASS")
    e = pd.read_csv(args.source / "panel_e_regional.tsv", sep="\t")
    f_metrics = pd.read_csv(args.source / "panel_f_patient_metrics.tsv", sep="\t")
    f_summary = pd.read_csv(args.source / "panel_f_summary.tsv", sep="\t")
    g = pd.read_csv(args.source / "panel_g_external_vascular.tsv.gz", sep="\t", low_memory=False)
    h = pd.read_csv(args.source / "panel_h_internal_her2.tsv.gz", sep="\t", low_memory=False)
    concordance = pd.read_csv(args.source / "panel_gh_concordance.tsv", sep="\t")
    scales = pd.read_csv(args.source / "panel_gh_scales.tsv", sep="\t").set_index("panel")
    ledger = pd.read_csv(args.images / "overview_ledger.tsv", sep="\t")
    if set(e.phenotype) != set(ORDER) or set(e.n_patients.astype(int)) != {30}:
        raise ValueError("Panel e source contract failed")
    if set(f_summary.n_patients.astype(int)) != {13} or not np.allclose(f_summary.positive_patient_fraction, 1.0):
        raise ValueError("Panel f source contract failed")
    for panel_name, frame, patients in (("g", g, {"pt71", "pt80"}), ("h", h, {"pt42", "pt65"})):
        if set(frame.patient_id.astype(str)) != patients:
            raise ValueError(f"Panel {panel_name} patient contract failed")
        if frame.duplicated(["patient_id", "tile_id"]).any():
            raise ValueError(f"Panel {panel_name} duplicate tile")
        if not np.isfinite(frame[["measured", "predicted", "x_min", "y_min", "x_max", "y_max"]].to_numpy(float)).all():
            raise ValueError(f"Panel {panel_name} non-finite value")

    outputs = {}
    fig, ax = plt.subplots(figsize=(3.5, 2.25))
    draw_panel_e(ax, e)
    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.19, top=0.84)
    outputs["panel_e"] = save_figure(fig, args.output / "Figure7e_regional")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(3.5, 2.25))
    draw_panel_f(ax, f_metrics, f_summary)
    fig.subplots_adjust(left=0.17, right=0.98, bottom=0.19, top=0.84)
    outputs["panel_f"] = save_figure(fig, args.output / "Figure7f_colocalization")
    plt.close(fig)
    outputs["panel_g"] = render_individual_map(
        args.output, "Figure7g_vascular_maps", g, concordance, ledger, args.images,
        "Vascular", "g", "cividis", float(scales.loc["g", "vmin"]), float(scales.loc["g", "vmax"]),
    )
    outputs["panel_h"] = render_individual_map(
        args.output, "Figure7h_her2_maps", h, concordance, ledger, args.images,
        "HER2", "h", "magma", float(scales.loc["h", "vmin"]), float(scales.loc["h", "vmax"]),
    )

    # Paste-ready e-h block matching the current two-column manuscript layout.
    fig = plt.figure(figsize=(7.2, 5.55))
    outer = fig.add_gridspec(2, 2, height_ratios=[0.78, 1.32], hspace=0.25, wspace=0.24,
                             left=0.055, right=0.985, top=0.97, bottom=0.07)
    draw_panel_e(fig.add_subplot(outer[0, 0]), e)
    draw_panel_f(fig.add_subplot(outer[0, 1]), f_metrics, f_summary)
    draw_map_panel(fig, outer[1, 0], g, concordance, ledger, args.images,
                   phenotype="Vascular", panel="g", cmap_name="cividis",
                   vmin=float(scales.loc["g", "vmin"]), vmax=float(scales.loc["g", "vmax"]))
    draw_map_panel(fig, outer[1, 1], h, concordance, ledger, args.images,
                   phenotype="HER2", panel="h", cmap_name="magma",
                   vmin=float(scales.loc["h", "vmin"]), vmax=float(scales.loc["h", "vmax"]))
    outputs["block_eh"] = save_figure(fig, args.output / "Figure7_eh")
    plt.close(fig)

    # Compact block for direct replacement of the lower portion of the current
    # full Figure 7 page; scientific coordinates remain aspect-locked.
    fig = plt.figure(figsize=(7.2, 3.95))
    outer = fig.add_gridspec(2, 2, height_ratios=[0.62, 1.0], hspace=0.18, wspace=0.24,
                             left=0.055, right=0.985, top=0.975, bottom=0.085)
    draw_panel_e(fig.add_subplot(outer[0, 0]), e)
    draw_panel_f(fig.add_subplot(outer[0, 1]), f_metrics, f_summary)
    draw_map_panel(fig, outer[1, 0], g, concordance, ledger, args.images,
                   phenotype="Vascular", panel="g", cmap_name="cividis",
                   vmin=float(scales.loc["g", "vmin"]), vmax=float(scales.loc["g", "vmax"]))
    draw_map_panel(fig, outer[1, 1], h, concordance, ledger, args.images,
                   phenotype="HER2", panel="h", cmap_name="magma",
                   vmin=float(scales.loc["h", "vmin"]), vmax=float(scales.loc["h", "vmax"]))
    outputs["block_eh_compact"] = save_figure(fig, args.output / "Figure7_eh_compact")
    plt.close(fig)

    # Page-fit block for replacing exactly the original e-h footprint in the
    # one-page Figure 7 PDF (the page itself retains whitespace below the art).
    fig = plt.figure(figsize=(7.2, 4.5))
    outer = fig.add_gridspec(2, 2, height_ratios=[0.68, 1.08], hspace=0.21, wspace=0.24,
                             left=0.055, right=0.985, top=0.975, bottom=0.080)
    draw_panel_e(fig.add_subplot(outer[0, 0]), e)
    draw_panel_f(fig.add_subplot(outer[0, 1]), f_metrics, f_summary)
    draw_map_panel(fig, outer[1, 0], g, concordance, ledger, args.images,
                   phenotype="Vascular", panel="g", cmap_name="cividis",
                   vmin=float(scales.loc["g", "vmin"]), vmax=float(scales.loc["g", "vmax"]))
    draw_map_panel(fig, outer[1, 1], h, concordance, ledger, args.images,
                   phenotype="HER2", panel="h", cmap_name="magma",
                   vmin=float(scales.loc["h", "vmin"]), vmax=float(scales.loc["h", "vmax"]))
    outputs["block_eh_pagefit"] = save_figure(fig, args.output / "Figure7_eh_pagefit")
    plt.close(fig)

    qa = {
        "status": "PASS", "backend": "Python/matplotlib", "font": "Arial",
        "map_points_rasterized": True, "text_vector_in_pdf_svg": True,
        "panel_e_n": 30, "panel_f_n": 13,
        "panel_g_samples": ["pt71", "pt80"], "panel_h_samples": ["pt65", "pt42"],
        "outputs": outputs,
    }
    (args.output / "render_QA.json").write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qa, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
