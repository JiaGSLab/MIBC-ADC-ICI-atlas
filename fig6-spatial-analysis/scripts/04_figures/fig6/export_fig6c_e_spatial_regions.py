#!/usr/bin/env python3
"""看全景 -> 读坐标 -> 定窗口 -> 预览裁剪效果。

对应 scimap_IT.ipynb 里的 notebook_preview_full,但接的是 B0 数据
(neighborhood_composition.parquet 的 CN_k10 / celltype / spatial_x / spatial_y),
并且裁剪预览直接调 WP3_1b.scatter_crop,所见即最终图里的样子。

命令行:
  # 1) 全景 + 自动候选窗口
  python preview_sample_regions.py pt53-pre
  # 2) 只把某几个 CN 上色,其余压灰(挑 6c 的 TLS / 6d 的血管时好用)
  python preview_sample_regions.py pt53-pre --highlight CN1 CN7
  # 3) 预览一个窗口(px),两行 phenotype / CN,和最终图一致
  python preview_sample_regions.py pt53-pre --roi 12000 15000 8000 11000

Notebook 里:
  from preview_sample_regions import overview, preview_crop, load_sample
  overview("pt53-pre", highlight=["CN1", "CN7"])
  preview_crop("pt53-pre", [12000, 15000, 8000, 11000])
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import matplotlib

if __name__ == "__main__":
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
CACHE = HERE / ".sample_cache"
OUT = HERE / "region_previews"
CAND_DIR = HERE / "panel_crop_candidates"
# 像素尺寸的唯一来源:与 WP3_1b.PIXEL_SIZE_UM 一致(写进 metadata 的也是这个值)
PIXEL_SIZE_UM = 0.4995            # µm / px
PIXELS_PER_UM = 1.0 / PIXEL_SIZE_UM

# CN 功能命名,画 legend 时带上,方便认
CN_SHORT = {
    "CN0": "Mac-rich", "CN1": "HER2 tumor-Mac", "CN2": "CD4/Treg", "CN3": "PanCK tumor",
    "CN4": "CD8 cytotoxic", "CN5": "CD21+B/TLS", "CN6": "HER2 dense", "CN7": "Vascular",
    "CN8": "B tumor-immune", "CN9": "Fibro-stromal",
}

# 图面显示名 —— 只影响 legend/标注,不动数据列。
# 左边是 B0 里的 celltype 值,右边是手稿图沿用的叫法。
DISPLAY_NAMES = {
    "B cell": "B cell",
    "CD163+ Macrophage": "CD163+ Macrophage",
    "CD163- Macrophage": "CD163- Macrophage",
    "CD21+ B cell": "CD21 B cell",
    "CD4 T cell": "CD4 T cell",
    "CD8 T cell": "CD8 T cell",
    "CD8 exhausted T cell": "CD8 Tex cell",
    "DC": "DC",
    "Vascular": "Endothelial",
    "Fibroblast / SMA+ stromal": "Fibroblast",
    "HER2+ Tumor": "HER2 Epithelial",
    "Treg": "Treg",
    "PanCK+ Tumor": "PanCK Epithelial",
}


# 字体:Arial 优先。服务器没装 Arial,渲染回落到 Liberation Sans
# (与 Arial 度量完全一致的替代字体,字宽相同,换字不会跑版)。
# SVG 里 svg.fonttype='none' 会把整个 fallback 列表写进 font-family,
# Arial 排第一 -> 在你本机 AI 里直接用真 Arial 渲染。
SANS_STACK = ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]


def apply_font(size: float = 7.0) -> None:
    """统一字体设置,文字保持可编辑矢量。"""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": SANS_STACK,
        "font.size": size,
        "pdf.fonttype": 42,     # Type42 = 真文字,不是轮廓
        "ps.fonttype": 42,
        "svg.fonttype": "none", # SVG 里保留 <text> + font-family
    })


def display_name(label: str, *, with_niche: bool = False) -> str:
    """图面显示名。CN 可选带上功能 niche 名。"""
    if label in DISPLAY_NAMES:
        return DISPLAY_NAMES[label]
    if with_niche and label in CN_SHORT:
        return f"{label} {CN_SHORT[label]}"
    return label


def _wp31b():
    spec = importlib.util.spec_from_file_location("wp31b", HERE / "build_fig6b_c_e_panel_sources.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["wp31b"] = module
    spec.loader.exec_module(module)
    return module


def load_sample(sample: str, *, refresh: bool = False) -> pd.DataFrame:
    """取单样本的 spatial 表,带本地缓存(首次约 8 s,之后瞬时)。"""
    CACHE.mkdir(exist_ok=True)
    cached = CACHE / f"{sample}.parquet"
    if cached.exists() and not refresh:
        return pd.read_parquet(cached)
    df = _wp31b().load_spatial_for_samples([sample])
    df.to_parquet(cached, index=False)
    return df


def _crop_subset(sample: str, roi: list[float]) -> tuple[pd.DataFrame, dict[str, float]]:
    """取窗口内的细胞 + crop dict。notebook 里做样式微调时直接用。"""
    xmin, xmax, ymin, ymax = [float(v) for v in roi]
    crop = {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax}
    df = load_sample(sample)
    sub = df[(df.spatial_x >= xmin) & (df.spatial_x <= xmax)
             & (df.spatial_y >= ymin) & (df.spatial_y <= ymax)]
    if sub.empty:
        raise ValueError(f"{sample} 在窗口 {roi} 内没有细胞,检查坐标")
    return sub, crop


def candidate_boxes(sample: str) -> pd.DataFrame:
    """已有的自动候选窗口(WP3.2 smart crop + 6c/6d crop_metadata)。"""
    parts = []
    for name, tag in [
        ("WP3.2_smart_crop_candidates.csv", "smart"),
        ("fig6c_spatial_immune_tls_crop_candidates_crop_metadata.csv", "6c_meta"),
        ("fig6d_spatial_vascular_axis_crop_candidates_crop_metadata.csv", "6d_meta"),
    ]:
        path = CAND_DIR / name
        if not path.exists():
            continue
        df = pd.read_csv(path)
        df = df[df["sample"].eq(sample)]
        if df.empty:
            continue
        rank_col = "candidate_rank" if "candidate_rank" in df.columns else "selection_rank"
        keep = df[["panel", rank_col, "xmin", "xmax", "ymin", "ymax"]].copy()
        keep.columns = ["panel", "rank", "xmin", "xmax", "ymin", "ymax"]
        keep["source"] = tag
        parts.append(keep)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        columns=["panel", "rank", "xmin", "xmax", "ymin", "ymax", "source"])


# --------------------------------------------------------------------------- #
def overview(
    sample: str,
    *,
    by: str = "CN_k10",
    highlight: list[str] | None = None,
    boxes: bool = True,
    max_points: int = 400_000,
    figsize: tuple[float, float] = (11.0, 9.0),
    save: bool = True,
):
    """全景散点 + 坐标网格 + 已有候选窗口。

    highlight 非空时,只有这些类别上色,其余压成浅灰 —— 挑区域时最好用。
    """
    apply_font()
    m = _wp31b()
    palette = m.CN_COLORS if by == "CN_k10" else m.CELLTYPE_COLORS
    df = load_sample(sample)
    if len(df) > max_points:
        df = df.sample(n=max_points, random_state=17)

    fig, ax = plt.subplots(figsize=figsize)
    if highlight:
        hl = set(highlight)
        bg = df[~df[by].isin(hl)]
        ax.scatter(bg["spatial_x"], bg["spatial_y"], s=1.0, color="#E2E2E2",
                   alpha=0.45, linewidths=0, rasterized=True, zorder=2)
        for label in sorted(hl):
            part = df[df[by].eq(label)]
            if part.empty:
                continue
            ax.scatter(part["spatial_x"], part["spatial_y"], s=4.0,
                       color=palette.get(label, "#888888"), alpha=0.95,
                       linewidths=0, rasterized=True, zorder=3, label=label)
        shown = sorted(hl)
    else:
        for label in sorted(df[by].unique()):
            part = df[df[by].eq(label)]
            ax.scatter(part["spatial_x"], part["spatial_y"], s=1.6,
                       color=palette.get(label, "#888888"), alpha=0.75,
                       linewidths=0, rasterized=True, label=label)
        shown = sorted(df[by].unique())

    if boxes:
        cand = candidate_boxes(sample)
        for r in cand.itertuples(index=False):
            ax.add_patch(Rectangle((r.xmin, r.ymin), r.xmax - r.xmin, r.ymax - r.ymin,
                                   fill=False, edgecolor="#111111", linewidth=1.0,
                                   linestyle="--", zorder=6))
            ax.text(r.xmin, r.ymin - 60, f"{r.panel}#{int(r.rank)} [{r.source}]",
                    fontsize=7, color="#111111", va="bottom", zorder=7,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))

    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.grid(True, linestyle="--", alpha=0.45, linewidth=0.6)
    ax.set_xlabel("spatial_x (px)   ·   1 px = %.4f µm" % PIXEL_SIZE_UM)
    ax.set_ylabel("spatial_y (px)")
    ax.set_title(f"{sample} — overview ({by})   n={len(df):,}", fontsize=11)
    handles = [Line2D([], [], marker="o", ls="", ms=5, color=palette.get(l, "#888888"),
                      label=display_name(l, with_niche=True)) for l in shown]
    ax.legend(handles=handles, fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5),
              frameon=False, title=by)
    fig.tight_layout()
    if save:
        OUT.mkdir(exist_ok=True)
        tag = "_".join(highlight) if highlight else "all"
        path = OUT / f"{sample}_overview_{by}_{tag}.png"
        fig.savefig(path, dpi=150)
        print("wrote", path)
    return fig, ax


def preview_crop(
    sample: str,
    roi: list[float],
    *,
    focus_cn: list[str],
    cell_highlights: list[str] | None = None,
    cn_color_highlights: list[str] | None = None,
    point_size: float = 1.65,
    other_keep_color: bool = True,
    style: dict | None = None,
    save: bool = True,
):
    """预览一个窗口:上 phenotype / 下 CN,两行都按 focus_cn 决定放大。

    focus_cn                本面板要突出的 CN(唯一的 focus 依据)。
    cell_highlights         默认 None = 所有 celltype 都上自己的色(参考图口径)。
    cn_color_highlights     默认 None = 所有 CN 都上自己的色。
    other_keep_color        背景点保留类别原色(默认 True),只是缩小压淡。
    """
    apply_font()
    m = _wp31b()
    sub, crop = _crop_subset(sample, roi)
    xmin, xmax, ymin, ymax = crop["xmin"], crop["xmax"], crop["ymin"], crop["ymax"]
    cn_hl = set(focus_cn)
    side_um = (xmax - xmin) * PIXEL_SIZE_UM

    common = dict(crop=crop, point_size=point_size, style=style,
                  other_keep_color=other_keep_color,
                  emphasis_col="CN_k10", emphasis_values=cn_hl)
    fig, axes = plt.subplots(2, 1, figsize=(4.6, 9.4))
    m.scatter_crop(axes[0], sub, color_col="celltype",
                   highlights=cell_highlights, palette=m.CELLTYPE_COLORS, **common)
    m.scatter_crop(axes[1], sub, color_col="CN_k10",
                   highlights=cn_color_highlights, palette=m.CN_COLORS, **common)
    axes[0].set_title(f"{sample}  n={len(sub):,}  {side_um:.0f}×{(ymax-ymin)*PIXEL_SIZE_UM:.0f} µm",
                      fontsize=9)
    axes[0].set_ylabel("phenotype", fontsize=8)
    axes[1].set_ylabel("CN", fontsize=8)

    counts = sub["CN_k10"].value_counts()
    txt = "  ".join(f"{c}:{counts.get(c, 0)}" for c in sorted(cn_hl))
    axes[1].text(0.01, -0.04, txt, transform=axes[1].transAxes, fontsize=7,
                 va="top", color="#444444")
    fig.tight_layout()
    if save:
        OUT.mkdir(exist_ok=True)
        path = OUT / f"{sample}_crop_{int(xmin)}_{int(ymin)}.png"
        fig.savefig(path, dpi=170, bbox_inches="tight")
        print("wrote", path)
    print(f"窗口内细胞 {len(sub):,};各重点 CN 计数: "
          + ", ".join(f"{c}={counts.get(c, 0)}" for c in sorted(cn_hl)))
    return fig, axes


def export_crop_panels(
    sample: str,
    roi: list[float],
    *,
    focus_cn: list[str],
    out_stem: str | None = None,
    cell_highlights: list[str] | None = None,
    cn_color_highlights: list[str] | None = None,
    panel_size_in: float = 1.6,
    point_size: float = 1.65,
    dpi: int = 600,
    other_keep_color: bool = True,
    style: dict | None = None,
    out_dir: Path | None = None,
    legend: bool = True,
    legend_scope: str = "all",
    cn_niche_names: bool = False,
    scalebar_um: float | None = 100.0,
) -> dict[str, Path]:
    """逐面板导出矢量图,给 Illustrator 组版用 —— 不出拼版、不出 TIFF。

    每个面板单独一个 PDF + SVG:画布纯净(无标题/坐标/边框、透明底),
    散点用 rasterized=True 内嵌成 600 dpi 位图,文字与其余元素仍是矢量。
    这样单文件通常只有几百 KB,AI 里放十几个也不会炸内存。

    两行的"重点细胞"都由 focus_cn 决定 —— phenotype 行放大的就是落在这些 CN 里的细胞,
    颜色仍按 celltype 上,所以读出来是"这个 niche 由什么细胞组成"。

    另出一份矢量 legend PDF。legend_scope="all"(默认)列出完整类别表,
    这样多个面板共用一份 legend、同一 CN 颜色跨面板一致;"roi" 则只列本窗口出现的。
    """
    apply_font()
    m = _wp31b()
    out_dir = Path(out_dir) if out_dir else (OUT / "ai_export")
    out_dir.mkdir(parents=True, exist_ok=True)
    sub, crop = _crop_subset(sample, roi)
    xmin, xmax, ymin, ymax = crop["xmin"], crop["xmax"], crop["ymin"], crop["ymax"]

    cn_hl = set(focus_cn)
    stem = out_stem or f"{sample}_{int(xmin)}_{int(ymin)}"

    written: dict[str, Path] = {}
    for kind, col, hl, palette in [
        ("phenotype", "celltype", cell_highlights, m.CELLTYPE_COLORS),
        ("CN", "CN_k10", cn_color_highlights, m.CN_COLORS),
    ]:
        fig = plt.figure(figsize=(panel_size_in, panel_size_in))
        ax = fig.add_axes([0, 0, 1, 1])
        m.scatter_crop(ax, sub, color_col=col, highlights=hl, palette=palette,
                       crop=crop, point_size=point_size, style=style,
                       other_keep_color=other_keep_color,
                       emphasis_col="CN_k10", emphasis_values=cn_hl)
        ax.set_axis_off()
        for spine in ax.spines.values():
            spine.set_visible(False)
        if scalebar_um:
            bar_px = scalebar_um / PIXEL_SIZE_UM
            x0 = xmax - bar_px - (xmax - xmin) * 0.05
            y0 = ymax - (ymax - ymin) * 0.06
            # 线是矢量路径,文字是真文字(pdf.fonttype=42 / svg.fonttype=none)
            ax.plot([x0, x0 + bar_px], [y0, y0], color="black", lw=1.1,
                    solid_capstyle="butt", zorder=20)
            ax.text(x0 + bar_px / 2, y0 - (ymax - ymin) * 0.015,
                    f"{scalebar_um:g} µm", ha="center", va="bottom",
                    fontsize=5.5, color="black", zorder=20)
        for ext in ("pdf", "svg"):
            path = out_dir / f"{stem}_{kind}.{ext}"
            fig.savefig(path, dpi=dpi, transparent=True, pad_inches=0)
            written[f"{kind}_{ext}"] = path
        plt.close(fig)

    if legend:
        if legend_scope == "roi":
            cn_labels = sorted(sub["CN_k10"].unique())
            ct_labels = sorted(sub["celltype"].unique())
            note = " (present in ROI)"
        else:   # "all" —— 固定顺序、固定配色,多个面板共用一份
            cn_labels = m.CN_ORDER
            ct_labels = [c for c in m.CELLTYPE_COLORS]
            note = ""
        fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.2))
        for ax, title, labels, palette, with_niche in [
            (axes[0], "CN", cn_labels, m.CN_COLORS, cn_niche_names),
            (axes[1], "Phenotype", ct_labels, m.CELLTYPE_COLORS, False),
        ]:
            handles = [Line2D([], [], marker="o", ls="", ms=6, color=palette.get(l, "#888888"),
                              label=display_name(l, with_niche=with_niche)) for l in labels]
            ax.legend(handles=handles, title=f"{title}{note}", loc="center",
                      frameon=False, fontsize=8)
            ax.axis("off")
        path = out_dir / f"{stem}_legend.pdf"
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        written["legend_pdf"] = path

    print(f"{sample} {roi} — {len(sub):,} cells")
    for key, path in written.items():
        print(f"  {key:<16} {path.stat().st_size/1024:8.1f} KB  {path}")
    return written


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("sample")
    ap.add_argument("--by", default="CN_k10", choices=["CN_k10", "celltype"])
    ap.add_argument("--highlight", nargs="*", default=None, help="只给这些类别上色")
    ap.add_argument("--roi", nargs=4, type=float, default=None,
                    metavar=("XMIN", "XMAX", "YMIN", "YMAX"))
    ap.add_argument("--no-boxes", action="store_true")
    ap.add_argument("--focus-cn", dest="focus_cn", nargs="*", default=None,
                    help="本面板要突出的 CN,决定两行的重点细胞")
    ap.add_argument("--export", action="store_true",
                    help="配合 --roi:逐面板导出矢量 PDF/SVG + legend,给 AI 组版")
    args = ap.parse_args()

    if args.roi and args.export:
        export_crop_panels(args.sample, args.roi, focus_cn=args.focus_cn or ["CN1", "CN7"])
    elif args.roi:
        preview_crop(args.sample, args.roi, focus_cn=args.focus_cn or ["CN1", "CN7"])
    else:
        overview(args.sample, by=args.by, highlight=args.highlight, boxes=not args.no_boxes)


if __name__ == "__main__":
    main()
