#!/usr/bin/env python3
"""Fig. 6a workflow schematic — BioRender-style flat vector, drawn in matplotlib.

替代 `fig6a_workflow_schematic_placeholder.*`。所有元素均为矢量路径,
导出 PDF/SVG 后可在 Illustrator 中逐元素编辑(字体已 outline-free, svg.fonttype='none')。

流程与 fig6_review_gallery.html 中 6a 的描述一致:
  3 切片 + 3 panel -> registration -> section-specific cells
  -> radius200 neighbourhood -> CN -> coupling -> prediction

Run:
  /data/home/chenwenhui/miniconda3/envs/scimap/bin/python fig6a_workflow_schematic.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, PathPatch, Polygon
from matplotlib.path import Path as MplPath

OUT_DIR = Path("/data/home/chenwenhui/fig6/figures_vector/main")
STEM = "fig6a_workflow_schematic"

# 与 WP3/R-F 一致的手稿配色
CR, NR = "#E4212B", "#2C5F9E"
CELLTYPE_COLORS = {
    "B cell": "#A52A2A", "CD163+ Macrophage": "#053E7A", "CD163- Macrophage": "#E46B4F",
    "CD21+ B cell": "#8B8B00", "CD4 T cell": "#B8860B", "CD8 T cell": "#CD661D",
    "CD8 exhausted T cell": "#8B5A2B", "DC": "#CDAD00", "Vascular": "#2F4F4F",
    "Fibroblast / SMA+ stromal": "#483D8B", "HER2+ Tumor": "#7A378B", "Treg": "#4A7088",
    "PanCK+ Tumor": "#CD7054",
}
CN_COLORS = ["#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",
             "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B"]

# 每个 stage 的 (浅底色, 描边色)
STAGE_TINTS = [
    ("#E9EFF7", "#3A6FA8"),
    ("#E9EFF7", "#3A6FA8"),
    ("#F1ECF7", "#7A5AA0"),
    ("#E9F3EE", "#3E8E6E"),
    ("#FDF2E7", "#CE8433"),
    ("#FBEBEC", "#C0383F"),
]

FIG_W, FIG_H = 7.48, 2.05          # 与占位图同宽(538.92 pt),略增高以容纳图标
XMAX = 100.0
YMAX = XMAX * FIG_H / FIG_W        # 保持单位为正方形,圆不变形

BOX_W, BOX_GAP = 14.9, 1.8
BOX_X0 = 0.8
BOX_Y0, BOX_Y1 = 2.6, 22.6
ICON_CY = 12.9                     # 图标区中心
TITLE_Y = 21.0
CAP_Y = 6.6


def rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(seed)


def blob(cx, cy, rx, ry, seed=0, n=90, wobble=0.10):
    """柔和的组织轮廓(BioRender 里那种不规则圆)。"""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    g = rng(seed)
    r = 1.0 + wobble * np.sin(3 * t + g.uniform(0, 6)) + 0.6 * wobble * np.sin(5 * t + g.uniform(0, 6))
    return np.column_stack([cx + rx * r * np.cos(t), cy + ry * r * np.sin(t)])


def add_blob(ax, cx, cy, rx, ry, seed, fc, ec, lw=0.7, alpha=1.0, ls="-", z=2):
    poly = Polygon(blob(cx, cy, rx, ry, seed), closed=True, facecolor=fc, edgecolor=ec,
                   linewidth=lw, alpha=alpha, linestyle=ls, zorder=z, joinstyle="round")
    ax.add_patch(poly)
    return poly


def points_in_blob(cx, cy, rx, ry, seed, n):
    """在椭圆内均匀采样(用作细胞点)。"""
    g = rng(seed)
    r = np.sqrt(g.uniform(0, 1, n))
    t = g.uniform(0, 2 * np.pi, n)
    return cx + rx * r * np.cos(t), cy + ry * r * np.sin(t)


def stage_box(ax, i):
    x0 = BOX_X0 + i * (BOX_W + BOX_GAP)
    fc, ec = STAGE_TINTS[i]
    ax.add_patch(FancyBboxPatch(
        (x0, BOX_Y0), BOX_W, BOX_Y1 - BOX_Y0,
        boxstyle="round,pad=0,rounding_size=1.4",
        facecolor=fc, edgecolor=ec, linewidth=0.6, zorder=1))
    return x0


def stage_text(ax, x0, title, caption, ec):
    cx = x0 + BOX_W / 2
    ax.text(cx, TITLE_Y, title, ha="center", va="top", fontsize=5.6, fontweight="bold",
            color=ec, linespacing=1.28, zorder=6)
    ax.text(cx, CAP_Y, caption, ha="center", va="top", fontsize=4.5, color="#4A4A4A",
            linespacing=1.35, zorder=6)


def chevron(ax, x_from, x_to, y):
    ax.add_patch(FancyArrowPatch(
        (x_from, y), (x_to, y), arrowstyle="-|>", mutation_scale=5.5,
        linewidth=0.9, color="#8A8F98", shrinkA=0, shrinkB=0, zorder=5))


# --------------------------------------------------------------------------- #
# 各 stage 图标
# --------------------------------------------------------------------------- #
def icon_sections(ax, cx, cy):
    """3 张连续切片 + 每张 3 个 panel 的色点。"""
    panel_cols = ["#3A6FA8", "#5FA88F", "#D08A3C"]
    for k, dx in enumerate([-2.5, 0.0, 2.5]):
        dy = 1.6 - 1.6 * k
        ax.add_patch(FancyBboxPatch(
            (cx + dx - 3.1, cy + dy - 2.35), 6.2, 4.7,
            boxstyle="round,pad=0,rounding_size=0.45",
            facecolor="#FFFFFF", edgecolor="#9AA3AE", linewidth=0.55, zorder=3 + k))
        add_blob(ax, cx + dx, cy + dy, 2.0, 1.35, seed=11 + k,
                 fc="#D7E3F0", ec="#7E93AC", lw=0.5, z=3.5 + k)
        for j, col in enumerate(panel_cols):
            ax.add_patch(Circle((cx + dx - 2.2 + 0.85 * j, cy + dy + 1.72), 0.28,
                                facecolor=col, edgecolor="none", zorder=4 + k))


def icon_registration(ax, cx, cy):
    """错位轮廓 -> 对齐叠合。"""
    add_blob(ax, cx - 3.4, cy + 0.55, 2.0, 1.5, seed=21, fc="none", ec="#B2BCC8", lw=0.7, ls=(0, (2, 1.6)))
    add_blob(ax, cx - 2.6, cy - 0.65, 2.0, 1.5, seed=22, fc="none", ec="#B2BCC8", lw=0.7, ls=(0, (2, 1.6)))
    ax.add_patch(FancyArrowPatch((cx - 0.55, cy), (cx + 0.85, cy), arrowstyle="-|>",
                                 mutation_scale=4.5, linewidth=0.75, color="#3A6FA8", zorder=5))
    add_blob(ax, cx + 3.2, cy, 2.15, 1.6, seed=21, fc="#D7E3F0", ec="#3A6FA8", lw=0.7, z=4)
    add_blob(ax, cx + 3.2, cy, 2.15, 1.6, seed=22, fc="none", ec="#C0383F", lw=0.6, ls=(0, (2.2, 1.4)), z=5)
    # 十字标记
    for (mx, my) in [(cx + 2.0, cy + 1.0), (cx + 4.4, cy - 0.9)]:
        ax.plot([mx - 0.45, mx + 0.45], [my, my], lw=0.55, color="#3A6FA8", zorder=6)
        ax.plot([mx, mx], [my - 0.45, my + 0.45], lw=0.55, color="#3A6FA8", zorder=6)


def icon_cells(ax, cx, cy):
    """分割 + crosstalk 校正后的细胞表型点。"""
    add_blob(ax, cx, cy, 4.5, 2.9, seed=31, fc="#FFFFFF", ec="#9C86B8", lw=0.7, z=3)
    cols = list(CELLTYPE_COLORS.values())
    xs, ys = points_in_blob(cx, cy, 3.8, 2.35, seed=32, n=110)
    g = rng(33)
    ax.scatter(xs, ys, s=1.9, c=[cols[i] for i in g.integers(0, len(cols), xs.size)],
               linewidths=0, zorder=4)


def icon_neighbourhood(ax, cx, cy):
    """r = 200 µm 邻域 -> CN 色块。"""
    lx = cx - 3.0
    cols = CN_COLORS
    xs, ys = points_in_blob(lx, cy, 3.3, 2.6, seed=41, n=95)
    g = rng(42)
    ax.scatter(xs, ys, s=1.5, c=[cols[i] for i in g.integers(0, 6, xs.size)],
               linewidths=0, alpha=0.55, zorder=3)
    ax.add_patch(Circle((lx, cy), 2.0, facecolor="none", edgecolor="#2C6B52",
                        linewidth=0.75, linestyle=(0, (2.2, 1.5)), zorder=5))
    ax.add_patch(Circle((lx, cy), 0.42, facecolor="#2C6B52", edgecolor="white",
                        linewidth=0.35, zorder=6))
    ax.annotate("", xy=(lx + 2.0, cy), xytext=(lx, cy),
                arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#2C6B52",
                                mutation_scale=4, shrinkA=0, shrinkB=0), zorder=6)
    ax.text(lx, cy + 2.25, "r = 200 µm", ha="center", va="bottom", fontsize=3.9,
            color="#2C6B52", zorder=7)
    # 右侧: k=10 CN 拼块
    ax.annotate("", xy=(cx + 1.35, cy), xytext=(cx + 0.35, cy),
                arrowprops=dict(arrowstyle="-|>", lw=0.7, color="#8A8F98",
                                mutation_scale=4.5), zorder=5)
    for j in range(10):
        r, c = divmod(j, 5)
        ax.add_patch(FancyBboxPatch(
            (cx + 1.9 + c * 1.02, cy + 0.15 - r * 1.02), 0.88, 0.88,
            boxstyle="round,pad=0,rounding_size=0.2",
            facecolor=CN_COLORS[j], edgecolor="white", linewidth=0.35, zorder=4))
    ax.text(cx + 4.45, cy - 1.35, "k = 10", ha="center", va="top", fontsize=4.0,
            color="#2C6B52", zorder=6)


def icon_coupling(ax, cx, cy):
    """成对(二元)的 CN block 接触:两块 territory 相接,一条 pair 边。

    coupling 定义在 CN_i x CN_j 上(WP1.4 的 CN_pair / z_score),不是多元团簇,
    所以这里只画一条二元边,不画三角形闭环。
    """
    terr = [(cx - 2.05, cy + 0.75, 2.5, 1.85, CN_COLORS[1], 51),
            (cx + 2.05, cy - 0.75, 2.5, 1.85, CN_COLORS[7], 52)]
    for tx, ty, rx, ry, col, sd in terr:
        add_blob(ax, tx, ty, rx, ry, seed=sd, fc=col, ec="white", lw=0.6, alpha=0.30, z=3)
        add_blob(ax, tx, ty, rx, ry, seed=sd, fc="none", ec=col, lw=0.65, z=4)
    (ax0, ay0), (ax1, ay1) = [(t[0], t[1]) for t in terr]
    ax.plot([ax0, ax1], [ay0, ay1], lw=0.9, color="#CE8433", zorder=5, solid_capstyle="round")
    for (nx, ny), col in zip([(ax0, ay0), (ax1, ay1)], [CN_COLORS[1], CN_COLORS[7]]):
        ax.add_patch(Circle((nx, ny), 0.46, facecolor=col, edgecolor="white",
                            linewidth=0.45, zorder=6))
    ax.text(cx - 2.05, cy + 2.95, "CN$_i$", ha="center", va="bottom", fontsize=4.2,
            color="#7A6030", zorder=6)
    ax.text(cx + 2.05, cy - 3.05, "CN$_j$", ha="center", va="top", fontsize=4.2,
            color="#7A6030", zorder=6)
    ax.text(cx, cy, "pair $z$", ha="center", va="center", fontsize=4.0,
            color="#CE8433", zorder=7,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.9, alpha=0.92))


def icon_prediction(ax, cx, cy):
    """嵌套 CV 的 ROC 曲线。"""
    w, h = 5.4, 4.4
    x0, y0 = cx - w / 2, cy - h / 2 - 0.2
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0,rounding_size=0.3",
                                facecolor="white", edgecolor="#C6CBD2", linewidth=0.5, zorder=3))
    ax.plot([x0, x0 + w], [y0, y0 + h], lw=0.5, color="#B7BCC4",
            linestyle=(0, (2, 1.5)), zorder=4)
    t = np.linspace(0, 1, 60)
    roc = 1 - (1 - t) ** 2.75
    ax.plot(x0 + t * w, y0 + roc * h, lw=1.0, color=CR, zorder=5, solid_capstyle="round")
    ax.text(cx, y0 - 0.55, "AUC 0.845", ha="center", va="top", fontsize=4.2,
            color=CR, fontweight="bold", zorder=6)
    for lbl, col, dx in [("CR", CR, -1.9), ("NR", NR, 0.55)]:
        ax.add_patch(Circle((cx + dx, y0 + h + 0.85), 0.3, facecolor=col,
                            edgecolor="none", zorder=6))
        ax.text(cx + dx + 0.5, y0 + h + 0.85, lbl, ha="left", va="center",
                fontsize=4.2, color="#3A3A3A", zorder=6)


STAGES = [
    dict(title="Serial sections\n× 3 mIHC panels",
         caption="48 patients · 82 sections\n3 consecutive slides",
         icon=icon_sections),
    dict(title="Cross-panel\nregistration",
         caption="VALIS rigid + non-rigid\n(no micro-registration)",
         icon=icon_registration),
    dict(title="Section-specific\ncells",
         caption="crosstalk-corrected\nphenotyping · 13 cell types\n23.1M cells",
         icon=icon_cells),
    dict(title="Radius-200\nneighbourhoods",
         caption="cell composition within\n200 µm · k-means\n10 CN territories",
         icon=icon_neighbourhood),
    dict(title="Pairwise CN\ncoupling",
         caption="Delaunay block adjacency\nper CN$_i$–CN$_j$ pair\n→ contact rate · RoE · z",
         icon=icon_coupling),
    dict(title="Response\nprediction",
         caption="pre-treatment CR vs NR\nfully nested CV",
         icon=icon_prediction),
]


def build() -> plt.Figure:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
        "pdf.fonttype": 42, "ps.fonttype": 42,
        # mathtext(CN$_i$ / $z$)也走同一套无衬线字体,避免 PDF 里混进 DejaVu
        "mathtext.fontset": "custom", "mathtext.rm": "sans", "mathtext.it": "sans:italic",
        "mathtext.bf": "sans:bold", "mathtext.default": "regular",
        "svg.fonttype": "none", "figure.facecolor": "white",
    })
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, XMAX)
    ax.set_ylim(0, YMAX)
    ax.set_axis_off()

    ax.text(0.8, YMAX - 0.6, "a", ha="left", va="top", fontsize=9, fontweight="bold",
            color="#111111", zorder=10)

    for i, st in enumerate(STAGES):
        x0 = stage_box(ax, i)
        cx = x0 + BOX_W / 2
        st["icon"](ax, cx, ICON_CY)
        stage_text(ax, x0, st["title"], st["caption"], STAGE_TINTS[i][1])
        if i < len(STAGES) - 1:
            chevron(ax, x0 + BOX_W + 0.35, x0 + BOX_W + BOX_GAP - 0.35, ICON_CY)
    return fig


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = build()
    for ext, dpi in [("pdf", 600), ("svg", 600), ("png", 400)]:
        fig.savefig(OUT_DIR / f"{STEM}.{ext}", dpi=dpi, facecolor="white")
        print("wrote", OUT_DIR / f"{STEM}.{ext}")
    plt.close(fig)


if __name__ == "__main__":
    main()
