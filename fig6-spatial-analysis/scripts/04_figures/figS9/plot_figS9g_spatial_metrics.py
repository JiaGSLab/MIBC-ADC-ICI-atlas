"""Render the four specimen-level spatial metrics used in Supplementary Fig. 9g."""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
import numpy as np
import pandas as pd

ROOT = Path("/data/home/chenwenhui/fig6")
SCRIPTS = Path(__file__).resolve().parents[1] / "fig6"
OUT = Path(os.environ.get("FIG6_FIGURE_OUTPUT", ROOT / "figures_vector" / "supplementary" / "figS9g_final"))
FONT_DIR = OUT / "fonts"
PROJECT = Path("/data/home/chenwenhui/mIHC/fig6_pm")


def load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / file)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def arial_8(fig):
    """Apply the requested typeface/size to every visible Matplotlib text artist."""
    for path in sorted(FONT_DIR.glob("*.ttf")):
        font_manager.fontManager.addfont(str(path))
    plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "pdf.fonttype": 42, "ps.fonttype": 42})
    for artist in fig.findobj(Text):
        if artist.get_visible():
            artist.set_fontfamily("Arial")
            artist.set_fontsize(8)


def save(fig, name: str, size: tuple[float, float]):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.set_size_inches(*size, forward=True)
    arial_8(fig)
    fig.savefig(OUT / name, format="pdf", facecolor="white")
    plt.close(fig)


def niche(rev):
    metrics = rev.build_sample_metrics()
    stats = [pd.read_csv(rev.DELIVERABLE_DIR / f"WP1.{i}_CRvsNR_stats.csv") for i in (4, 5)]
    wp14, wp15 = stats
    specs = [
        ("cd21b_same_type_knn_fraction", "CD21+B coherence", rev.stat_row_from_table(wp15, {"analysis_set": "primary_threshold30k_panelaware", "feature_class": "coherence", "metric": "same_type_knn_fraction", "feature_id": "CD21+ B cell"})),
        ("coupling_z_CN1_CN7", "CN1-CN7 coupling z", rev.stat_row_from_table(wp14, {"analysis_set": rev.PRIMARY_PROP_ANALYSIS, "feature_class": "adjacency", "metric": "z_score", "feature_id": "CN1__CN7"})),
        ("coupling_z_CN3_CN7", "CN3-CN7 coupling z", rev.stat_row_from_table(wp14, {"analysis_set": rev.PRIMARY_PROP_ANALYSIS, "feature_class": "adjacency", "metric": "z_score", "feature_id": "CN3__CN7"})),
        ("coupling_z_CN4_CN7", "CN4-CN7 coupling z", rev.stat_row_from_table(wp14, {"analysis_set": rev.PRIMARY_PROP_ANALYSIS, "feature_class": "adjacency", "metric": "z_score", "feature_id": "CN4__CN7"})),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.1))
    rng = np.random.default_rng(31)
    for ax, (key, label, stat) in zip(axes, specs):
        for xpos, group in enumerate(["CR", "NR"]):
            vals = metrics.loc[metrics["Outcome_short"].eq(group), key].dropna().astype(float).to_numpy()
            ax.scatter(np.full(len(vals), xpos) + rng.normal(0, .035, len(vals)), vals, color=rev.GROUP_COLORS[group], s=15, edgecolor="white", linewidth=.25, alpha=.82)
            if len(vals):
                mean = float(vals.mean()); boot = [float(np.mean(rng.choice(vals, len(vals), replace=True))) for _ in range(1000)]
                low, high = np.quantile(boot, [.025, .975]); ax.errorbar([xpos], [mean], yerr=[[mean-low], [high-mean]], fmt="_", color="black", ms=14, lw=1, capsize=2.5)
        ax.set_title(label); ax.set_xticks([0, 1]); ax.set_xticklabels(["CR", "NR"]); ax.grid(axis="y", color="#E7E7E7", linewidth=.6)
        ax.text(.5, .98, f"d = {float(stat['cohens_d_NR_minus_CR']):.2f}\np = {float(stat['wilcoxon_mannwhitney_p']):.3g}", transform=ax.transAxes, ha="center", va="top", bbox=dict(facecolor="white", edgecolor="#DDDDDD", alpha=.88, pad=2))
    axes[0].set_ylabel("Per-sample value"); fig.text(.01, .98, "g", weight="bold", va="top")
    fig.subplots_adjust(left=.055, right=.997, top=.77, bottom=.23, wspace=.50)
    save(fig, "figS9g_spatial_metrics.pdf", (7.2, 2.100))


def main():
    os.environ["FIG6_PROJECT_DIR"] = str(PROJECT)
    assembly = load("compact_assembly", "assemble_fig6_panels.py")
    revision = load("compact_revision", "build_fig6b_c_e_panel_sources.py")
    assembly.configure_matplotlib({**assembly.CONFIG, "font_family": "Arial"})
    niche(revision)


if __name__ == "__main__":
    main()
