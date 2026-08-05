#!/usr/bin/env python3
"""R-F Fig. 6 figure rework from frozen P1 and WP3 deliverables only.

No analysis is run here.  This script reads frozen CSV/parquet outputs, reuses the
existing WP3 spatial plotting helpers, and writes new canonical R-F exports.
Run in the scimap environment:
  python deliverables/R_F_figure_rework.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


PROJECT = Path("/data/home/chenwenhui/mIHC/fig6_pm")
DELIVERABLES = PROJECT / "deliverables"
CELLINFO = Path("/data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh")
REV_DIR = CELLINFO / "fig6_panels_rev"
SUPP_DIR = CELLINFO / "fig6_panels_supp"
P1 = DELIVERABLES / "P1"
CR, NR = "#E4212B", "#2C5F9E"
CN_COLORS = ["#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B", "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B"]


def configure() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.65,
        "axes.edgecolor": "#333333", "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white",
    })


def wp3_revision_module():
    path = Path(__file__).with_name("build_fig6b_c_e_panel_sources.py")
    if not path.exists():
        path = DELIVERABLES / "WP3_1b_fig6_revision.py"
    spec = importlib.util.spec_from_file_location("wp3_revision", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def tagged(df: pd.DataFrame, source_table: str, panel_part: str) -> pd.DataFrame:
    out = df.copy()
    out.insert(0, "panel_part", panel_part)
    out.insert(0, "source_table", source_table)
    return out


def combine(*frames: pd.DataFrame) -> pd.DataFrame:
    return pd.concat([f for f in frames if f is not None and not f.empty], ignore_index=True, sort=False)


def export(fig: plt.Figure, panel_dir: Path, basename: str, source: pd.DataFrame) -> dict[str, str]:
    panel_dir.mkdir(parents=True, exist_ok=True)
    paths = {ext: panel_dir / f"{basename}.{ext}" for ext in ("svg", "pdf", "png", "tiff")}
    source_path = panel_dir / f"{basename}_source_data.csv"
    source.to_csv(source_path, index=False)
    fig.savefig(paths["svg"], bbox_inches="tight", dpi=600)
    fig.savefig(paths["pdf"], bbox_inches="tight", dpi=600)
    fig.savefig(paths["png"], bbox_inches="tight", dpi=300)
    fig.savefig(paths["tiff"], bbox_inches="tight", dpi=600)
    plt.close(fig)
    return {**{key: str(value) for key, value in paths.items()}, "source_data": str(source_path)}


def source_path(relative: str) -> Path:
    return P1 / relative


def validate_frozen_inputs() -> pd.DataFrame:
    requirements = [
        ("6c", "R-A1_sample_balanced_cn/sample_cn5_metrics.csv"),
        ("6c", "../WP1.3_sample_proportion_long.csv"),
        ("6f", "R-A5b_composite_ladder/auc_summary.csv"),
        ("6f", "R-A5b_composite_ladder/paired_bootstrap_auc_differences.csv"),
        ("S11", "R-A1_sample_balanced_cn/center_matching_by_seed.csv"),
        ("S11", "R-A2_radius_sensitivity/functional_niche_matching.csv"),
        ("S11", "R-A2_radius_sensitivity/stability_verdicts.csv"),
        ("S12", "../WP2.2_leakage_comparison.csv"),
        ("S12", "R-A5_fair_auc_technical_controls/auc_summary.csv"),
        ("S12", "R-A5_fair_auc_technical_controls/control_availability.csv"),
        ("S13", "R-A1_sample_balanced_cn/headline_coupling_effects.csv"),
        ("S13", "R-A2_radius_sensitivity/outcome_statistics.csv"),
        ("S13", "R-A3_slice_registration_aware/complete/twoD_vs_slice_comparison.csv"),
        ("S13", "R-A4_stronger_coupling_null/headline_pair_group_results.csv"),
    ]
    rows = []
    for panel_id, relative in requirements:
        path = (P1 / relative).resolve()
        rows.append({"panel_id": panel_id, "source_table": str(path), "exists": path.exists(), "validated": path.exists() and path.stat().st_size > 0})
    return pd.DataFrame(rows)


def load_composite_ladder() -> tuple[pd.DataFrame, pd.DataFrame]:
    ladder = pd.read_csv(source_path("R-A5b_composite_ladder/auc_summary.csv"))
    order = ["celltype_proportion", "CN_proportion", "CN_coupling"]
    ladder = ladder.set_index("model_id").loc[order].reset_index()
    ladder["display_label"] = ["Cell-type proportions", "CN proportions", "Niche coupling"]
    ladder["color"] = ["#777777", "#4A7088", CR]
    differences = pd.read_csv(source_path("R-A5b_composite_ladder/paired_bootstrap_auc_differences.csv"))
    return ladder, differences


def leakage_labels() -> list[str]:
    return ["ISOS", "spatial logistic", "CN proportion", "cell type"]


def s13_column_labels() -> list[str]:
    return ["sample-\nbalanced", "75/100/125\num", "slice-\naware", "stronger\nnull"]


def select_6c_samples() -> pd.DataFrame:
    cn5 = pd.read_csv(source_path("R-A1_sample_balanced_cn/sample_cn5_metrics.csv"))
    frozen = cn5.loc[cn5["seed"].eq("frozen_reference")].copy()
    cr = frozen.loc[(frozen["sample"].eq("pt41-pre")) & frozen["CN5_proportion"].ge(0.01)].copy()
    if len(cr) != 1:
        raise RuntimeError("pt41-pre must remain the frozen CN5-positive CR representative")
    cr["selection_rule"] = "CN5 proportion >=1% with 1,578 CN5 cells and a contiguous lymphoid domain; selected over invalid pt63-pre (CN5 <1%)."
    nr = frozen.loc[frozen["sample"].eq("pt50-pre")].copy()
    if len(nr) != 1:
        raise RuntimeError("pt50-pre NR contrast not found")
    nr["selection_rule"] = "Sparse fragmented immune-infiltrate contrast; frozen reference has zero CN5 cells."
    selected = pd.concat([cr, nr], ignore_index=True)
    selected["selection_rank"] = [1, 2]
    return selected


def draw_bracket(ax, x: float, y1: float, y2: float, label: str) -> None:
    ax.plot([x, x + 0.012, x + 0.012, x], [y1, y1, y2, y2], color="#333333", lw=0.75, clip_on=False)
    ax.text(x + 0.018, (y1 + y2) / 2, label, va="center", fontsize=6.2)


def draw_ladder(ax, ladder: pd.DataFrame, differences: pd.DataFrame, compact: bool = False) -> None:
    ypos = np.arange(len(ladder))
    ax.axvline(0.5, color="#999999", ls="--", lw=0.8)
    for y, row in zip(ypos, ladder.itertuples(index=False)):
        ax.hlines(y, row.auc_bootstrap95_low, row.auc_bootstrap95_high, color=row.color, lw=1.9)
        ax.scatter(row.auc, y, s=30, color=row.color, zorder=3)
        ax.text(row.auc_bootstrap95_high + 0.015, y, f"{row.auc:.2f} [{row.auc_bootstrap95_low:.2f}, {row.auc_bootstrap95_high:.2f}]", va="center", fontsize=6.1)
    labels = ["Cell type", "CN proportion", "Niche coupling"] if compact else ladder["display_label"].tolist()
    ax.set_yticks(ypos, labels)
    ax.invert_yaxis()
    ax.set_xlim(0.20, 1.37)
    ax.set_xlabel("LOOCV AUC (bootstrap 95% CI)")
    ax.set_title("Fair fold-internal signed composite (<=3 features)", fontsize=7.5 if compact else 8.4)
    by_name = differences.set_index("comparison")
    draw_bracket(ax, 1.08, 0, 1, f"Delta=0.28; FDR={by_name.loc['CN_proportion - celltype_proportion', 'bootstrap_sign_p_fdr_bh']:.3f} **")
    draw_bracket(ax, 1.19, 0, 2, f"Delta=0.41; FDR={by_name.loc['CN_coupling - celltype_proportion', 'bootstrap_sign_p_fdr_bh']:.3f} **")
    draw_bracket(ax, 1.30, 1, 2, f"Delta=0.13; FDR={by_name.loc['CN_coupling - CN_proportion', 'bootstrap_sign_p_fdr_bh']:.3f} n.s.")


def panel_6f() -> tuple[dict[str, str], pd.DataFrame]:
    ladder, differences = load_composite_ladder()
    sensitivity = pd.read_csv(source_path("R-A5_fair_auc_technical_controls/auc_summary.csv"))
    sensitivity = sensitivity.loc[sensitivity["model_id"].isin(["celltype_proportion", "CN_proportion", "CN_coupling"])].copy()
    fig, ax = plt.subplots(figsize=(6.0, 2.8))
    draw_ladder(ax, ladder, differences)
    ax.text(0.01, -0.35, "Sensitivity (L2-logistic): cell type 0.49, CN proportion 0.68, coupling 0.80; same ordering.", transform=ax.transAxes, fontsize=6.2)
    fig.subplots_adjust(left=0.27, right=0.98, top=0.88, bottom=0.29)
    source = combine(tagged(ladder, "R-A5b_composite_ladder/auc_summary.csv", "fair_composite_ladder"), tagged(differences, "R-A5b_composite_ladder/paired_bootstrap_auc_differences.csv", "paired_bootstrap"), tagged(sensitivity, "R-A5_fair_auc_technical_controls/auc_summary.csv", "L2_logistic_sensitivity"))
    return export(fig, REV_DIR, "fig6f_fair_composite_ladder_revised", source), source


def panel_6c() -> tuple[dict[str, str], pd.DataFrame]:
    old = wp3_revision_module()
    selected = select_6c_samples()
    spatial = old.load_spatial_for_samples(selected["sample"].tolist())
    fig, axes = plt.subplots(2, 2, figsize=(5.8, 5.0), gridspec_kw={"hspace": 0.05, "wspace": 0.05})
    source_parts = [tagged(selected, "R-A1_sample_balanced_cn/sample_cn5_metrics.csv", "representative_selection")]
    for col, row in selected.reset_index(drop=True).iterrows():
        sample = row["sample"]
        sample_df = spatial.loc[spatial["sample"].eq(sample)].copy()
        spec = old.spatial_panel_specs("6c", row["Outcome_short"])
        crop = old.select_crop_window(sample_df, crop_size_px=int(spec["crop_size_px"]), cn_targets=list(spec["cn_targets"]), cell_targets=list(spec["cell_targets"]), cn_weights=dict(spec["cn_weights"]), cell_weights=dict(spec["cell_weights"]), require_first_two=bool(spec["require_first_two"]), interface_weight=float(spec["interface_weight"]))
        crop_df = old.crop_dataframe(sample_df, crop).copy()
        crop_df["panel_sample"] = sample
        crop_df["Outcome_short"] = row["Outcome_short"]
        crop_df["selection_rule"] = row["selection_rule"]
        crop_df["crop_xmin"], crop_df["crop_xmax"] = crop["xmin"], crop["xmax"]
        crop_df["crop_ymin"], crop_df["crop_ymax"] = crop["ymin"], crop["ymax"]
        source_parts.append(tagged(crop_df, "neighborhood_radius200_B0/neighborhood_composition.parquet", "spatial_crop_coordinates"))
        old.scatter_crop(axes[0, col], crop_df, color_col="celltype", highlights=set(spec["cell_highlights"]), palette=old.CELLTYPE_COLORS, crop=crop, point_size=1.65)
        old.scatter_crop(axes[1, col], crop_df, color_col="CN_k10", highlights=set(spec["cn_highlights"]), palette=old.CN_COLORS, crop=crop, point_size=1.65)
        old.add_locator(axes[0, col], sample_df, crop, CR if row["Outcome_short"] == "CR" else NR)
        old.draw_scale_bar(axes[1, col], crop)
        axes[0, col].set_title(f"{sample} ({row['Outcome_short']})", color=CR if row["Outcome_short"] == "CR" else NR, fontsize=8)
        axes[1, col].text(0.02, 0.04, f"CN5={row['CN5_proportion']:.2%}; n={int(row['CN5_cells'])}", transform=axes[1, col].transAxes, fontsize=6.1, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1})
    axes[0, 0].set_ylabel("Phenotype")
    axes[1, 0].set_ylabel("CN")
    fig.suptitle("6c. CR lymphoid domain versus sparse fragmented NR immune infiltrate", x=0.02, ha="left", fontsize=9, fontweight="bold")
    handles = [Line2D([0], [0], marker="o", color="w", label="CD21+ B cell", markerfacecolor=old.CELLTYPE_COLORS["CD21+ B cell"], markersize=5), Line2D([0], [0], marker="o", color="w", label="CN5 CD21+B/TLS-like", markerfacecolor=old.CN_COLORS["CN5"], markersize=5), Line2D([0], [0], marker="o", color="w", label="other cells", markerfacecolor="#D9D9D9", markersize=5)]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=6.2, bbox_to_anchor=(0.5, -0.01))
    fig.subplots_adjust(bottom=0.12, top=0.90, left=0.08, right=0.99)
    source = combine(*source_parts)
    return export(fig, REV_DIR, "fig6c_spatial_immune_tls_revised", source), source


def panel_s11() -> tuple[dict[str, str], pd.DataFrame]:
    k_match = pd.read_csv(DELIVERABLES / "WP1.2_k10_k11_matching.csv")
    seed = pd.read_csv(source_path("R-A1_sample_balanced_cn/center_matching_by_seed.csv"))
    radius = pd.read_csv(source_path("R-A2_radius_sensitivity/functional_niche_matching.csv"))
    verdict = pd.read_csv(source_path("R-A2_radius_sensitivity/stability_verdicts.csv"))
    seed_plot = seed.loc[~seed["seed"].eq("frozen_reference")].groupby("matched_current_CN", as_index=False)["cosine_similarity"].mean()
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 5.3), gridspec_kw={"hspace": 0.55, "wspace": 0.42})
    axes[0, 0].hist(k_match["best_cosine_similarity"], bins=np.linspace(0.9, 1, 15), color="#053E7A")
    axes[0, 0].set(title="Existing k10-vs-k11 centre matching", xlabel="best cosine", ylabel="clusters")
    x = np.arange(len(seed_plot)); colors = [CR if cn == "CN5" else CN_COLORS[int(cn[2:])] for cn in seed_plot["matched_current_CN"]]
    axes[0, 1].bar(x, seed_plot["cosine_similarity"], color=colors); axes[0, 1].set_xticks(x, seed_plot["matched_current_CN"]); axes[0, 1].set_ylim(0, 1.05); axes[0, 1].set(title="Sample-balanced centres (5 seeds)", ylabel="mean cosine")
    axes[0, 1].text(0.02, 0.04, "CN5 lowest: 0.51", transform=axes[0, 1].transAxes, fontsize=6.5)
    radius_plot = radius.loc[radius["radius_um"].isin([75.0, 125.0]) & radius["matched_frozen_CN"].str.match(r"^CN\d+$")].copy()
    for cn, sub in radius_plot.groupby("matched_frozen_CN"):
        axes[1, 0].plot(sub["radius_um"], sub["cosine_similarity"], marker="o", lw=0.7, color=CR if cn == "CN5" else "#999999", alpha=0.8)
    axes[1, 0].axhline(0.6, color="#777777", ls="--", lw=0.7); axes[1, 0].set(title="Radius sensitivity of functional niches", xlabel="radius (um)", ylabel="matched-centroid cosine", ylim=(0, 1.05))
    v = verdict.loc[verdict["finding"].str.startswith(("CN5_", "coupling_RoE_"))].copy(); y = np.arange(len(v)); axes[1, 1].scatter(np.where(v["verdict"].eq("preserved"), 1, 0), y, s=35, c=np.where(v["verdict"].eq("preserved"), "#2F4F4F", CR)); axes[1, 1].set_yticks(y, v["finding"].str.replace("coupling_RoE_", "", regex=False)); axes[1, 1].set_xticks([0, 1], ["not stable", "preserved"]); axes[1, 1].set(title="Outcome robustness verdicts")
    fig.suptitle("S11. CN clustering robustness", x=0.01, ha="left", fontsize=10, fontweight="bold")
    source = combine(tagged(k_match, "WP1.2_k10_k11_matching.csv", "existing_k10_k11"), tagged(seed, "R-A1_sample_balanced_cn/center_matching_by_seed.csv", "sample_balanced"), tagged(radius, "R-A2_radius_sensitivity/functional_niche_matching.csv", "radius_matching"), tagged(verdict, "R-A2_radius_sensitivity/stability_verdicts.csv", "stability_verdicts"))
    return export(fig, SUPP_DIR, "S11_cn_robustness_revised", source), source


def panel_s12() -> tuple[dict[str, str], pd.DataFrame]:
    leak = pd.read_csv(DELIVERABLES / "WP2.2_leakage_comparison.csv")
    ladder, differences = load_composite_ladder()
    controls = pd.read_csv(source_path("R-A5_fair_auc_technical_controls/auc_summary.csv"))
    availability = pd.read_csv(source_path("R-A5_fair_auc_technical_controls/control_availability.csv"))
    technical_pairs = pd.read_csv(source_path("R-A5_fair_auc_technical_controls/paired_bootstrap_auc_differences.csv"))
    fig, axes = plt.subplots(2, 2, figsize=(8.1, 5.4), gridspec_kw={"hspace": 0.58, "wspace": 0.48})
    leak_plot = leak.melt(id_vars="comparison", value_vars=["wp2_1_auc", "wp2_2_auc"], var_name="analysis", value_name="auc")
    for i, name in enumerate(leak["comparison"]):
        vals = leak_plot.loc[leak_plot["comparison"].eq(name), "auc"].to_numpy(); axes[0, 0].plot(vals, [i, i], marker="o", color="#555555")
    axes[0, 0].set(yticks=np.arange(len(leak)), yticklabels=leakage_labels(), xlim=(0, 1), xlabel="AUC", title="Leakage demo retained")
    draw_ladder(axes[0, 1], ladder, differences, compact=True)
    tech = controls.loc[controls["model_id"].str.startswith("technical_")].copy(); labels = ["cell count", "defined fraction", "panel version"]
    axes[1, 0].bar(np.arange(len(tech)), tech["auc"], color=["#556B2F", "#4A7088", "#8B5A2B"]); axes[1, 0].set(xticks=np.arange(len(tech)), xticklabels=labels, ylim=(0, 1), ylabel="AUC", title="Technical / panel controls")
    comp = technical_pairs.loc[technical_pairs["model_left"].eq("CN_coupling")].copy(); y = np.arange(len(comp)); axes[1, 1].hlines(y, comp["bootstrap95_low"], comp["bootstrap95_high"], color="#555555"); axes[1, 1].scatter(comp["auc_difference_left_minus_right"], y, color=CR); axes[1, 1].axvline(0, color="#999999", ls="--", lw=0.7); axes[1, 1].set(yticks=y, yticklabels=comp["model_right"].str.replace("technical_", "", regex=False), xlabel="coupling minus comparator AUC", title="Coupling is not a technical proxy")
    fig.text(0.51, 0.02, "Tissue area, HER2-IHC, and clinical covariates unavailable and not proxied.", fontsize=6.2)
    fig.suptitle("S12. Prediction integrity and technical controls", x=0.01, ha="left", fontsize=10, fontweight="bold")
    source = combine(tagged(leak, "WP2.2_leakage_comparison.csv", "leakage_demo"), tagged(ladder, "R-A5b_composite_ladder/auc_summary.csv", "fair_composite_ladder"), tagged(differences, "R-A5b_composite_ladder/paired_bootstrap_auc_differences.csv", "fair_paired_bootstrap"), tagged(controls, "R-A5_fair_auc_technical_controls/auc_summary.csv", "technical_controls"), tagged(technical_pairs, "R-A5_fair_auc_technical_controls/paired_bootstrap_auc_differences.csv", "technical_paired_bootstrap"), tagged(availability, "R-A5_fair_auc_technical_controls/control_availability.csv", "control_availability"))
    return export(fig, SUPP_DIR, "S12_prediction_integrity_revised", source), source


def panel_s13() -> tuple[dict[str, str], pd.DataFrame]:
    a1 = pd.read_csv(source_path("R-A1_sample_balanced_cn/headline_coupling_effects.csv"))
    a2 = pd.read_csv(source_path("R-A2_radius_sensitivity/outcome_statistics.csv"))
    a3 = pd.read_csv(source_path("R-A3_slice_registration_aware/complete/twoD_vs_slice_comparison.csv"))
    a4 = pd.read_csv(source_path("R-A4_stronger_coupling_null/headline_pair_group_results.csv"))
    pairs = ["CN1__CN7", "CN3__CN7", "CN4__CN7"]
    labels = ["CN1-CN7 (CR)", "CN3-CN7 (NR)", "CN4-CN7 (CR)"]
    expected = np.array([-1, 1, -1])  # sign of NR-minus-CR Cohen d expected by the protocol
    a1m = a1.loc[(a1["seed"].ne("frozen_reference")) & a1["feature"].isin(pairs)].groupby("feature")["cohens_d_NR_minus_CR"].median().reindex(pairs).to_numpy()
    a2m = a2.loc[a2["endpoint"].str.startswith("coupling_RoE_")].assign(pair=lambda d: d["endpoint"].str.replace("coupling_RoE_", "", regex=False)).groupby("pair")["effect_size_cohens_d_NR_minus_CR"].median().reindex(pairs).to_numpy()
    a3m = a3.groupby("coupling_pair")["same_sign_fraction"].median().reindex(pairs).to_numpy()
    a4s = a4.loc[a4["null_method"].eq("degree_blocksize_stratified")].set_index("CN_pair").reindex(pairs)["cohens_d_RoE_NR_minus_CR"].to_numpy()
    matrix = np.vstack([np.sign(a1m) * expected, np.sign(a2m) * expected, a3m, np.sign(a4s) * expected]).T
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.65), gridspec_kw={"wspace": 0.55})
    im = axes[0].imshow(matrix, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
    axes[0].set(xticks=np.arange(4), xticklabels=s13_column_labels(), yticks=np.arange(3), yticklabels=labels, title="Direction preserved across perturbations")
    for i in range(3):
        for j in range(4): axes[0].text(j, i, "preserved" if matrix[i, j] > 0.5 else "weak", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, label="agreement")
    plot = a4.loc[a4["null_method"].isin(["global_label_permutation", "degree_blocksize_stratified"])].copy(); plot["pair_label"] = plot["CN_pair"].map(dict(zip(pairs, labels))); plot["color"] = np.where(plot["CN_pair"].eq("CN4__CN7"), "#E08A3C", CR)
    for i, (_, row) in enumerate(plot.iterrows()):
        axes[1].plot([row["cohens_d_bootstrap95_low"], row["cohens_d_bootstrap95_high"]], [i, i], color="#444444", lw=1); axes[1].scatter(row["cohens_d_RoE_NR_minus_CR"], i, c=row["color"], s=28); axes[1].text(0.03, i + 0.18, f"{row['null_method'].replace('_', ' ')}; FDR={row['FDR_BH_group_empirical_p']:.3f}", fontsize=5.8)
    axes[1].axvline(0, color="#999999", ls="--", lw=0.7); axes[1].set(yticks=np.arange(len(plot)), yticklabels=plot["pair_label"], xlabel="Cohen d (NR - CR)", title="R-A4 stronger null (3 pre-specified pairs)")
    fig.text(0.01, 0.01, "R-A4 FDR is confirmatory within the three pre-specified pairs; cohort-wide discovery FDR: none <0.10. CD8/CN4 passes stratified null but not global (p=0.071), the weakest result.", fontsize=6.1)
    fig.suptitle("S13. Spatial-coupling robustness summary", x=0.01, ha="left", fontsize=10, fontweight="bold")
    source = combine(tagged(a1, "R-A1_sample_balanced_cn/headline_coupling_effects.csv", "sample_balanced"), tagged(a2, "R-A2_radius_sensitivity/outcome_statistics.csv", "radius_sensitivity"), tagged(a3, "R-A3_slice_registration_aware/complete/twoD_vs_slice_comparison.csv", "slice_aware"), tagged(a4, "R-A4_stronger_coupling_null/headline_pair_group_results.csv", "stronger_null"))
    return export(fig, SUPP_DIR, "S13_spatial_coupling_robustness", source), source


def main_assembly(revised_6c: dict[str, str], revised_6f: dict[str, str]) -> tuple[dict[str, str], pd.DataFrame]:
    panels = [
        ("6a", REV_DIR / "fig6a_workflow_schematic_placeholder.png"), ("6b", REV_DIR / "fig6b_cn_composition_heatmap.png"),
        ("6c", Path(revised_6c["png"])), ("6d", REV_DIR / "fig6d_spatial_vascular_axis_crop_candidates.png"),
        ("6e", REV_DIR / "fig6e_tls_presence_and_vascular_coupling.png"), ("6f", Path(revised_6f["png"])), ("6g", REV_DIR / "fig6g_prediction_roc_nested.png"),
    ]
    fig, axes = plt.subplots(4, 2, figsize=(13, 14), constrained_layout=True)
    records = []
    for ax, (panel, path) in zip(axes.flat, panels):
        ax.imshow(mpimg.imread(path)); ax.axis("off"); ax.set_title(panel, loc="left", fontweight="bold")
        records.append({"panel": panel, "rendered_from_png": str(path), "revised": panel in {"6c", "6f"}})
    axes.flat[-1].axis("off")
    source = pd.DataFrame(records)
    return export(fig, REV_DIR, "fig6_main_assembly_revised", source), source


def run() -> dict[str, object]:
    configure()
    checks = validate_frozen_inputs()
    if not checks["validated"].all(): raise RuntimeError("Missing frozen input(s)")
    out6c, source6c = panel_6c()
    out6f, source6f = panel_6f()
    out11, source11 = panel_s11()
    out12, source12 = panel_s12()
    out13, source13 = panel_s13()
    assembly, assembly_source = main_assembly(out6c, out6f)
    outputs = [("6c", "fig6c_spatial_immune_tls_revised", out6c, len(source6c)), ("6f", "fig6f_fair_composite_ladder_revised", out6f, len(source6f)), ("S11", "S11_cn_robustness_revised", out11, len(source11)), ("S12", "S12_prediction_integrity_revised", out12, len(source12)), ("S13", "S13_spatial_coupling_robustness", out13, len(source13)), ("assembly", "fig6_main_assembly_revised", assembly, len(assembly_source))]
    index_rows = []
    for panel_id, basename, files, rows in outputs:
        index_rows.append({"panel_id": panel_id, "basename": basename, "source_rows": rows, **files})
    index = pd.DataFrame(index_rows)
    index_path = DELIVERABLES / "R_F_figure_rework_index.csv"; index.to_csv(index_path, index=False)
    audit = []
    for _, row in index.iterrows():
        for kind in ["svg", "pdf", "png", "tiff", "source_data"]:
            p = Path(row[kind]); audit.append({"panel_id": row["panel_id"], "kind": kind, "path": str(p), "exists": p.exists(), "bytes": p.stat().st_size if p.exists() else 0})
    audit_path = DELIVERABLES / "R_F_figure_rework_export_audit.csv"; pd.DataFrame(audit).to_csv(audit_path, index=False)
    metadata = {"work_package": "R-F", "backend": "Python/matplotlib", "analysis_rerun": False, "frozen_input_validation": str(DELIVERABLES / "R_F_frozen_input_validation.csv"), "index": str(index_path), "audit": str(audit_path), "canonical_revised_outputs": index[["panel_id", "basename", "png", "pdf", "source_data"]].to_dict(orient="records"), "representative_selection": "pt41-pre CR (CN5=1.58%, 1,578 cells); pt50-pre NR contrast; pt63-pre excluded because CN5 <1%."}
    checks.to_csv(DELIVERABLES / "R_F_frozen_input_validation.csv", index=False)
    (DELIVERABLES / "R_F_figure_rework_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
