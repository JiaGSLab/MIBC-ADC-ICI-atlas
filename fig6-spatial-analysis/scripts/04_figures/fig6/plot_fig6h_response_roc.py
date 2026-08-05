"""Auditable comparative ROC: P1 fold-internal composite versus three singles.

Conclusion: the held-out composite score can be compared fairly with each
individual coupling score only when every score uses its own training-fold
direction, imputation, and scaling.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
import numpy as np
import pandas as pd
from sklearn.metrics import auc, roc_auc_score, roc_curve

PROJECT = Path("/data/home/chenwenhui/mIHC/fig6_pm")
P1 = PROJECT / "deliverables" / "P1" / "R-A5b_composite_ladder"
OUT = Path("/data/home/chenwenhui/fig6/figures_vector/main/candidate_6g_comparative_roc_ci_20260731")
FONT_DIR = Path("/data/home/chenwenhui/fig6/figures_vector/main/candidate_compact_20260728/fonts")


def protocol_module():
    path = P1 / "R_A5b_composite_ladder.py"
    spec = importlib.util.spec_from_file_location("ra5b_protocol", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def apply_ariel_8(fig):
    for path in FONT_DIR.glob("*.ttf"):
        font_manager.fontManager.addfont(str(path))
    plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
    for artist in fig.findobj(Text):
        if artist.get_visible():
            artist.set_fontfamily("Arial")
            artist.set_fontsize(8)


def main():
    protocol = protocol_module()
    models, y, _cohort, _manifest = protocol.build_inputs()
    requested = {
        "CN1-CN7 coupling": "CN_coupling__CN1__CN7",
        "CN3-CN7 coupling": "CN_coupling__CN3__CN7",
        "CN4-CN7 coupling": "CN_coupling__CN4__CN7",
    }
    missing = set(requested.values()).difference(models["CN_coupling"].columns)
    if missing:
        raise ValueError(f"Requested coupling predictors missing from locked P1 inputs: {sorted(missing)}")

    composite = pd.read_csv(P1 / "fold_predictions.csv")
    composite = composite.loc[composite["model_id"].eq("CN_coupling")].copy()
    composite["model_label"] = "Composite coupling"
    # Reuse the protocol's exact LOOCV function with one candidate column.  This
    # preserves fold-specific CR-minus-NR directions and keeps only held-out scores.
    singles = []
    directions = []
    for label, feature in requested.items():
        pred, _selected, direct = protocol.loocv_signed_composite(models["CN_coupling"][[feature]], y, feature, cap=1)
        pred["model_label"] = label
        singles.append(pred)
        direct["model_label"] = label
        directions.append(direct)
    scores = pd.concat([composite, *singles], ignore_index=True, sort=False)
    scores["score_kind"] = "held_out_fold_score"
    score_labels = {"CN_coupling": "Composite coupling", **{feature: label for label, feature in requested.items()}}
    summary_rows, curve_rows = [], []
    palette = {"Composite coupling": "#E4212B", "CN1-CN7 coupling": "#1F77B4", "CN3-CN7 coupling": "#2F7D4A", "CN4-CN7 coupling": "#9467BD"}
    for model_id, label in score_labels.items():
        part = scores.loc[scores["model_id"].eq(model_id)].sort_values("sample")
        if len(part) != len(y) or not part["sample"].is_unique:
            raise ValueError(f"{label} does not contain exactly one held-out score per locked sample")
        truth = part.set_index("sample").loc[y.index, "truth"].astype(int)
        heldout = part.set_index("sample").loc[y.index, "score"].astype(float)
        if not truth.equals(y):
            raise ValueError(f"{label} truth labels do not match the locked P1 cohort")
        observed = float(roc_auc_score(truth, heldout))
        fpr, tpr, _ = roc_curve(truth, heldout)
        summary_rows.append({"model_id": model_id, "model_label": label, "auc": observed, "n_held_out_scores": len(part), "protocol": "P1 LOOCV training-fold direction, imputation, scaling; held-out score only"})
        curve_rows.extend({"model_id": model_id, "model_label": label, "fpr": x, "tpr": z, "auc": observed} for x, z in zip(fpr, tpr))
    summary = pd.DataFrame(summary_rows)
    reference = pd.read_csv(P1 / "auc_summary.csv").set_index("model_id")
    composite_auc = float(summary.loc[summary.model_id.eq("CN_coupling"), "auc"].iloc[0])
    if abs(composite_auc - float(reference.loc["CN_coupling", "auc"])) > 1e-12:
        raise AssertionError("Composite fold-score AUC does not reproduce frozen P1 AUC")

    # The composite alone receives the same deterministic held-out ROC bootstrap
    # treatment used by the current Fig. 6g: 2,000 sample-level resamples with
    # fixed seed 20260708; single-feature curves deliberately receive no band.
    composite_part = scores.loc[scores["model_id"].eq("CN_coupling")].set_index("sample").loc[y.index]
    truth = composite_part["truth"].to_numpy(int)
    heldout = composite_part["score"].to_numpy(float)
    grid = np.linspace(0, 1, 151)
    rng = np.random.default_rng(20260708)
    boot = []
    for _ in range(2000):
        idx = rng.choice(len(truth), len(truth), replace=True)
        if len(np.unique(truth[idx])) < 2:
            continue
        fpr, tpr, _ = roc_curve(truth[idx], heldout[idx])
        interpolated = np.interp(grid, fpr, tpr)
        interpolated[0], interpolated[-1] = 0, 1
        boot.append(interpolated)
    ci = np.vstack(boot)
    ci_frame = pd.DataFrame({"source_type": "composite_bootstrap_roc_95ci", "model_id": "CN_coupling", "model_label": "Composite coupling", "fpr": grid, "tpr_low": np.quantile(ci, .025, axis=0), "tpr_high": np.quantile(ci, .975, axis=0), "bootstrap_resamples": 2000, "bootstrap_seed": 20260708, "ci_method": "sample-level bootstrap of 36 held-out composite scores; interpolation grid 0..1 by 1/150"})
    OUT.mkdir(parents=True, exist_ok=True)
    score_audit = scores.copy(); score_audit["source_type"] = "held_out_fold_score"
    summary_audit = summary.copy(); summary_audit["source_type"] = "auc_summary"; summary_audit["ci_method"] = "Composite only: 2000 sample-level bootstrap resamples of 36 held-out scores, seed 20260708; no CI bands for single predictors."
    audit = pd.concat([score_audit, summary_audit, ci_frame], ignore_index=True, sort=False)
    audit.to_csv(OUT / "fig6g_comparative_foldinternal_roc_source_data.csv", index=False)
    pd.concat(directions, ignore_index=True).to_csv(OUT / "fig6g_comparative_single_predictor_fold_directions.csv", index=False)

    fig, ax = plt.subplots(figsize=(1.9845, 2.065), facecolor="white")
    ax.plot([0, 1], [0, 1], color="#999999", lw=0.8, ls="--", zorder=1)
    ax.fill_between(ci_frame.fpr, ci_frame.tpr_low, ci_frame.tpr_high, color=palette["Composite coupling"], alpha=.28, linewidth=0, zorder=2)
    for row in summary_rows:
        curve = pd.DataFrame([r for r in curve_rows if r["model_id"] == row["model_id"]])
        alpha = 1.0 if row["model_id"] == "CN_coupling" else .55
        ax.plot(curve.fpr, curve.tpr, color=palette[row["model_label"]], alpha=alpha, lw=1.6, label=f"{row['model_label']}\nAUC {row['auc']:.3f}", zorder=3)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
    ax.set_title("Held-out LOOCV ROC", pad=5)
    ax.text(-.18, 1.06, "g", transform=ax.transAxes, fontweight="bold", va="top")
    ax.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#DDDDDD", framealpha=.94, handlelength=1.2, borderpad=.35, labelspacing=.25)
    fig.subplots_adjust(left=.31, right=.97, top=.84, bottom=.24)
    apply_ariel_8(fig)
    fig.savefig(OUT / "fig6g_comparative_foldinternal_roc_composite_ci.pdf", format="pdf", facecolor="white")
    plt.close(fig)
    print(summary.to_csv(index=False))


if __name__ == "__main__":
    main()
