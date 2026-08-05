#!/usr/bin/env python3
"""Render the quantitative response panels (Figure 7b-d)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve


CR = "#D55E00"
NR = "#0072B2"
BLACK = "#333333"
GREY = "#8A8A8A"


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial"],
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def clean_axis(axis: plt.Axes) -> None:
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(length=2.5, width=0.7, pad=2)


def load_tables(evaluation: Path, score_id: str) -> tuple[pd.Series, pd.DataFrame, pd.Series]:
    metrics = pd.read_csv(evaluation / "metrics.tsv", sep="\t")
    classifications = pd.read_csv(evaluation / "patient_classifications.tsv", sep="\t")
    confusion = pd.read_csv(evaluation / "confusion_matrices.tsv", sep="\t")
    metric = metrics.loc[metrics.score_id.eq(score_id)]
    patients = classifications.loc[classifications.score_id.eq(score_id)].copy()
    matrix = confusion.loc[confusion.score_id.eq(score_id)]
    if len(metric) != 1 or len(matrix) != 1 or patients.patient_id.nunique() != len(patients):
        raise ValueError("response evaluation tables do not have one unambiguous selected score")
    return metric.iloc[0], patients, matrix.iloc[0]


def render(evaluation: Path, output: Path, score_id: str) -> dict[str, str]:
    set_style()
    metric, patients, matrix = load_tables(evaluation, score_id)
    patients = patients.sort_values(["score", "patient_id"], kind="mergesort").reset_index(drop=True)
    fig = plt.figure(figsize=(7.0, 2.25))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.45, 0.72], wspace=0.70)

    roc_axis = fig.add_subplot(grid[0, 0])
    fpr, tpr, _ = roc_curve(patients.response.to_numpy(int), patients.score.to_numpy(float))
    roc_axis.plot(fpr, tpr, color=BLACK, lw=1.6)
    roc_axis.plot([0, 1], [0, 1], color="#C7C7C7", lw=0.8, ls="--")
    roc_axis.set(xlim=(-0.02, 1.02), ylim=(-0.02, 1.02), xlabel="1 − specificity", ylabel="Sensitivity")
    roc_axis.set_xticks([0, 0.5, 1])
    roc_axis.set_yticks([0, 0.5, 1])
    roc_axis.set_title("Locked temporal validation", pad=5)
    roc_axis.text(
        0.97,
        0.08,
        f"AUC {float(metric.roc_auc):.3f}\n"
        f"95% CI {float(metric.roc_auc_ci025):.2f}–{float(metric.roc_auc_ci975):.2f}\n"
        f"exact $p$ = {float(metric.exact_p_auc_one_sided):.3f}",
        ha="right",
        va="bottom",
        fontsize=7,
    )
    clean_axis(roc_axis)
    roc_axis.text(-0.43, 1.12, "b", transform=roc_axis.transAxes, fontsize=11, fontweight="bold")

    score_axis = fig.add_subplot(grid[0, 1])
    x = np.arange(len(patients))
    for response, color, marker, label in ((1, CR, "o", "CR"), (0, NR, "s", "NR")):
        keep = patients.response.eq(response).to_numpy()
        score_axis.scatter(
            x[keep], patients.score.to_numpy(float)[keep], s=21, color=color,
            marker=marker, edgecolor="white", linewidth=0.45, label=label, zorder=3,
        )
    threshold = float(patients.threshold.iloc[0])
    score_axis.axhline(threshold, color=BLACK, lw=0.8, ls=(0, (3, 2)))
    score_axis.text(len(patients) - 0.2, threshold - 0.012, "development threshold", ha="right", color=GREY, fontsize=6.5)
    score_axis.set_xticks(x, patients.patient_id.astype(str), rotation=90)
    score_axis.set_ylabel("M0 response score")
    score_axis.set_title("One score per patient", pad=5)
    score_axis.legend(frameon=False, ncol=2, loc="upper center", handletextpad=0.25, columnspacing=0.8)
    padding = max(0.02, 0.12 * float(patients.score.max() - patients.score.min()))
    score_axis.set_ylim(float(patients.score.min() - padding), float(patients.score.max() + padding))
    clean_axis(score_axis)
    score_axis.text(-0.35, 1.12, "c", transform=score_axis.transAxes, fontsize=11, fontweight="bold")

    confusion_axis = fig.add_subplot(grid[0, 2])
    values = np.array([[int(matrix.tn), int(matrix.fp)], [int(matrix.fn), int(matrix.tp)]])
    confusion_axis.imshow(values, cmap=mpl.colors.ListedColormap(["#F3F3F3", "#B8D8E8"]), vmin=0, vmax=max(1, values.max()))
    for (row, column), value in np.ndenumerate(values):
        confusion_axis.text(column, row, str(value), ha="center", va="center", fontsize=9)
    confusion_axis.set_xticks([0, 1], ["NR", "CR"])
    confusion_axis.set_yticks([0, 1], ["NR", "CR"])
    confusion_axis.set_xlabel("Prediction")
    confusion_axis.set_ylabel("True")
    confusion_axis.set_title("Confusion matrix", pad=5)
    confusion_axis.tick_params(length=0)
    confusion_axis.text(-0.60, 1.12, "d", transform=confusion_axis.transAxes, fontsize=11, fontweight="bold")

    output.mkdir(parents=True, exist_ok=True)
    stem = output / "Figure7_bcd_response"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    result = {suffix: str(stem.with_suffix(suffix)) for suffix in (".pdf", ".svg", ".png")}
    (output / "Figure7_bcd_response.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--score-id", default="prov_score")
    args = parser.parse_args()
    render(args.evaluation, args.output, args.score_id)


if __name__ == "__main__":
    main()
