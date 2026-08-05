#!/usr/bin/env python3
"""Outcome evaluation for frozen H&E response scores.

The validation scores must be produced before this script is given access to
validation outcomes. Thresholds are derived from development OOF predictions
only. AUC/AP confidence intervals use class-stratified bootstrap resampling,
and one-sided P values enumerate all label allocations with the observed class
balance (1,287 allocations for 8 CR and 5 NR patients).
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_SEED = 20260727
DEFAULT_BOOTSTRAPS = 10_000


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def auc(y: np.ndarray, score: np.ndarray) -> float:
    labels = np.asarray(y, dtype=int)
    values = np.asarray(score, dtype=float)
    positive = values[labels == 1][:, None]
    negative = values[labels == 0][None, :]
    if len(positive) == 0 or negative.shape[1] == 0 or not np.isfinite(values).all():
        raise ValueError("invalid binary metric input")
    wins = (positive > negative).sum() + 0.5 * (positive == negative).sum()
    return float(wins / (len(positive) * negative.shape[1]))


def average_precision(y: np.ndarray, score: np.ndarray) -> float:
    labels = np.asarray(y, dtype=int)
    values = np.asarray(score, dtype=float)
    total_positive = int(labels.sum())
    if total_positive == 0 or not np.isfinite(values).all():
        raise ValueError("invalid average-precision input")
    true_positive = retrieved = 0
    previous_recall = result = 0.0
    for threshold in sorted(set(values.tolist()), reverse=True):
        mask = values == threshold
        true_positive += int(labels[mask].sum())
        retrieved += int(mask.sum())
        recall = true_positive / total_positive
        precision = true_positive / retrieved
        result += (recall - previous_recall) * precision
        previous_recall = recall
    return float(result)


def confusion(y: np.ndarray, score: np.ndarray, threshold: float) -> dict[str, float | int]:
    labels = np.asarray(y, dtype=int)
    predicted = np.asarray(score, dtype=float) >= float(threshold)
    tp = int((predicted & (labels == 1)).sum())
    fn = int((~predicted & (labels == 1)).sum())
    tn = int((~predicted & (labels == 0)).sum())
    fp = int((predicted & (labels == 0)).sum())
    return {
        "tp": tp,
        "fn": fn,
        "tn": tn,
        "fp": fp,
        "sensitivity": tp / (tp + fn),
        "specificity": tn / (tn + fp),
        "accuracy": (tp + tn) / len(labels),
        "ppv": tp / (tp + fp) if tp + fp else math.nan,
        "npv": tn / (tn + fn) if tn + fn else math.nan,
    }


def development_youden(y: np.ndarray, score: np.ndarray) -> dict[str, float | int]:
    rows = []
    for threshold in sorted(set(np.asarray(score, dtype=float).tolist())):
        row = confusion(y, score, threshold)
        row["threshold"] = threshold
        row["youden_j"] = row["sensitivity"] + row["specificity"] - 1
        rows.append(row)
    return max(
        rows,
        key=lambda row: (
            round(float(row["youden_j"]), 15),
            round(float(row["specificity"]), 15),
            -float(row["threshold"]),
        ),
    )


def stratified_bootstrap(
    y: np.ndarray,
    score: np.ndarray,
    *,
    iterations: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    labels = np.asarray(y, dtype=int)
    values = np.asarray(score, dtype=float)
    positive = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    rng = np.random.default_rng(seed)
    auc_values = np.empty(iterations, dtype=float)
    ap_values = np.empty(iterations, dtype=float)
    for index in range(iterations):
        sampled = np.r_[
            rng.choice(positive, len(positive), replace=True),
            rng.choice(negative, len(negative), replace=True),
        ]
        auc_values[index] = auc(labels[sampled], values[sampled])
        ap_values[index] = average_precision(labels[sampled], values[sampled])
    return auc_values, ap_values


def exact_one_sided_p(y: np.ndarray, score: np.ndarray) -> tuple[float, float, int]:
    labels = np.asarray(y, dtype=int)
    values = np.asarray(score, dtype=float)
    observed_auc = auc(labels, values)
    observed_ap = average_precision(labels, values)
    allocation_auc = []
    allocation_ap = []
    for positive_indices in itertools.combinations(range(len(labels)), int(labels.sum())):
        allocation = np.zeros(len(labels), dtype=int)
        allocation[list(positive_indices)] = 1
        allocation_auc.append(auc(allocation, values))
        allocation_ap.append(average_precision(allocation, values))
    auc_values = np.asarray(allocation_auc)
    ap_values = np.asarray(allocation_ap)
    return (
        float(np.mean(auc_values >= observed_auc - 1e-15)),
        float(np.mean(ap_values >= observed_ap - 1e-15)),
        int(len(auc_values)),
    )


def load_inputs(
    scores_path: Path,
    outcomes_path: Path,
    development_path: Path,
    score_columns: list[str],
    expected_validation: int,
    expected_development: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.read_csv(scores_path, sep="\t")
    outcomes = pd.read_csv(outcomes_path, sep="\t")
    development = pd.read_csv(development_path, sep="\t")
    for name, frame in (("scores", scores), ("outcomes", outcomes), ("development", development)):
        if "patient_id" not in frame or frame.patient_id.astype(str).duplicated().any():
            raise ValueError(f"{name} must contain one row per patient_id")
    if "response" not in outcomes or "response" not in development:
        raise ValueError("binary response is required in outcome and development tables")
    missing = set(score_columns) - set(scores) | set(score_columns) - set(development)
    if missing:
        raise ValueError(f"score columns missing: {sorted(missing)}")
    joined = outcomes[["patient_id", "response"]].merge(
        scores[["patient_id", *(["sample_id"] if "sample_id" in scores else []), *score_columns]],
        on="patient_id",
        how="inner",
        validate="one_to_one",
    )
    if len(joined) != expected_validation or len(development) != expected_development:
        raise ValueError("cohort size differs from the frozen analysis contract")
    for frame in (joined, development):
        if set(frame.response.astype(int).unique()) != {0, 1}:
            raise ValueError("response must be binary")
        if not np.isfinite(frame[score_columns].to_numpy(float)).all():
            raise ValueError("non-finite score")
    return joined.sort_values("patient_id").reset_index(drop=True), development.sort_values("patient_id").reset_index(drop=True)


def evaluate(
    joined: pd.DataFrame,
    development: pd.DataFrame,
    score_columns: list[str],
    *,
    bootstraps: int,
    seed: int,
) -> dict[str, pd.DataFrame]:
    y = joined.response.to_numpy(int)
    development_y = development.response.to_numpy(int)
    metric_rows = []
    threshold_rows = []
    confusion_rows = []
    classification_rows = []
    influence_rows = []
    for score_id in score_columns:
        score = joined[score_id].to_numpy(float)
        development_score = development[score_id].to_numpy(float)
        observed_auc = auc(y, score)
        observed_ap = average_precision(y, score)
        boot_auc, boot_ap = stratified_bootstrap(
            y,
            score,
            iterations=bootstraps,
            seed=seed,
        )
        p_auc, p_ap, allocation_count = exact_one_sided_p(y, score)
        metric_rows.append(
            {
                "score_id": score_id,
                "n": len(y),
                "n_positive": int(y.sum()),
                "n_negative": int((y == 0).sum()),
                "development_oof_auc": auc(development_y, development_score),
                "roc_auc": observed_auc,
                "roc_auc_ci025": float(np.quantile(boot_auc, 0.025)),
                "roc_auc_ci975": float(np.quantile(boot_auc, 0.975)),
                "average_precision": observed_ap,
                "average_precision_ci025": float(np.quantile(boot_ap, 0.025)),
                "average_precision_ci975": float(np.quantile(boot_ap, 0.975)),
                "exact_p_auc_one_sided": p_auc,
                "exact_p_average_precision_one_sided": p_ap,
                "n_exact_allocations": allocation_count,
                "n_bootstrap": bootstraps,
                "bootstrap_seed": seed,
            }
        )
        selected = development_youden(development_y, development_score)
        threshold = float(selected["threshold"])
        threshold_rows.append(
            {
                "score_id": score_id,
                "threshold_source": "development_OOF_only",
                "threshold_rule": "max_Youden_J; tie=higher_specificity_then_lowest_threshold",
                "threshold": threshold,
                "development_youden_j": selected["youden_j"],
                "development_sensitivity": selected["sensitivity"],
                "development_specificity": selected["specificity"],
            }
        )
        observed_confusion = confusion(y, score, threshold)
        confusion_rows.append({"score_id": score_id, "threshold": threshold, **observed_confusion})
        for patient_id, response, value in zip(joined.patient_id.astype(str), y, score):
            classification_rows.append(
                {
                    "score_id": score_id,
                    "patient_id": patient_id,
                    "response": int(response),
                    "score": float(value),
                    "threshold": threshold,
                    "predicted_response": int(value >= threshold),
                }
            )
        for index, patient_id in enumerate(joined.patient_id.astype(str)):
            keep = np.arange(len(y)) != index
            influence_rows.append(
                {
                    "score_id": score_id,
                    "deleted_patient_id": patient_id,
                    "full_auc": observed_auc,
                    "leave_one_out_auc": auc(y[keep], score[keep]),
                    "full_average_precision": observed_ap,
                    "leave_one_out_average_precision": average_precision(y[keep], score[keep]),
                }
            )
    metrics = pd.DataFrame(metric_rows)
    metrics["descriptive_auc_rank"] = metrics.roc_auc.rank(method="min", ascending=False).astype(int)
    metrics["descriptive_best_auc_tie"] = metrics.roc_auc.eq(metrics.roc_auc.max())
    return {
        "metrics": metrics,
        "thresholds": pd.DataFrame(threshold_rows),
        "confusion": pd.DataFrame(confusion_rows),
        "classifications": pd.DataFrame(classification_rows),
        "influence": pd.DataFrame(influence_rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate frozen response scores once outcomes are unlocked.")
    parser.add_argument("--validation-scores", type=Path, required=True)
    parser.add_argument("--validation-outcomes", type=Path, required=True)
    parser.add_argument("--development-oof", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--score-columns", default="prov_score")
    parser.add_argument("--expected-validation", type=int, default=13)
    parser.add_argument("--expected-development", type=int, default=35)
    parser.add_argument("--bootstraps", type=int, default=DEFAULT_BOOTSTRAPS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    score_columns = [column.strip() for column in args.score_columns.split(",") if column.strip()]
    joined, development = load_inputs(
        args.validation_scores,
        args.validation_outcomes,
        args.development_oof,
        score_columns,
        args.expected_validation,
        args.expected_development,
    )
    tables = evaluate(joined, development, score_columns, bootstraps=args.bootstraps, seed=args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    names = {
        "metrics": "metrics.tsv",
        "thresholds": "threshold_ledger.tsv",
        "confusion": "confusion_matrices.tsv",
        "classifications": "patient_classifications.tsv",
        "influence": "leave_one_out_influence.tsv",
    }
    for key, filename in names.items():
        tables[key].to_csv(args.output_dir / filename, sep="\t", index=False, float_format="%.17g")
    metadata = {
        "status": "COMPLETE",
        "validation_outcomes_used_only_in_this_evaluation_step": True,
        "threshold_source": "development OOF only",
        "score_columns": score_columns,
        "seed": args.seed,
        "bootstraps": args.bootstraps,
        "input_sha256": {
            "validation_scores": sha256_file(args.validation_scores),
            "validation_outcomes": sha256_file(args.validation_outcomes),
            "development_oof": sha256_file(args.development_oof),
        },
    }
    (args.output_dir / "evaluation_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
