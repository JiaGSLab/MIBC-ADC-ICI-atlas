#!/usr/bin/env python3
"""Audit and aggregate full-coverage regional out-of-fold predictions."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from .core import combined_epithelial, safe_spearman
except ImportError:  # direct script execution
    from core import combined_epithelial, safe_spearman


BASE_MARKERS = ("Vascular", "PanCK", "HER2")
COMBINED = "Epithelial (PanCK + HER2)"
MARKERS = (*BASE_MARKERS, COMBINED)
REGIMES = ("primary", "sensitivity_2048")
KEY = [
    "regime", "outer_fold", "patient_id", "sample_id", "tile_id",
    "x_min", "y_min", "x_max", "y_max", "control",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def add_combined(frame: pd.DataFrame) -> pd.DataFrame:
    subset = frame.loc[frame.marker.isin(["PanCK", "HER2"])].copy()
    actual = subset.pivot(index=KEY, columns="marker", values="actual").reindex(columns=["PanCK", "HER2"])
    predicted = subset.pivot(index=KEY, columns="marker", values="prediction").reindex(columns=["PanCK", "HER2"])
    if actual.isna().any().any() or predicted.isna().any().any():
        raise ValueError("PanCK-signal channel pivot is incomplete")
    output = actual.reset_index()[KEY].copy()
    output["marker"] = COMBINED
    output["actual"] = combined_epithelial(actual.PanCK.to_numpy(float), actual.HER2.to_numpy(float))
    output["prediction"] = combined_epithelial(predicted.PanCK.to_numpy(float), predicted.HER2.to_numpy(float))
    return output


def metrics_from_predictions(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, group in frame.groupby(
        ["regime", "outer_fold", "patient_id", "sample_id", "marker", "control"], sort=True
    ):
        regime, fold, patient, sample, marker, control = keys
        rows.append(
            {
                "regime": regime,
                "outer_fold": int(fold),
                "patient_id": str(patient),
                "sample_id": str(sample),
                "marker": marker,
                "control": control,
                "n_tiles": int(len(group)),
                "spearman": safe_spearman(group.actual.to_numpy(float), group.prediction.to_numpy(float)),
            }
        )
    return pd.DataFrame(rows)


def bootstrap_ci(values: np.ndarray, *, seed: int, n_resamples: int) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    rng = np.random.default_rng(int(seed))
    samples = rng.choice(values, size=(int(n_resamples), len(values)), replace=True)
    medians = np.median(samples, axis=1)
    return float(np.quantile(medians, 0.025)), float(np.quantile(medians, 0.975))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    root = Path(config["output_root"])
    real_parts, control_parts, sampling_parts = [], [], []
    metadata = []
    for regime in REGIMES:
        for fold in range(5):
            folder = root / regime / f"fold_{fold:02d}"
            status = json.loads((folder / "complete.json").read_text(encoding="utf-8"))
            if status.get("status") != "PASS" or int(status.get("outer_fold")) != fold:
                raise ValueError(f"fold status failed: {folder}")
            if set(status["train_patients"]) & set(status["heldout_patients"]):
                raise ValueError(f"patient leakage recorded: {folder}")
            real_parts.append(pd.read_parquet(folder / "oof_predictions.parquet"))
            control_parts.append(pd.read_parquet(folder / "control_oof_predictions.parquet"))
            sampling = pd.read_csv(folder / "training_sampling_ledger.tsv", sep="\t")
            if not sampling.all_bins_represented.astype(bool).all():
                raise ValueError(f"sampling coverage failed: {folder}")
            sampling.insert(0, "regime", regime)
            sampling.insert(1, "outer_fold", fold)
            sampling_parts.append(sampling)
            metadata.append(status)
    real = pd.concat(real_parts, ignore_index=True)
    controls = pd.concat(control_parts, ignore_index=True)
    all_base = pd.concat([real, controls], ignore_index=True)
    if all_base.duplicated(["regime", "patient_id", "marker", "control", "tile_id"]).any():
        raise ValueError("duplicate OOF tile prediction")
    combined = add_combined(all_base)
    all_predictions = pd.concat([all_base, combined], ignore_index=True)
    metrics = metrics_from_predictions(all_predictions)
    expected_controls = {"none", *[f"within_patient_shuffle_{int(seed)}" for seed in config["control_seeds"]]}
    if set(metrics.control.astype(str)) != expected_controls:
        raise ValueError("control set mismatch")

    summary_rows = []
    for regime in REGIMES:
        for marker in MARKERS:
            subset = metrics.loc[metrics.regime.eq(regime) & metrics.marker.eq(marker)]
            observed = subset.loc[subset.control.eq("none"), "spearman"].to_numpy(float)
            controls_by_seed = subset.loc[~subset.control.eq("none")].groupby("control").spearman.median()
            lower, upper = bootstrap_ci(
                observed,
                seed=int(config["bootstrap_seed"]),
                n_resamples=int(config["bootstrap_resamples"]),
            )
            median = float(np.nanmedian(observed))
            maximum_control = float(controls_by_seed.max())
            summary_rows.append(
                {
                    "regime": regime,
                    "marker": marker,
                    "n_patients": int(len(observed)),
                    "median_patient_spearman": median,
                    "bootstrap_95ci_lower": lower,
                    "bootstrap_95ci_upper": upper,
                    "max_shuffled_median_spearman": maximum_control,
                    "passes_descriptive_gate": bool(median >= 0.30 and median > maximum_control),
                    "positive_patient_fraction": float(np.mean(observed > 0)),
                    "median_tiles_per_patient": float(
                        subset.loc[subset.control.eq("none"), "n_tiles"].median()
                    ),
                }
            )
    summary = pd.DataFrame(summary_rows)
    sensitivity = metrics.loc[metrics.control.eq("none")].pivot_table(
        index=["patient_id", "marker"], columns="regime", values="spearman"
    ).reset_index()
    sensitivity["primary_minus_sensitivity"] = sensitivity.primary - sensitivity.sensitivity_2048

    root.mkdir(parents=True, exist_ok=True)
    paths = {
        "oof_predictions.parquet": real,
        "control_oof_predictions.parquet": controls,
        "combined_epithelial_oof.parquet": combined,
        "patient_metrics.tsv": metrics,
        "regional_summary.tsv": summary,
        "sensitivity.tsv": sensitivity,
        "training_sampling_ledger.tsv": pd.concat(sampling_parts, ignore_index=True),
    }
    for name, frame in paths.items():
        path = root / name
        if path.suffix == ".parquet":
            frame.to_parquet(path, index=False)
        else:
            frame.to_csv(path, sep="\t", index=False)
    panel_e = summary.loc[
        summary.regime.eq("primary") & summary.marker.isin(["Vascular", COMBINED, "HER2"])
    ].copy()
    panel_e["phenotype"] = panel_e.marker.map(
        {"Vascular": "Vascular", COMBINED: "PanCK signal", "HER2": "HER2"}
    )
    panel_e["analysis_status"] = np.where(panel_e.marker.eq(COMBINED), "post_hoc_combined", "prespecified_component")
    panel_e.to_csv(root / "panel_e_regional.tsv", sep="\t", index=False)
    status = {
        "status": "PASS",
        "evaluation_uncapped": True,
        "outcome_access": False,
        "regimes": list(REGIMES),
        "markers": list(MARKERS),
        "outputs": {name: sha256_file(root / name) for name in paths},
        "panel_e_regional.tsv": sha256_file(root / "panel_e_regional.tsv"),
    }
    (root / "aggregate_complete.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
