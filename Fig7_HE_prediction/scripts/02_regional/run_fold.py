#!/usr/bin/env python3
"""Run one leakage-controlled fold of the full-coverage regional analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

try:
    from .core import MARKERS, load_patient, permute_within_patient, safe_spearman, spatially_balanced_indices
except ImportError:  # direct script execution
    from core import MARKERS, load_patient, permute_within_patient, safe_spearman, spatially_balanced_indices


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_config(path: Path) -> dict:
    config = json.loads(path.read_text(encoding="utf-8"))
    required = {
        "development_manifest", "development_folds", "regional_targets", "output_root",
        "expected_feature_dim", "primary_training_cap_per_patient",
        "sensitivity_training_cap_per_patient", "spatial_bin_um", "pca_dim",
        "ridge_alpha", "model_seed", "control_seeds", "target_columns",
    }
    missing = required - set(config)
    if missing:
        raise ValueError(f"regional config keys missing: {sorted(missing)}")
    return config


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--fold", type=int, required=True, choices=range(5))
    parser.add_argument("--regime", choices=("primary", "sensitivity_2048"), required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    start = time.time()
    config = load_config(args.config)
    cap = int(
        config["primary_training_cap_per_patient"]
        if args.regime == "primary"
        else config["sensitivity_training_cap_per_patient"]
    )
    output = Path(config["output_root"]) / args.regime / f"fold_{args.fold:02d}"
    output.mkdir(parents=True, exist_ok=True)

    development = pd.read_csv(
        config["development_manifest"], sep="\t", dtype=str,
        usecols=lambda column: column in {"patient_id", "sample_id", "mpp", "prov_path"},
    )
    folds = pd.read_csv(config["development_folds"], sep="\t", dtype={"patient_id": str})
    targets = pd.read_parquet(config["regional_targets"])
    measured = set(folds.patient_id.astype(str))
    development = development.loc[development.patient_id.astype(str).isin(measured)].copy()
    if set(development.patient_id.astype(str)) != measured:
        raise ValueError("fold manifest and development manifest patient sets differ")
    fold_map = folds.set_index(folds.patient_id.astype(str)).outer_fold.astype(int).to_dict()

    patients = {}
    for row in development.sort_values("patient_id").itertuples(index=False):
        patient = load_patient(
            row,
            fold_map[str(row.patient_id)],
            targets,
            target_columns=config["target_columns"],
            expected_dim=int(config["expected_feature_dim"]),
            bin_um=float(config["spatial_bin_um"]),
        )
        patients[patient.patient_id] = patient

    train_ids = sorted(pid for pid, patient in patients.items() if patient.outer_fold != args.fold)
    held_ids = sorted(pid for pid, patient in patients.items() if patient.outer_fold == args.fold)
    if args.smoke:
        train_ids = train_ids[:5]
        held_ids = held_ids[:1]
        cap = min(cap, 128)
    if not train_ids or not held_ids or set(train_ids) & set(held_ids):
        raise ValueError("invalid train/held-out patient split")

    indices = {
        pid: spatially_balanced_indices(patients[pid], cap, int(config["model_seed"]))
        for pid in train_ids
    }
    sampling_rows = []
    for pid in train_ids:
        patient = patients[pid]
        selected = indices[pid]
        all_bins = set(map(tuple, patient.bins))
        selected_bins = set(map(tuple, patient.bins[selected]))
        sampling_rows.append(
            {
                "patient_id": pid,
                "eligible_tiles": patient.n_tiles,
                "selected_tiles": int(len(selected)),
                "total_spatial_bins": len(all_bins),
                "selected_spatial_bins": len(selected_bins),
                "all_bins_represented": selected_bins == all_bins,
                "eligible_x_min": int(patient.coords[:, 0].min()),
                "eligible_x_max": int(patient.coords[:, 2].max()),
                "eligible_y_min": int(patient.coords[:, 1].min()),
                "eligible_y_max": int(patient.coords[:, 3].max()),
                "selected_x_min": int(patient.coords[selected, 0].min()),
                "selected_x_max": int(patient.coords[selected, 2].max()),
                "selected_y_min": int(patient.coords[selected, 1].min()),
                "selected_y_max": int(patient.coords[selected, 3].max()),
            }
        )
    sampling = pd.DataFrame(sampling_rows)
    if not sampling.all_bins_represented.all():
        raise ValueError("spatially balanced sampler missed an occupied physical bin")

    train_x = np.concatenate([patients[pid].features[indices[pid]] for pid in train_ids])
    train_y = np.concatenate([patients[pid].targets[indices[pid]] for pid in train_ids])
    train_patient = np.concatenate([np.repeat(pid, len(indices[pid])) for pid in train_ids])
    held_x = np.concatenate([patients[pid].features for pid in held_ids])
    held_y = np.concatenate([patients[pid].targets for pid in held_ids])
    held_offsets = np.cumsum([0] + [patients[pid].n_tiles for pid in held_ids])

    scaler = StandardScaler(copy=True).fit(train_x)
    scaled = scaler.transform(train_x)
    pca = PCA(
        n_components=int(config["pca_dim"]),
        svd_solver="randomized",
        random_state=int(config["model_seed"]),
    ).fit(scaled)
    z_train = pca.transform(scaled)
    z_held = pca.transform(scaler.transform(held_x))

    prediction_rows = []
    control_rows = []
    metric_rows = []
    heads = {}
    for marker_index, marker in enumerate(MARKERS):
        values = train_y[:, marker_index].astype(np.float64)
        finite = np.isfinite(values)
        head = Ridge(alpha=float(config["ridge_alpha"])).fit(z_train[finite], values[finite])
        heads[marker] = head
        predictions = [("none", head.predict(z_held))]
        for seed in map(int, config["control_seeds"]):
            shuffled = permute_within_patient(values, train_patient, marker=marker, seed=seed)
            valid = np.isfinite(shuffled)
            control_head = Ridge(alpha=float(config["ridge_alpha"])).fit(z_train[valid], shuffled[valid])
            predictions.append((f"within_patient_shuffle_{seed}", control_head.predict(z_held)))
        for control, prediction in predictions:
            for index, patient_id in enumerate(held_ids):
                patient = patients[patient_id]
                lo, hi = int(held_offsets[index]), int(held_offsets[index + 1])
                actual = held_y[lo:hi, marker_index].astype(float)
                predicted = prediction[lo:hi].astype(float)
                metric_rows.append(
                    {
                        "regime": args.regime,
                        "outer_fold": args.fold,
                        "patient_id": patient_id,
                        "sample_id": patient.sample_id,
                        "marker": marker,
                        "control": control,
                        "n_tiles": int(np.isfinite(actual).sum()),
                        "spearman": safe_spearman(actual, predicted),
                    }
                )
                table = pd.DataFrame(
                    {
                        "regime": args.regime,
                        "outer_fold": args.fold,
                        "patient_id": patient_id,
                        "sample_id": patient.sample_id,
                        "tile_id": patient.tile_id,
                        "x_min": patient.coords[:, 0],
                        "y_min": patient.coords[:, 1],
                        "x_max": patient.coords[:, 2],
                        "y_max": patient.coords[:, 3],
                        "marker": marker,
                        "control": control,
                        "actual": actual,
                        "prediction": predicted,
                    }
                )
                finite_held = np.isfinite(table.actual) & np.isfinite(table.prediction)
                (prediction_rows if control == "none" else control_rows).append(table.loc[finite_held].copy())

    outputs = {
        "oof_predictions.parquet": pd.concat(prediction_rows, ignore_index=True),
        "control_oof_predictions.parquet": pd.concat(control_rows, ignore_index=True),
        "patient_metrics.tsv": pd.DataFrame(metric_rows),
        "training_sampling_ledger.tsv": sampling,
    }
    for name, frame in outputs.items():
        path = output / name
        if path.suffix == ".parquet":
            frame.to_parquet(path, index=False)
        else:
            frame.to_csv(path, sep="\t", index=False)
    bundle_path = output / "fold_bundle.joblib"
    joblib.dump(
        {
            "scaler": scaler,
            "pca": pca,
            "heads": heads,
            "markers": MARKERS,
            "train_cap": cap,
            "bin_um": float(config["spatial_bin_um"]),
            "train_patients": train_ids,
            "heldout_patients": held_ids,
            "outcome_access": False,
        },
        bundle_path,
        compress=3,
    )
    metadata = {
        "status": "PASS",
        "regime": args.regime,
        "outer_fold": args.fold,
        "train_cap": cap,
        "train_patients": train_ids,
        "heldout_patients": held_ids,
        "training_rows": int(len(z_train)),
        "heldout_rows_uncapped": int(len(z_held)),
        "outcome_access": False,
        "smoke": bool(args.smoke),
        "elapsed_seconds": float(time.time() - start),
        "runtime": {"python": sys.version},
        "outputs": {
            **{name: sha256_file(output / name) for name in outputs},
            "fold_bundle.joblib": sha256_file(bundle_path),
        },
    }
    (output / "complete.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

