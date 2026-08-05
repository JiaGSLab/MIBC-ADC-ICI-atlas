#!/usr/bin/env python3
"""Refit regional heads on measured patients and infer every temporal tile."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

try:
    from .core import MARKERS, feature_columns, load_patient, spatially_balanced_indices
except ImportError:  # direct script execution
    from core import MARKERS, feature_columns, load_patient, spatially_balanced_indices


IDENTITY = ["patient_id", "sample_id", "tile_id", "x_min", "y_min", "x_max", "y_max"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    output = Path(config["output_root"]) / ("full_refit_smoke" if args.smoke else "full_refit")
    map_root = output / "temporal_maps"
    map_root.mkdir(parents=True, exist_ok=True)

    development = pd.read_csv(
        config["development_manifest"], sep="\t", dtype=str,
        usecols=lambda column: column in {"patient_id", "sample_id", "mpp", "prov_path"},
    )
    temporal = pd.read_csv(
        config["temporal_manifest"], sep="\t", dtype=str,
        usecols=lambda column: column in {"patient_id", "sample_id", "mpp", "prov_path"},
    )
    folds = pd.read_csv(config["development_folds"], sep="\t", dtype={"patient_id": str})
    targets = pd.read_parquet(config["regional_targets"])
    measured = set(folds.patient_id.astype(str))
    development = development.loc[development.patient_id.astype(str).isin(measured)].copy()
    if set(development.patient_id.astype(str)) != measured:
        raise ValueError("fold manifest and measured development patients differ")
    if args.smoke:
        development = development.sort_values("patient_id").head(5)
        temporal = temporal.sort_values("patient_id").head(1)
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

    cap = min(int(config["primary_training_cap_per_patient"]), 128) if args.smoke else int(
        config["primary_training_cap_per_patient"]
    )
    indices = {
        pid: spatially_balanced_indices(patient, cap, int(config["model_seed"]))
        for pid, patient in patients.items()
    }
    sampling_rows = []
    for pid, patient in patients.items():
        selected = indices[pid]
        all_bins = set(map(tuple, patient.bins))
        selected_bins = set(map(tuple, patient.bins[selected]))
        sampling_rows.append(
            {
                "patient_id": pid,
                "eligible_valid_tiles": patient.n_tiles,
                "selected_tiles": int(len(selected)),
                "total_spatial_bins": len(all_bins),
                "selected_spatial_bins": len(selected_bins),
                "all_bins_represented": selected_bins == all_bins,
            }
        )
    sampling = pd.DataFrame(sampling_rows)
    if not sampling.all_bins_represented.all():
        raise ValueError("full-refit sampler missed an occupied spatial bin")

    train_x = np.concatenate([patients[pid].features[indices[pid]] for pid in sorted(patients)])
    train_y = np.concatenate([patients[pid].targets[indices[pid]] for pid in sorted(patients)])
    scaler = StandardScaler(copy=True).fit(train_x)
    scaled = scaler.transform(train_x)
    pca = PCA(
        n_components=int(config["pca_dim"]),
        svd_solver="randomized",
        random_state=int(config["model_seed"]),
    ).fit(scaled)
    z = pca.transform(scaled)
    heads = {}
    for marker_index, marker in enumerate(MARKERS):
        values = train_y[:, marker_index].astype(np.float64)
        finite = np.isfinite(values)
        heads[marker] = Ridge(alpha=float(config["ridge_alpha"])).fit(z[finite], values[finite])

    bundle = {
        "scaler": scaler,
        "pca": pca,
        "heads": heads,
        "markers": MARKERS,
        "train_cap_per_patient": cap,
        "sampling": "deterministic_physical-bin_round-robin",
        "development_patients": sorted(patients),
        "temporal_patients_used_for_fit_or_selection": [],
        "outcome_access": False,
    }
    bundle_path = output / "regional_full_refit.joblib"
    joblib.dump(bundle, bundle_path, compress=3)
    sampling_path = output / "training_sampling_ledger.tsv"
    sampling.to_csv(sampling_path, sep="\t", index=False)

    ledger_rows = []
    for row in temporal.sort_values("patient_id").itertuples(index=False):
        source = Path(str(row.prov_path))
        frame = pd.read_parquet(source)
        frame["tile_id"] = frame.tile_id.astype(str)
        if frame.tile_id.duplicated().any():
            raise ValueError(f"duplicate temporal tile: {row.sample_id}")
        columns = feature_columns(frame.columns, int(config["expected_feature_dim"]))
        matrix = frame[columns].to_numpy(dtype=np.float32, copy=True)
        if not np.isfinite(matrix).all():
            raise ValueError(f"non-finite temporal embedding: {row.sample_id}")
        latent = pca.transform(scaler.transform(matrix))
        output_frame = frame[IDENTITY].copy()
        for marker in MARKERS:
            output_frame[f"{marker}_predicted"] = heads[marker].predict(latent).astype(np.float32)
        path = map_root / f"{row.sample_id}.parquet"
        output_frame.to_parquet(path, index=False)
        if len(output_frame) != len(frame):
            raise RuntimeError(f"temporal inference was capped: {row.sample_id}")
        ledger_rows.append(
            {
                "patient_id": str(row.patient_id),
                "sample_id": str(row.sample_id),
                "input_sha256": sha256_file(source),
                "all_foreground_tiles_inferred": True,
                "tile_count": int(len(output_frame)),
                "output_file": str(path.name),
                "output_sha256": sha256_file(path),
            }
        )
    ledger = pd.DataFrame(ledger_rows)
    ledger_path = output / "temporal_map_ledger.tsv"
    ledger.to_csv(ledger_path, sep="\t", index=False)
    status = {
        "status": "PASS",
        "development_patients": int(len(development)),
        "temporal_patients": int(len(temporal)),
        "training_valid_rows_selected": int(sampling.selected_tiles.sum()),
        "temporal_rows_inferred": int(ledger.tile_count.sum()),
        "temporal_inference_uncapped": True,
        "outcome_access": False,
        "smoke": bool(args.smoke),
        "outputs": {
            "bundle": sha256_file(bundle_path),
            "sampling": sha256_file(sampling_path),
            "temporal_ledger": sha256_file(ledger_path),
        },
    }
    (output / "complete.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

