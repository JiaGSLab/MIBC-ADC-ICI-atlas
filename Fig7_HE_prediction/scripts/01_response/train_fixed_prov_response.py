#!/usr/bin/env python3
"""Fit the frozen Prov-GigaPath response model used in Figure 7.

The model specification was fixed after development-cohort screening:
patient-level mean pooling, StandardScaler, two-component PCA and an L2
logistic regression with C=0.01.  Validation outcomes are deliberately not
accepted by this script; they are joined only by ``evaluate_external_cohort``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


SEED = 20260725
EXPECTED_DIM = 1536
PCA_DIM = 2
LOGISTIC_C = 0.01
OUTCOME_LIKE = re.compile(r"(^|_)(response|outcome|label|target)(_|$)", re.I)


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_config(path: Path) -> dict:
    config = json.loads(path.read_text(encoding="utf-8"))
    required = {
        "development_manifest",
        "development_folds",
        "validation_manifest",
        "output_root",
        "expected_development_patients",
        "expected_validation_patients",
    }
    missing = required - set(config)
    if missing:
        raise ValueError(f"response config keys missing: {sorted(missing)}")
    return config


def _embedding_columns(path: Path) -> list[str]:
    import pyarrow.parquet as pq

    names = pq.ParquetFile(path).schema_arrow.names
    columns = [name for name in names if re.fullmatch(r"emb_\d+", str(name))]
    expected = [f"emb_{index}" for index in range(EXPECTED_DIM)]
    if columns != expected:
        raise ValueError(f"unexpected Prov embedding columns in {path}: {len(columns)}")
    for required in ("patient_id", "sample_id", "timepoint"):
        if required not in names:
            raise ValueError(f"embedding is missing {required}: {path}")
    return columns


def load_patient_mean(row: object) -> tuple[str, np.ndarray, dict]:
    import pyarrow.parquet as pq

    path = Path(str(row.prov_path))
    expected_hash = str(getattr(row, "prov_sha256", "")).strip()
    observed_hash = sha256_file(path)
    if expected_hash and observed_hash != expected_hash:
        raise ValueError(f"Prov embedding SHA-256 mismatch: {path}")
    columns = _embedding_columns(path)
    table = pq.read_table(path, columns=["patient_id", "sample_id", "timepoint", *columns])
    identity = table.select(["patient_id", "sample_id", "timepoint"]).to_pandas()
    patient_id = str(row.patient_id)
    sample_id = str(row.sample_id)
    if set(identity.patient_id.astype(str)) != {patient_id}:
        raise ValueError(f"patient identity mismatch: {path}")
    if set(identity.sample_id.astype(str)) != {sample_id}:
        raise ValueError(f"sample identity mismatch: {path}")
    if set(identity.timepoint.astype(str)) != {"pre"}:
        raise ValueError(f"non-pretreatment tile in {path}")
    values = table.select(columns).to_pandas().to_numpy(dtype=np.float64, copy=False)
    if values.shape[0] == 0 or values.shape[1] != EXPECTED_DIM or not np.isfinite(values).all():
        raise ValueError(f"invalid embedding matrix: {path} / {values.shape}")
    mean = values.mean(axis=0)
    return patient_id, mean, {
        "patient_id": patient_id,
        "sample_id": sample_id,
        "embedding_file": path.name,
        "embedding_sha256": observed_hash,
        "tile_count": int(len(values)),
        "embedding_dim": EXPECTED_DIM,
    }


def load_means(frame: pd.DataFrame, workers: int) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    means: dict[str, np.ndarray] = {}
    ledger: list[dict] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(load_patient_mean, row): str(row.patient_id) for row in frame.itertuples(index=False)}
        for future in as_completed(futures):
            patient_id, mean, record = future.result()
            if patient_id in means:
                raise ValueError(f"duplicate patient embedding: {patient_id}")
            means[patient_id] = mean
            ledger.append(record)
    if set(means) != set(frame.patient_id.astype(str)):
        raise ValueError("not every manifest patient produced a Prov mean")
    return means, pd.DataFrame(ledger).sort_values("patient_id").reset_index(drop=True)


def fit_model(features: np.ndarray, labels: np.ndarray) -> dict[str, object]:
    scaler = StandardScaler()
    pca = PCA(n_components=PCA_DIM, svd_solver="full", random_state=SEED)
    logistic = LogisticRegression(
        penalty="l2",
        C=LOGISTIC_C,
        solver="liblinear",
        random_state=SEED,
        max_iter=10_000,
    )
    reduced = pca.fit_transform(scaler.fit_transform(features))
    logistic.fit(reduced, labels)
    return {"scaler": scaler, "pca": pca, "logistic": logistic}


def decision(model: dict[str, object], features: np.ndarray) -> np.ndarray:
    scaled = model["scaler"].transform(features)
    reduced = model["pca"].transform(scaled)
    return np.asarray(model["logistic"].decision_function(reduced), dtype=float)


def sigmoid(logits: np.ndarray) -> np.ndarray:
    logits = np.asarray(logits, dtype=float)
    return np.where(
        logits >= 0,
        1.0 / (1.0 + np.exp(-logits)),
        np.exp(logits) / (1.0 + np.exp(logits)),
    )


def read_cohorts(config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    dev = pd.read_csv(config["development_manifest"], sep="\t", dtype=str)
    folds = pd.read_csv(config["development_folds"], sep="\t", dtype={"patient_id": str})
    validation = pd.read_csv(config["validation_manifest"], sep="\t", dtype=str)
    forbidden = [column for column in validation.columns if OUTCOME_LIKE.search(str(column))]
    if forbidden:
        raise ValueError(f"validation manifest contains outcome-like columns: {forbidden}")
    required_manifest = {"patient_id", "sample_id", "timepoint", "prov_path"}
    if not required_manifest.issubset(dev) or not required_manifest.issubset(validation):
        raise ValueError("manifest is missing a required identity/path column")
    required_folds = {"patient_id", "response", "fold"}
    if not required_folds.issubset(folds):
        raise ValueError("development fold table is incomplete")
    dev = dev.drop(columns=[column for column in ("response", "fold") if column in dev]).merge(
        folds[["patient_id", "response", "fold"]], on="patient_id", validate="one_to_one"
    )
    expected_dev = int(config["expected_development_patients"])
    expected_validation = int(config["expected_validation_patients"])
    if len(dev) != expected_dev or dev.patient_id.nunique() != expected_dev:
        raise ValueError("development cohort size/identity mismatch")
    if len(validation) != expected_validation or validation.patient_id.nunique() != expected_validation:
        raise ValueError("validation cohort size/identity mismatch")
    if set(dev.patient_id) & set(validation.patient_id):
        raise ValueError("development-validation patient overlap")
    if set(dev.response.astype(int).unique()) != {0, 1} or set(dev.fold.astype(int).unique()) != set(range(5)):
        raise ValueError("development response/fold contract failed")
    if set(dev.timepoint) != {"pre"} or set(validation.timepoint) != {"pre"}:
        raise ValueError("only pretreatment samples are permitted")
    sort_key = lambda series: series.str.extract(r"(\d+)")[0].astype(int)
    return (
        dev.sort_values("patient_id", key=sort_key).reset_index(drop=True),
        validation.sort_values("patient_id", key=sort_key).reset_index(drop=True),
    )


def run(config: dict, workers: int) -> None:
    output = Path(config["output_root"])
    output.mkdir(parents=True, exist_ok=True)
    development, validation = read_cohorts(config)
    combined = pd.concat([development, validation], ignore_index=True, sort=False)
    means, ledger = load_means(combined, workers)

    dev_ids = development.patient_id.astype(str).to_numpy()
    val_ids = validation.patient_id.astype(str).to_numpy()
    labels = development.response.astype(int).to_numpy()
    folds = development.fold.astype(int).to_numpy()
    dev_x = np.vstack([means[patient] for patient in dev_ids])
    val_x = np.vstack([means[patient] for patient in val_ids])

    oof_logits = np.full(len(development), np.nan, dtype=float)
    audit_rows: list[dict] = []
    for fold in range(5):
        held = folds == fold
        train = ~held
        train_ids = dev_ids[train].tolist()
        held_ids = dev_ids[held].tolist()
        if not held.any() or set(train_ids) & set(held_ids):
            raise ValueError(f"invalid patient split in fold {fold}")
        model = fit_model(dev_x[train], labels[train])
        oof_logits[held] = decision(model, dev_x[held])
        audit_rows.append(
            {
                "fold": fold,
                "train_patient_ids": train_ids,
                "held_patient_ids": held_ids,
                "held_excluded_from_scaler_pca_logistic": True,
            }
        )
    if not np.isfinite(oof_logits).all():
        raise ValueError("OOF predictions are incomplete")

    final_model = fit_model(dev_x, labels)
    validation_logits = decision(final_model, val_x)
    development_scores = pd.DataFrame(
        {
            "patient_id": dev_ids,
            "response": labels,
            "fold": folds,
            "prov_logit": oof_logits,
            "prov_score": sigmoid(oof_logits),
        }
    )
    validation_scores = pd.DataFrame(
        {
            "patient_id": val_ids,
            "sample_id": validation.sample_id.astype(str).to_numpy(),
            "prov_score": sigmoid(validation_logits),
        }
    )
    development_scores.to_csv(output / "development_oof_scores.tsv", sep="\t", index=False, float_format="%.17g")
    validation_scores.to_csv(
        output / "locked_validation_scores_outcome_free.tsv", sep="\t", index=False, float_format="%.17g"
    )
    ledger.to_csv(output / "input_embedding_ledger.tsv", sep="\t", index=False)
    (output / "fold_leakage_audit.json").write_text(
        json.dumps(audit_rows, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    bundle = {
        "models": {"prov": final_model},
        "selected_candidates": {
            "prov": {
                "candidate_id": "prov__mean__pca2__C0.01",
                "representation": "mean",
                "pca_dim": PCA_DIM,
                "C": LOGISTIC_C,
            }
        },
        "seed": SEED,
        "development_patient_ids": dev_ids.tolist(),
        "validation_patient_ids_used_for_fit_or_selection": [],
        "outcome_blind_validation": True,
    }
    joblib.dump(bundle, output / "fixed_prov_response_bundle.joblib", compress=3)
    artifacts = [
        output / "development_oof_scores.tsv",
        output / "locked_validation_scores_outcome_free.tsv",
        output / "input_embedding_ledger.tsv",
        output / "fold_leakage_audit.json",
        output / "fixed_prov_response_bundle.joblib",
    ]
    complete = {
        "status": "PASS",
        "model": "Prov mean -> StandardScaler -> PCA2 -> L2 logistic C=0.01",
        "seed": SEED,
        "development_patients": len(development),
        "validation_patients": len(validation),
        "validation_outcomes_accessed": False,
        "outputs": {path.name: sha256_file(path) for path in artifacts},
    }
    (output / "complete.json").write_text(json.dumps(complete, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def verify(config: dict) -> None:
    output = Path(config["output_root"])
    complete = json.loads((output / "complete.json").read_text(encoding="utf-8"))
    if complete.get("status") != "PASS" or complete.get("validation_outcomes_accessed") is not False:
        raise ValueError("response completion contract failed")
    for filename, expected in complete["outputs"].items():
        if sha256_file(output / filename) != expected:
            raise ValueError(f"response artifact hash mismatch: {filename}")
    development = pd.read_csv(output / "development_oof_scores.tsv", sep="\t")
    validation = pd.read_csv(output / "locked_validation_scores_outcome_free.tsv", sep="\t")
    if len(development) != int(config["expected_development_patients"]):
        raise ValueError("development output size mismatch")
    if len(validation) != int(config["expected_validation_patients"]):
        raise ValueError("validation output size mismatch")
    if OUTCOME_LIKE.search(" ".join(validation.columns)):
        raise ValueError("validation score output contains an outcome-like column")
    print("VERIFY PASS: response hashes, cohort sizes and outcome-blind validation schema")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("run", "verify"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    if not 1 <= args.workers <= 64:
        raise SystemExit("workers must be between 1 and 64")
    config = load_config(args.config)
    if args.stage == "run":
        run(config, args.workers)
    else:
        verify(config)


if __name__ == "__main__":
    main()
