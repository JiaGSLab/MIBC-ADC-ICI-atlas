"""Core primitives for full-coverage regional phenotype modelling.

Training is spatially balanced and capped per patient; evaluation and temporal
inference are deliberately uncapped. All joins use tile ID and exact H&E
level-0 coordinates.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


COORDS = ["x_min", "y_min", "x_max", "y_max"]
IDENTITY = ["tile_id", *COORDS]
MARKERS = ("Vascular", "PanCK", "HER2")
DEFAULT_TARGET_COLUMNS = {
    "Vascular": "Vascular_log1p_cells_per_mm2",
    "PanCK": "PanCK+ Tumor_log1p_cells_per_mm2",
    "HER2": "HER2+ Tumor_log1p_cells_per_mm2",
}


@dataclass
class PatientData:
    patient_id: str
    sample_id: str
    outer_fold: int
    mpp: float
    features: np.ndarray
    tile_id: np.ndarray
    coords: np.ndarray
    bins: np.ndarray
    targets: np.ndarray

    @property
    def n_tiles(self) -> int:
        return int(len(self.tile_id))


def seed32(*tokens: object) -> int:
    payload = "|".join(map(str, tokens)).encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:8], 16) % (2**32)


def feature_columns(columns: list[str] | pd.Index, expected_dim: int = 1536) -> list[str]:
    selected = sorted(
        (str(column) for column in columns if str(column).startswith("emb_")),
        key=lambda value: int(value.split("_")[-1]),
    )
    expected = [f"emb_{index}" for index in range(int(expected_dim))]
    if selected != expected:
        raise ValueError(f"embedding feature schema mismatch: {len(selected)} != {expected_dim}")
    return selected


def load_patient(
    manifest_row: object,
    outer_fold: int,
    targets: pd.DataFrame,
    *,
    target_columns: dict[str, str] | None = None,
    expected_dim: int = 1536,
    bin_um: float = 2000.0,
) -> PatientData:
    target_columns = dict(target_columns or DEFAULT_TARGET_COLUMNS)
    if tuple(target_columns) != MARKERS:
        raise ValueError(f"target marker order must be {MARKERS}")
    patient_id = str(manifest_row.patient_id)
    sample_id = str(manifest_row.sample_id)
    mpp = float(manifest_row.mpp)
    if not np.isfinite(mpp) or mpp <= 0:
        raise ValueError(f"invalid mpp for {sample_id}: {mpp}")
    path = Path(str(manifest_row.prov_path))
    if not path.is_file():
        raise FileNotFoundError(path)
    embedding = pd.read_parquet(path)
    fcols = feature_columns(embedding.columns, expected_dim)
    embedding["tile_id"] = embedding.tile_id.astype(str)
    if embedding.tile_id.duplicated().any():
        raise ValueError(f"duplicate embedding tile IDs: {sample_id}")
    identity = embedding[IDENTITY].copy()
    identity["embedding_row"] = np.arange(len(identity), dtype=np.int64)

    columns = ["sample_id", *IDENTITY, *target_columns.values()]
    if "geometry_valid" in targets.columns:
        columns.append("geometry_valid")
    sample_targets = targets.loc[targets.sample_id.astype(str).eq(sample_id), columns].copy()
    sample_targets["tile_id"] = sample_targets.tile_id.astype(str)
    if "geometry_valid" in sample_targets.columns:
        valid = sample_targets.geometry_valid
        if valid.dtype != bool:
            valid = valid.astype(str).str.lower().isin(["true", "1"])
        sample_targets = sample_targets.loc[valid].copy()
    values = sample_targets[list(target_columns.values())].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    sample_targets = sample_targets.loc[np.isfinite(values).all(axis=1)].copy()
    if sample_targets.duplicated(IDENTITY).any() or sample_targets.tile_id.duplicated().any():
        raise ValueError(f"duplicate target identity: {sample_id}")

    joined = identity.merge(
        sample_targets.drop(columns=["sample_id", "geometry_valid"], errors="ignore"),
        on=IDENTITY,
        how="inner",
        validate="one_to_one",
    )
    if joined.empty:
        raise ValueError(f"no exact foreground/target overlap: {sample_id}")
    if len(joined) != len(sample_targets):
        raise ValueError(
            f"valid H&E-native target lacks exact foreground embedding: "
            f"{sample_id} / {len(joined)} != {len(sample_targets)}"
        )
    joined = joined.sort_values(["y_min", "x_min", "tile_id"], kind="mergesort").reset_index(drop=True)
    row_index = joined.embedding_row.to_numpy(np.int64)
    x = embedding.iloc[row_index][fcols].to_numpy(dtype=np.float32, copy=True)
    if x.shape != (len(joined), int(expected_dim)) or not np.isfinite(x).all():
        raise ValueError(f"invalid feature matrix: {sample_id} / {x.shape}")
    coords = joined[COORDS].to_numpy(np.int32, copy=True)
    centres = np.column_stack(((coords[:, 0] + coords[:, 2]) / 2, (coords[:, 1] + coords[:, 3]) / 2))
    bins = np.floor(centres * mpp / float(bin_um)).astype(np.int32)
    marker_targets = np.column_stack(
        [pd.to_numeric(joined[target_columns[marker]], errors="coerce").to_numpy(np.float32) for marker in MARKERS]
    )
    return PatientData(
        patient_id=patient_id,
        sample_id=sample_id,
        outer_fold=int(outer_fold),
        mpp=mpp,
        features=x,
        tile_id=joined.tile_id.astype(str).to_numpy(),
        coords=coords,
        bins=bins,
        targets=marker_targets,
    )


def spatially_balanced_indices(patient: PatientData, cap: int, seed: int) -> np.ndarray:
    """Round-robin occupied physical bins after deterministic within-bin shuffling."""
    if int(cap) <= 0:
        raise ValueError("training cap must be positive")
    groups: dict[tuple[int, int], list[int]] = {}
    for index, region in enumerate(patient.bins):
        groups.setdefault((int(region[0]), int(region[1])), []).append(index)
    queues: dict[tuple[int, int], np.ndarray] = {}
    for region in sorted(groups):
        indices = np.asarray(sorted(groups[region], key=lambda i: str(patient.tile_id[i])), dtype=np.int64)
        rng = np.random.default_rng(seed32("regional-balanced", seed, patient.patient_id, *region))
        queues[region] = indices[rng.permutation(len(indices))]
    maximum = min(int(cap), patient.n_tiles)
    selected: list[int] = []
    cursor = 0
    while len(selected) < maximum:
        advanced = False
        for region in sorted(queues):
            queue = queues[region]
            if cursor < len(queue) and len(selected) < maximum:
                selected.append(int(queue[cursor]))
                advanced = True
        if not advanced:
            break
        cursor += 1
    result = np.asarray(selected, dtype=np.int64)
    if len(result) != maximum or len(np.unique(result)) != len(result):
        raise RuntimeError(f"balanced sampling failed: {patient.patient_id}")
    return result


def permute_within_patient(
    values: np.ndarray,
    patient_ids: np.ndarray,
    *,
    marker: str,
    seed: int,
) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64).copy()
    patient_ids = patient_ids.astype(str)
    for patient_id in sorted(set(patient_ids)):
        indices = np.flatnonzero(patient_ids == patient_id)
        finite = indices[np.isfinite(result[indices])]
        rng = np.random.default_rng(seed32("regional-control", seed, marker, patient_id))
        result[finite] = result[finite][rng.permutation(len(finite))]
    return result


def safe_spearman(actual: np.ndarray, predicted: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    finite = np.isfinite(actual) & np.isfinite(predicted)
    if int(finite.sum()) < 3 or np.std(actual[finite]) == 0 or np.std(predicted[finite]) == 0:
        return float("nan")
    return float(spearmanr(actual[finite], predicted[finite]).statistic)


def combined_epithelial(log_panck: np.ndarray, log_her2: np.ndarray) -> np.ndarray:
    """Combine overlapping marker-density channels; this is not a unique-cell count."""
    return np.log1p(
        np.maximum(0.0, np.expm1(np.asarray(log_panck, dtype=float)))
        + np.maximum(0.0, np.expm1(np.asarray(log_her2, dtype=float)))
    )

