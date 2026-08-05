#!/usr/bin/env python3
"""Build Figure 7e-h source data from the formal full-coverage analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import rankdata, spearmanr


CONFIG_ENV = "FIG7_SPATIAL_SOURCE_CONFIG"
RESULTS: Path
TEMPORAL_MANIFEST: Path
DEVELOPMENT_MANIFEST: Path
MEASURED_FOLDS: Path
M0_BUNDLE: Path
M0_SCORES: Path
TISSUE_ROOT: Path
EXTERNAL_MEASURED: Path
TEMPORAL_RESPONSE: Path
IDENTITY = ["patient_id", "sample_id", "tile_id", "x_min", "y_min", "x_max", "y_max"]
PHENOTYPES = ("Vascular", "Epithelial (PanCK + HER2)", "HER2")
BOOTSTRAP_SEED = 20260720
BOOTSTRAP_N = 10_000


def configure(path: Path) -> dict:
    """Load all private data locations from a path-free JSON configuration."""
    config = json.loads(path.read_text(encoding="utf-8"))
    required = {
        "regional_output_root",
        "temporal_manifest",
        "development_manifest",
        "measured_fold_manifest",
        "response_model_bundle",
        "response_validation_scores",
        "tissue_fraction_root",
        "external_measured_tile_density",
        "temporal_response_table",
    }
    missing = required - set(config)
    if missing:
        raise ValueError(f"spatial source config keys missing: {sorted(missing)}")
    global RESULTS, TEMPORAL_MANIFEST, DEVELOPMENT_MANIFEST, MEASURED_FOLDS
    global M0_BUNDLE, M0_SCORES, TISSUE_ROOT, EXTERNAL_MEASURED, TEMPORAL_RESPONSE
    RESULTS = Path(config["regional_output_root"])
    TEMPORAL_MANIFEST = Path(config["temporal_manifest"])
    DEVELOPMENT_MANIFEST = Path(config["development_manifest"])
    MEASURED_FOLDS = Path(config["measured_fold_manifest"])
    M0_BUNDLE = Path(config["response_model_bundle"])
    M0_SCORES = Path(config["response_validation_scores"])
    TISSUE_ROOT = Path(config["tissue_fraction_root"])
    EXTERNAL_MEASURED = Path(config["external_measured_tile_density"])
    TEMPORAL_RESPONSE = Path(config["temporal_response_table"])
    return config


if os.environ.get(CONFIG_ENV):
    configure(Path(os.environ[CONFIG_ENV]))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def feature_columns(frame: pd.DataFrame) -> list[str]:
    columns = sorted(
        (str(column) for column in frame.columns if str(column).startswith("emb_")),
        key=lambda value: int(value.split("_")[-1]),
    )
    if columns != [f"emb_{index}" for index in range(1536)]:
        raise ValueError(f"Prov feature dimension/order mismatch: {len(columns)}")
    return columns


def backproject_m0(model: dict[str, object]) -> tuple[np.ndarray, float]:
    scaler = model["scaler"]
    pca = model["pca"]
    logistic = model["logistic"]
    beta_scaled = pca.components_.T @ logistic.coef_.ravel()
    beta_raw = beta_scaled / scaler.scale_
    intercept = (
        float(logistic.intercept_[0])
        - float(pca.mean_ @ beta_scaled)
        - float((scaler.mean_ / scaler.scale_) @ beta_scaled)
    )
    return np.asarray(beta_raw, dtype=np.float64), intercept


def safe_rho(actual: np.ndarray, predicted: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    finite = np.isfinite(actual) & np.isfinite(predicted)
    if int(finite.sum()) < 3 or np.std(actual[finite]) == 0 or np.std(predicted[finite]) == 0:
        return float("nan")
    return float(spearmanr(actual[finite], predicted[finite]).statistic)


def spatial_metrics(contribution: np.ndarray, marker: np.ndarray, *, seed: int, permutations: int) -> dict:
    contribution = np.asarray(contribution, dtype=float)
    marker = np.asarray(marker, dtype=float)
    finite = np.isfinite(contribution) & np.isfinite(marker)
    contribution = contribution[finite]
    marker = marker[finite]
    if len(contribution) < 3:
        raise ValueError("Insufficient finite tiles for spatial metric")
    cr = rankdata(contribution).astype(float)
    mr = rankdata(marker).astype(float)
    cr -= cr.mean()
    mr -= mr.mean()
    denominator = float(np.linalg.norm(cr) * np.linalg.norm(mr))
    if denominator <= 0:
        raise ValueError("Constant spatial variable")
    observed = float(np.dot(cr, mr) / denominator)
    low = contribution <= np.quantile(contribution, 0.10)
    high = contribution >= np.quantile(contribution, 0.90)
    top_bottom = float(marker[high].mean() - marker[low].mean())
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(int(permutations)):
        null = float(np.dot(cr, rng.permutation(mr)) / denominator)
        exceed += int(null >= observed - 1e-15)
    return {
        "n_tiles": int(len(contribution)),
        "spearman": observed,
        "top_minus_bottom_decile": top_bottom,
        "permutation_p_one_sided": float((exceed + 1) / (int(permutations) + 1)),
        "n_permutations": int(permutations),
        "permutation_seed": int(seed),
    }


def bootstrap_median(values: np.ndarray, seed: int) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    rng = np.random.default_rng(seed)
    boot = np.median(values[rng.integers(0, len(values), size=(BOOTSTRAP_N, len(values)))], axis=1)
    return float(np.median(values)), float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def temporal_worker(row: dict, output_root: str, permutations: int, sample_index: int) -> dict:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    patient = str(row["patient_id"])
    sample = str(row["sample_id"])
    embeddings = pd.read_parquet(str(row["prov_path"]))
    embeddings["tile_id"] = embeddings.tile_id.astype(str)
    if embeddings.tile_id.duplicated().any():
        raise ValueError(f"Duplicate temporal embedding tile: {sample}")
    columns = feature_columns(embeddings)
    matrix = embeddings[columns].to_numpy(dtype=np.float64, copy=False)
    bundle = joblib.load(M0_BUNDLE)
    beta, intercept = backproject_m0(bundle["models"]["prov"])
    contribution = matrix @ beta
    reconstructed = float(expit(intercept + float(contribution.mean())))
    frozen = pd.read_csv(M0_SCORES, sep="\t")
    expected = float(frozen.loc[frozen.sample_id.astype(str).eq(sample), "prov_score"].iloc[0])
    error = abs(reconstructed - expected)
    if error > 1e-10:
        raise ValueError(f"M0 reconstruction mismatch: {sample} / {error}")

    maps = pd.read_parquet(RESULTS / "full_refit/temporal_maps" / f"{sample}.parquet")
    maps["tile_id"] = maps.tile_id.astype(str)
    if maps.tile_id.duplicated().any() or set(maps.tile_id) != set(embeddings.tile_id):
        raise ValueError(f"Regional-map/embedding tile mismatch: {sample}")
    maps = maps.set_index("tile_id").loc[embeddings.tile_id].reset_index()
    if not np.array_equal(maps[["x_min", "y_min", "x_max", "y_max"]].to_numpy(), embeddings[["x_min", "y_min", "x_max", "y_max"]].to_numpy()):
        raise ValueError(f"Regional-map/embedding coordinate mismatch: {sample}")
    fraction = pd.read_csv(TISSUE_ROOT / f"{sample}.tsv", sep="\t", usecols=["tile_id", "tissue_fraction_prov"])
    fraction["tile_id"] = fraction.tile_id.astype(str)
    if fraction.tile_id.duplicated().any():
        raise ValueError(f"Duplicate tissue fraction: {sample}")
    fraction = fraction.set_index("tile_id").loc[embeddings.tile_id].reset_index()
    tissue = pd.to_numeric(fraction.tissue_fraction_prov, errors="coerce").to_numpy(float)
    if not np.isfinite(tissue).all() or np.any(tissue < 0.01 - 1e-12):
        raise ValueError(f"Foreground threshold mismatch: {sample}")

    tiles = embeddings[IDENTITY].copy()
    tiles["tissue_fraction_prov"] = tissue
    tiles["m0_tile_contribution"] = contribution.astype(np.float32)
    tiles["Vascular_predicted"] = maps["Vascular_predicted"].to_numpy(np.float32)
    tiles["PanCK_predicted"] = maps["PanCK_predicted"].to_numpy(np.float32)
    tiles["HER2_predicted"] = maps["HER2_predicted"].to_numpy(np.float32)
    tiles["Epithelial (PanCK + HER2)_predicted"] = np.log1p(
        np.maximum(0.0, np.expm1(tiles.PanCK_predicted.to_numpy(float)))
        + np.maximum(0.0, np.expm1(tiles.HER2_predicted.to_numpy(float)))
    ).astype(np.float32)
    output_path = Path(output_root) / "temporal_tile_scores" / f"{sample}.parquet"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tiles.to_parquet(output_path, index=False)
    rows = []
    for marker_index, marker in enumerate(PHENOTYPES):
        rows.append(
            {
                "patient_id": patient,
                "sample_id": sample,
                "phenotype": marker,
                **spatial_metrics(
                    tiles.m0_tile_contribution.to_numpy(float),
                    tiles[f"{marker}_predicted"].to_numpy(float),
                    seed=20260727 + sample_index * 10 + marker_index,
                    permutations=permutations,
                ),
            }
        )
    return {
        "metrics": rows,
        "reconstruction": {
            "patient_id": patient, "sample_id": sample, "n_tiles": int(len(tiles)),
            "reconstructed_M0_score": reconstructed, "frozen_M0_score": expected,
            "absolute_error": error,
        },
        "ledger": {
            "patient_id": patient, "sample_id": sample,
            "embedding_path": str(row["prov_path"]), "embedding_sha256": sha256_file(Path(str(row["prov_path"]))),
            "regional_map_path": str(RESULTS / "full_refit/temporal_maps" / f"{sample}.parquet"),
            "all_foreground_tiles_used": True, "tile_count": int(len(tiles)),
            "output_path": str(output_path), "output_sha256": sha256_file(output_path),
        },
    }


def build_panel_e(output: Path) -> None:
    summary = pd.read_csv(RESULTS / "regional_summary.tsv", sep="\t")
    predictions = pd.read_parquet(RESULTS / "oof_predictions.parquet", columns=["regime", "tile_id", "marker"])
    combined = pd.read_parquet(RESULTS / "combined_epithelial_oof.parquet", columns=["regime", "tile_id", "marker"])
    selected = summary.loc[summary.regime.eq("primary") & summary.marker.isin(PHENOTYPES)].copy()
    selected["phenotype"] = selected.marker.map(
        {"Vascular": "Vascular", "Epithelial (PanCK + HER2)": "PanCK signal", "HER2": "HER2"}
    )
    selected["analysis_status"] = np.where(selected.marker.eq("Epithelial (PanCK + HER2)"), "post_hoc_combined", "prespecified_component")
    counts = predictions.loc[predictions.regime.eq("primary")].groupby("marker").tile_id.count().to_dict()
    counts["Epithelial (PanCK + HER2)"] = int(combined.loc[combined.regime.eq("primary"), "tile_id"].count())
    selected["n_valid_oof_tiles"] = selected.marker.map(counts).astype(int)
    selected.to_csv(output / "panel_e_regional.tsv", sep="\t", index=False)
    summary.loc[summary.regime.eq("primary")].to_csv(output / "panel_e_all_components_audit.tsv", sep="\t", index=False)


def build_panel_f(output: Path, workers: int, permutations: int) -> None:
    temporal = pd.read_csv(TEMPORAL_MANIFEST, sep="\t", dtype=str).sort_values("patient_id").reset_index(drop=True)
    if any(str(column).lower() in {"response", "outcome", "label", "target"} for column in temporal.columns):
        raise ValueError("Temporal outcome guard failed")
    results = []
    with ProcessPoolExecutor(max_workers=min(int(workers), len(temporal))) as pool:
        futures = {
            pool.submit(temporal_worker, row._asdict(), str(output), int(permutations), index): str(row.sample_id)
            for index, row in enumerate(temporal.itertuples(index=False))
        }
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"PANEL_F_PASS {futures[future]}", flush=True)
    metrics = pd.DataFrame([row for result in results for row in result["metrics"]]).sort_values(["phenotype", "patient_id"])
    reconstruction = pd.DataFrame([result["reconstruction"] for result in results]).sort_values("patient_id")
    ledger = pd.DataFrame([result["ledger"] for result in results]).sort_values("patient_id")
    if reconstruction.absolute_error.max() > 1e-10 or not ledger.all_foreground_tiles_used.all():
        raise ValueError("Panel f M0/full-foreground audit failed")
    summaries = []
    for index, phenotype in enumerate(PHENOTYPES):
        group = metrics.loc[metrics.phenotype.eq(phenotype)]
        median, lower, upper = bootstrap_median(group.spearman.to_numpy(float), BOOTSTRAP_SEED + index)
        summaries.append(
            {
                "phenotype": phenotype,
                "n_patients": int(group.patient_id.nunique()),
                "median_patient_spearman": median,
                "bootstrap_95ci_lower": lower,
                "bootstrap_95ci_upper": upper,
                "positive_patient_fraction": float(np.mean(group.spearman.to_numpy(float) > 0)),
                "median_top_minus_bottom_decile": float(group.top_minus_bottom_decile.median()),
                "analysis_status": "post_hoc_combined" if phenotype.startswith("Epithelial") else "prespecified_component",
            }
        )
    metrics.to_csv(output / "panel_f_patient_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(output / "panel_f_summary.tsv", sep="\t", index=False)
    reconstruction.to_csv(output / "panel_f_M0_reconstruction.tsv", sep="\t", index=False)
    ledger.to_csv(output / "panel_f_input_ledger.tsv", sep="\t", index=False)


def build_panel_g(output: Path) -> pd.DataFrame:
    measured_all = pd.read_csv(EXTERNAL_MEASURED, sep="\t", low_memory=False)
    responses = pd.read_csv(TEMPORAL_RESPONSE, sep="\t", dtype={"patient_id": str}).set_index("patient_id")
    pieces, concordance = [], []
    for patient in ("pt71", "pt80"):
        sample = f"{patient}-pre"
        valid_geometry = measured_all.geometry_valid
        if valid_geometry.dtype != bool:
            valid_geometry = valid_geometry.astype(str).str.lower().isin(["true", "1"])
        measured = measured_all.loc[
            measured_all.sample_id.astype(str).eq(sample)
            & measured_all.marker.astype(str).eq("Vascular")
            & valid_geometry
        ].copy()
        maps = pd.read_parquet(RESULTS / "full_refit/temporal_maps" / f"{sample}.parquet")
        measured["tile_id"] = measured.tile_id.astype(str)
        maps["tile_id"] = maps.tile_id.astype(str)
        joined = measured.drop(columns=["predicted"], errors="ignore").merge(
            maps[["tile_id", "x_min", "y_min", "x_max", "y_max", "Vascular_predicted"]],
            on=["tile_id", "x_min", "y_min", "x_max", "y_max"], how="inner", validate="one_to_one",
        ).rename(columns={"Vascular_predicted": "predicted"})
        if len(joined) != len(measured):
            raise ValueError(f"External Vascular geometry mismatch: {sample} / {len(joined)} != {len(measured)}")
        joined["response"] = int(responses.loc[patient, "response"])
        joined["response_label"] = joined.response.map({0: "NR", 1: "CR"})
        joined["source"] = "external_measured_same_section_plus_fullcoverage_prediction"
        concordance.append(
            {"panel": "g", "patient_id": patient, "sample_id": sample, "phenotype": "Vascular",
             "n_tiles": int(len(joined)), "spearman": safe_rho(joined.measured, joined.predicted)}
        )
        pieces.append(joined)
    panel = pd.concat(pieces, ignore_index=True)
    panel.to_csv(output / "panel_g_external_vascular.tsv.gz", sep="\t", index=False)
    return pd.DataFrame(concordance)


def build_panel_h(output: Path) -> pd.DataFrame:
    oof = pd.read_parquet(RESULTS / "oof_predictions.parquet")
    her2 = oof.loc[oof.regime.eq("primary") & oof.marker.eq("HER2") & oof.control.eq("none")].copy()
    folds = pd.read_csv(MEASURED_FOLDS, sep="\t", dtype={"patient_id": str})
    manifest = pd.read_csv(DEVELOPMENT_MANIFEST, sep="\t", dtype={"patient_id": str})
    response = manifest.loc[manifest.patient_id.astype(str).isin(set(folds.patient_id.astype(str))), ["patient_id", "response"]].copy()
    response["response"] = response.response.astype(int)
    rows = []
    for (patient, sample), group in her2.groupby(["patient_id", "sample_id"], sort=True):
        rows.append(
            {"patient_id": str(patient), "sample_id": str(sample), "n_tiles": int(len(group)),
             "spearman": safe_rho(group.actual, group.prediction)}
        )
    selection = pd.DataFrame(rows).merge(response, on="patient_id", validate="one_to_one")
    selection["response_label"] = selection.response.map({0: "NR", 1: "CR"})
    selection["group_median_spearman"] = selection.groupby("response").spearman.transform("median")
    selection["distance_to_group_median"] = (selection.spearman - selection.group_median_spearman).abs()
    selection["selected"] = False
    for group in (0, 1):
        candidates = selection.loc[selection.response.eq(group) & selection.n_tiles.ge(1000)].sort_values(
            ["distance_to_group_median", "n_tiles", "patient_id"], ascending=[True, False, True], kind="mergesort"
        )
        selection.loc[candidates.index[0], "selected"] = True
    selection["selection_rule"] = "within NR/CR, closest to group-median corrected OOF HER2 Spearman among n>=1000; ties larger n then patient ID"
    selection.to_csv(output / "panel_h_internal_case_selection.tsv", sep="\t", index=False)
    selected = selection.loc[selection.selected].copy()
    pieces, concordance = [], []
    for row in selected.sort_values("response").itertuples(index=False):
        frame = her2.loc[her2.patient_id.astype(str).eq(str(row.patient_id))].copy()
        frame["response"] = int(row.response)
        frame["response_label"] = str(row.response_label)
        frame["measured"] = frame.actual.astype(float)
        frame["predicted"] = frame.prediction.astype(float)
        frame["source"] = "patient_heldout_fullcoverage_map"
        pieces.append(frame)
        concordance.append(
            {"panel": "h", "patient_id": str(row.patient_id), "sample_id": str(row.sample_id),
             "phenotype": "HER2", "n_tiles": int(row.n_tiles), "spearman": float(row.spearman)}
        )
    panel = pd.concat(pieces, ignore_index=True)
    panel.to_csv(output / "panel_h_internal_her2.tsv.gz", sep="\t", index=False)
    return pd.DataFrame(concordance)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--workers", type=int, default=13)
    parser.add_argument("--permutations", type=int, default=999)
    args = parser.parse_args()
    if not 1 <= args.workers <= 64:
        raise ValueError("workers must be between 1 and 64")
    if args.permutations < 1:
        raise ValueError("permutations must be positive")
    resolved_config = args.config.resolve()
    os.environ[CONFIG_ENV] = str(resolved_config)
    config = configure(resolved_config)
    output = args.output or Path(config.get("source_output_root", RESULTS / "figure7_source_data"))
    output.mkdir(parents=True, exist_ok=True)
    build_panel_e(output)
    build_panel_f(output, args.workers, args.permutations)
    concordance = pd.concat([build_panel_g(output), build_panel_h(output)], ignore_index=True)
    concordance.to_csv(output / "panel_gh_concordance.tsv", sep="\t", index=False)
    scales = []
    for panel_name, path, phenotype in (
        ("g", output / "panel_g_external_vascular.tsv.gz", "Vascular"),
        ("h", output / "panel_h_internal_her2.tsv.gz", "HER2"),
    ):
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        values = np.concatenate([frame.measured.to_numpy(float), frame.predicted.to_numpy(float)])
        values = values[np.isfinite(values)]
        scales.append(
            {"panel": panel_name, "phenotype": phenotype, "vmin": float(np.quantile(values, 0.01)),
             "vmax": float(np.quantile(values, 0.99)), "n_values": int(len(values)),
             "scale_rule": "q01-q99 finite measured+predicted union across displayed cases"}
        )
    pd.DataFrame(scales).to_csv(output / "panel_gh_scales.tsv", sep="\t", index=False)
    direct_inputs = [
        RESULTS / "regional_summary.tsv", RESULTS / "oof_predictions.parquet",
        RESULTS / "control_oof_predictions.parquet", RESULTS / "full_refit/complete.json",
        M0_BUNDLE, M0_SCORES, EXTERNAL_MEASURED, TEMPORAL_MANIFEST, DEVELOPMENT_MANIFEST,
    ]
    hashes = pd.DataFrame({"path": list(map(str, direct_inputs)), "sha256": [sha256_file(path) for path in direct_inputs]})
    hashes.to_csv(output / "direct_input_hashes.tsv", sep="\t", index=False)
    output_files = sorted(path for path in output.iterdir() if path.is_file())
    status = {
        "status": "PASS", "panel_e_patients": 30, "panel_f_patients": 13,
        "panel_g_patients": ["pt71", "pt80"],
        "panel_h_patients": sorted(concordance.loc[concordance.panel.eq("h"), "patient_id"].tolist()),
        "panel_f_permutations_per_patient_phenotype": int(args.permutations),
        "outcome_used_for_model_fit": False,
        "regional_training_target": "coordinate-corrected valid H&E-native Panel-2 targets",
        "outputs": {path.name: sha256_file(path) for path in output_files},
    }
    (output / "complete.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
