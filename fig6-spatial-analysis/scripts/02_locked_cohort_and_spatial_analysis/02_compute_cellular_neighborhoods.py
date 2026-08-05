#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute sample-wise radius-neighborhood celltype composition and cluster valid cells.

Designed for the output produced by phenotype_filtered_h5ad_no_csv_v2.py when
--write_combined_table was used. The input can be:
    - all_samples_cell_phenotypes.pkl
    - all_samples_cell_phenotypes.parquet
    - a CSV file if needed

Core behavior
-------------
1. Group cells by sample.
2. Use spatial_x / spatial_y as registered coordinates.
3. For every cell, find neighbors within a fixed radius.
4. Remove the center cell itself from its neighborhood.
5. Count neighbor celltypes and convert counts to proportions.
6. Record neighbor_count for every cell.
7. Mark valid cells using min/max neighbor count filters.
8. Cluster only valid cells using their neighborhood composition vectors.

Outputs
-------
output_dir/
    neighborhood_composition.pkl or .parquet
    neighborhood_count_summary.csv
    k_selection_scores.csv
    neighborhood_cluster_centers.csv
    run_params.csv
    sample_neighborhood_manifest.csv
    sample_neighborhoods/
    plots/
        neighbor_count_hist_overall.png
        neighbor_count_boxplot_by_sample.png
        k_selection_calinski_harabasz.png

Notes
-----
- Neighborhood is computed within each sample only. Multiple samples can be processed in parallel.
- Clustering is global by default, so CN labels are comparable across samples.
- Center cells are not counted as their own neighbors.
"""

from __future__ import annotations

import argparse
import re
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.metrics import calinski_harabasz_score
from sklearn.preprocessing import StandardScaler


def log(msg: str):
    print(msg, flush=True)


def safe_name(x: str) -> str:
    x = str(x)
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", x)
    return x.strip("_") or "sample"


def read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in [".pkl", ".pickle"]:
        return pd.read_pickle(path)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported input format: {path}")


def read_written_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in [".pkl", ".pickle"]:
        return pd.read_pickle(path)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported intermediate table: {path}")


def write_table(df: pd.DataFrame, path_without_suffix: Path, prefer: str = "parquet") -> Path:
    path_without_suffix = Path(path_without_suffix)
    if prefer == "parquet":
        out = path_without_suffix.with_suffix(".parquet")
        try:
            df.to_parquet(out, index=False)
            return out
        except Exception as e:
            log(f"[WARN] Could not write parquet ({e}); falling back to pickle.")
    out = path_without_suffix.with_suffix(".pkl")
    df.to_pickle(out)
    return out


def write_small_csv(df: pd.DataFrame, path_without_suffix: Path) -> Path:
    """
    Small summary/statistic tables are written as CSV for quick inspection.
    """
    out = Path(path_without_suffix).with_suffix(".csv")
    df.to_csv(out, index=False)
    return out


def prepare_input_table(
    df: pd.DataFrame,
    sample_col: str = "sample",
    cellid_col: str = "cellid",
    x_col: str = "spatial_x",
    y_col: str = "spatial_y",
    phenotype_col: str = "celltype",
    drop_undefined: bool = False,
) -> pd.DataFrame:
    required = [sample_col, x_col, y_col, phenotype_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    out = df.copy()

    if cellid_col not in out.columns:
        log(f"[WARN] cellid column '{cellid_col}' not found. Using dataframe index as cellid.")
        out[cellid_col] = out.index.astype(str)

    out[sample_col] = out[sample_col].astype(str)
    out[cellid_col] = out[cellid_col].astype(str)
    out[phenotype_col] = out[phenotype_col].astype(str)

    out[x_col] = pd.to_numeric(out[x_col], errors="coerce")
    out[y_col] = pd.to_numeric(out[y_col], errors="coerce")

    before = len(out)
    out = out[np.isfinite(out[x_col]) & np.isfinite(out[y_col])].copy()
    if len(out) < before:
        log(f"[WARN] Removed {before - len(out)} cells with non-finite coordinates.")

    if drop_undefined:
        before = len(out)
        out = out[out[phenotype_col] != "undefined"].copy()
        log(f"[INFO] Removed {before - len(out)} undefined cells before neighborhood counting.")

    if out[cellid_col].duplicated().any():
        n_dup = int(out[cellid_col].duplicated().sum())
        log(f"[WARN] {n_dup} duplicated cellid values detected. Making them unique.")
        ranks = out.groupby(cellid_col, observed=False).cumcount().astype(str)
        dup_mask = out[cellid_col].duplicated(keep=False)
        out.loc[dup_mask, cellid_col] = out.loc[dup_mask, cellid_col] + "__dup" + ranks.loc[dup_mask]

    return out.reset_index(drop=True)


def compute_sample_neighborhood_composition(
    sample_df: pd.DataFrame,
    phenotype_categories: Sequence[str],
    radius: float,
    sample_col: str = "sample",
    cellid_col: str = "cellid",
    x_col: str = "spatial_x",
    y_col: str = "spatial_y",
    phenotype_col: str = "celltype",
    chunk_size: int = 20000,
) -> pd.DataFrame:
    """
    Compute radius-neighborhood phenotype proportions for one sample.

    Center cell itself is removed from its own neighborhood.
    """
    if sample_df.empty:
        return pd.DataFrame()

    sample_name = str(sample_df[sample_col].iloc[0])
    coords = sample_df[[x_col, y_col]].to_numpy(dtype=np.float64)
    phenotypes = sample_df[phenotype_col].astype(str).to_numpy()

    ph_to_code = {p: i for i, p in enumerate(phenotype_categories)}
    codes = np.array([ph_to_code.get(p, -1) for p in phenotypes], dtype=np.int32)
    n_ph = len(phenotype_categories)
    n = len(sample_df)

    tree = cKDTree(coords)

    count_matrix = np.zeros((n, n_ph), dtype=np.float32)
    neighbor_count = np.zeros(n, dtype=np.int32)

    log(f"[NEIGHBOR] sample={sample_name}, cells={n}, radius={radius}")

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        query_coords = coords[start:end]
        neighbor_lists = tree.query_ball_point(query_coords, r=radius)

        for offset, neigh in enumerate(neighbor_lists):
            i = start + offset
            if len(neigh) == 0:
                continue

            neigh_arr = np.asarray(neigh, dtype=np.int64)
            neigh_arr = neigh_arr[neigh_arr != i]  # remove self

            if neigh_arr.size == 0:
                continue

            neigh_codes = codes[neigh_arr]
            neigh_codes = neigh_codes[neigh_codes >= 0]
            if neigh_codes.size == 0:
                continue

            counts = np.bincount(neigh_codes, minlength=n_ph).astype(np.float32)
            count_matrix[i, :] = counts
            neighbor_count[i] = int(counts.sum())

    prop_matrix = count_matrix.copy()
    denom = neighbor_count.astype(np.float32)
    valid = denom > 0
    prop_matrix[valid, :] = prop_matrix[valid, :] / denom[valid, None]

    comp_cols = [f"neigh_prop__{p}" for p in phenotype_categories]
    count_cols = [f"neigh_count__{p}" for p in phenotype_categories]

    out = sample_df[[sample_col, cellid_col, x_col, y_col, phenotype_col]].copy()
    out["neighbor_count"] = neighbor_count

    prop_df = pd.DataFrame(prop_matrix, columns=comp_cols, index=out.index)
    count_df = pd.DataFrame(count_matrix, columns=count_cols, index=out.index)

    return pd.concat([out, prop_df, count_df], axis=1)


def _sample_worker(
    sample: str,
    sample_df: pd.DataFrame,
    phenotype_order: List[str],
    radius: float,
    sample_col: str,
    cellid_col: str,
    x_col: str,
    y_col: str,
    phenotype_col: str,
    chunk_size: int,
    sample_output_dir: str,
    prefer_output: str,
) -> Dict:
    """
    Worker for one sample. Writes one intermediate table to disk and returns only metadata.
    This avoids returning large DataFrames through multiprocessing pipes.
    """
    try:
        res = compute_sample_neighborhood_composition(
            sample_df=sample_df.reset_index(drop=True),
            phenotype_categories=phenotype_order,
            radius=radius,
            sample_col=sample_col,
            cellid_col=cellid_col,
            x_col=x_col,
            y_col=y_col,
            phenotype_col=phenotype_col,
            chunk_size=chunk_size,
        )
        out_base = Path(sample_output_dir) / safe_name(sample)
        out_path = write_table(res, out_base, prefer=prefer_output)
        return {
            "sample": sample,
            "status": "ok",
            "n_cells": int(len(sample_df)),
            "output_path": str(out_path),
            "message": "",
        }
    except Exception as e:
        return {
            "sample": sample,
            "status": "error",
            "n_cells": int(len(sample_df)),
            "output_path": "",
            "message": str(e) + "\n" + traceback.format_exc(),
        }


def compute_neighborhood_composition_by_sample_parallel(
    df: pd.DataFrame,
    radius: float,
    output_dir: Path,
    sample_col: str = "sample",
    cellid_col: str = "cellid",
    x_col: str = "spatial_x",
    y_col: str = "spatial_y",
    phenotype_col: str = "celltype",
    phenotype_order: Optional[List[str]] = None,
    chunk_size: int = 20000,
    workers: int = 1,
    prefer_output: str = "parquet",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute sample-wise radius-neighborhood composition.

    If workers > 1, samples are processed in parallel. Each worker writes an
    intermediate per-sample result to output_dir/sample_neighborhoods/ and
    returns only file metadata to reduce IPC overhead.
    """
    if phenotype_order is None:
        phenotype_order = sorted(df[phenotype_col].dropna().astype(str).unique().tolist())

    sample_output_dir = Path(output_dir) / "sample_neighborhoods"
    sample_output_dir.mkdir(parents=True, exist_ok=True)

    sample_items = [(str(s), sub.copy()) for s, sub in df.groupby(sample_col, sort=True)]
    log(f"[INFO] samples for neighborhood computation: {len(sample_items)}")
    log(f"[INFO] neighborhood workers: {workers}")

    records = []
    if workers <= 1:
        for sample, sample_df in sample_items:
            rec = _sample_worker(
                sample=sample,
                sample_df=sample_df,
                phenotype_order=phenotype_order,
                radius=radius,
                sample_col=sample_col,
                cellid_col=cellid_col,
                x_col=x_col,
                y_col=y_col,
                phenotype_col=phenotype_col,
                chunk_size=chunk_size,
                sample_output_dir=str(sample_output_dir),
                prefer_output=prefer_output,
            )
            records.append(rec)
            log(f"[{rec['status'].upper()}] {sample}: {rec.get('output_path', '')}")
            if rec["status"] != "ok":
                log(rec["message"])
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futures = {
                ex.submit(
                    _sample_worker,
                    sample,
                    sample_df,
                    phenotype_order,
                    radius,
                    sample_col,
                    cellid_col,
                    x_col,
                    y_col,
                    phenotype_col,
                    chunk_size,
                    str(sample_output_dir),
                    prefer_output,
                ): sample
                for sample, sample_df in sample_items
            }
            for fut in as_completed(futures):
                sample = futures[fut]
                try:
                    rec = fut.result()
                except Exception as e:
                    rec = {
                        "sample": sample,
                        "status": "worker_error",
                        "n_cells": np.nan,
                        "output_path": "",
                        "message": str(e),
                    }
                records.append(rec)
                log(f"[{rec['status'].upper()}] {sample}: {rec.get('output_path', '')}")
                if rec["status"] != "ok":
                    log(rec["message"])

    manifest = pd.DataFrame(records).sort_values("sample").reset_index(drop=True)
    manifest_path = write_small_csv(manifest, Path(output_dir) / "sample_neighborhood_manifest")
    log(f"[OUTPUT] sample manifest: {manifest_path}")

    ok_paths = manifest.loc[manifest["status"] == "ok", "output_path"].tolist()
    if not ok_paths:
        raise RuntimeError("No sample-level neighborhood tables were successfully generated.")

    log(f"[COMBINE] reading {len(ok_paths)} sample-level neighborhood tables")
    comp_df = pd.concat([read_written_table(Path(p)) for p in ok_paths], axis=0, ignore_index=True)
    return comp_df, manifest


def add_valid_neighborhood_flag(
    comp_df: pd.DataFrame,
    min_neighbors: int = 10,
    max_neighbors: Optional[int] = None,
    min_neighbors_quantile: Optional[float] = None,
    max_neighbors_quantile: Optional[float] = None,
) -> Tuple[pd.DataFrame, Dict]:
    out = comp_df.copy()
    counts = out["neighbor_count"].astype(float)

    effective_min = int(min_neighbors)
    effective_max = max_neighbors

    if min_neighbors_quantile is not None:
        qmin = int(np.nanquantile(counts, min_neighbors_quantile))
        effective_min = max(effective_min, qmin)

    if max_neighbors_quantile is not None:
        qmax = int(np.nanquantile(counts, max_neighbors_quantile))
        effective_max = qmax if effective_max is None else min(effective_max, qmax)

    valid = counts >= effective_min
    if effective_max is not None:
        valid &= counts <= effective_max

    out["valid_neighborhood"] = valid.astype(bool)

    params = {
        "min_neighbors": int(min_neighbors),
        "max_neighbors": None if max_neighbors is None else int(max_neighbors),
        "min_neighbors_quantile": min_neighbors_quantile,
        "max_neighbors_quantile": max_neighbors_quantile,
        "effective_min_neighbors": int(effective_min),
        "effective_max_neighbors": None if effective_max is None else int(effective_max),
        "n_cells": int(len(out)),
        "n_valid_cells": int(valid.sum()),
        "valid_fraction": float(valid.mean()) if len(valid) else np.nan,
    }
    return out, params


def summarize_neighbor_counts(comp_df: pd.DataFrame, sample_col: str = "sample") -> pd.DataFrame:
    def q(x, p):
        return np.nanquantile(x, p) if len(x) else np.nan

    records = []
    for sample, sub in comp_df.groupby(sample_col, sort=True):
        x = sub["neighbor_count"].astype(float).to_numpy()
        records.append({
            "sample": sample,
            "n_cells": len(sub),
            "neighbor_count_min": np.nanmin(x) if len(x) else np.nan,
            "neighbor_count_q01": q(x, 0.01),
            "neighbor_count_q05": q(x, 0.05),
            "neighbor_count_q10": q(x, 0.10),
            "neighbor_count_median": q(x, 0.50),
            "neighbor_count_q90": q(x, 0.90),
            "neighbor_count_q95": q(x, 0.95),
            "neighbor_count_q99": q(x, 0.99),
            "neighbor_count_max": np.nanmax(x) if len(x) else np.nan,
            "neighbor_count_mean": np.nanmean(x) if len(x) else np.nan,
            "n_valid_cells": int(sub["valid_neighborhood"].sum()) if "valid_neighborhood" in sub.columns else np.nan,
            "valid_fraction": float(sub["valid_neighborhood"].mean()) if "valid_neighborhood" in sub.columns else np.nan,
        })

    return pd.DataFrame(records)


def plot_neighbor_count_diagnostics(comp_df: pd.DataFrame, output_dir: Path, sample_col: str = "sample"):
    output_dir.mkdir(parents=True, exist_ok=True)
    counts = comp_df["neighbor_count"].astype(float)

    plt.figure(figsize=(7, 5))
    plt.hist(counts, bins=100)
    plt.xlabel("Neighbor count within radius")
    plt.ylabel("Cell count")
    plt.title("Overall neighbor count distribution")
    plt.tight_layout()
    plt.savefig(output_dir / "neighbor_count_hist_overall.png", dpi=180)
    plt.close()

    samples = sorted(comp_df[sample_col].astype(str).unique().tolist())
    data = [comp_df.loc[comp_df[sample_col].astype(str) == s, "neighbor_count"].to_numpy() for s in samples]

    plt.figure(figsize=(max(10, len(samples) * 0.35), 5))
    plt.boxplot(data, tick_labels=samples, showfliers=False)
    plt.xticks(rotation=90)
    plt.ylabel("Neighbor count within radius")
    plt.title("Neighbor count by sample")
    plt.tight_layout()
    plt.savefig(output_dir / "neighbor_count_boxplot_by_sample.png", dpi=180)
    plt.close()


def _subsample_for_score(X: np.ndarray, random_state: int, max_cells_for_k_selection: int) -> np.ndarray:
    n = X.shape[0]
    if n > max_cells_for_k_selection:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(n, size=max_cells_for_k_selection, replace=False)
        log(f"[K SCORE] using subsample {max_cells_for_k_selection}/{n} valid cells")
        return X[idx]
    return X


def score_k_values(
    X: np.ndarray,
    k_values: List[int],
    random_state: int = 0,
    max_cells_for_k_selection: int = 100000,
) -> pd.DataFrame:
    X_eval = _subsample_for_score(X, random_state, max_cells_for_k_selection)
    records = []

    for k in k_values:
        k = int(k)
        if k < 2 or k >= X_eval.shape[0]:
            records.append({"k": k, "calinski_harabasz_score": np.nan, "status": "invalid_k"})
            continue

        km = KMeans(n_clusters=k, n_init=20, random_state=random_state)
        labels = km.fit_predict(X_eval)

        if len(np.unique(labels)) < 2:
            score = np.nan
        else:
            score = calinski_harabasz_score(X_eval, labels)

        records.append({"k": k, "calinski_harabasz_score": score, "status": "ok"})
        log(f"[K SCORE] k={k}, CH={score:.4f}")

    return pd.DataFrame(records)


def choose_k_by_ch_score(
    X: np.ndarray,
    k_min: int,
    k_max: int,
    random_state: int = 0,
    max_cells_for_k_selection: int = 100000,
) -> Tuple[int, pd.DataFrame]:
    n = X.shape[0]
    if n < k_min + 1:
        raise ValueError(f"Too few valid cells for k selection: n={n}, k_min={k_min}")

    k_scores = score_k_values(
        X,
        k_values=list(range(k_min, k_max + 1)),
        random_state=random_state,
        max_cells_for_k_selection=max_cells_for_k_selection,
    )

    valid_scores = k_scores[np.isfinite(k_scores["calinski_harabasz_score"])]
    if valid_scores.empty:
        raise RuntimeError("Could not select k. All CH scores are invalid.")

    best_k = int(valid_scores.sort_values("calinski_harabasz_score", ascending=False).iloc[0]["k"])
    return best_k, k_scores


def _parse_k_values(k_values) -> Optional[List[int]]:
    if k_values is None:
        return None

    import re as _re

    if isinstance(k_values, str):
        raw = _re.split(r"[,\s]+", k_values.strip())
    else:
        raw = []
        for x in k_values:
            raw.extend(_re.split(r"[,\s]+", str(x).strip()))

    vals = sorted(set(int(x) for x in raw if str(x).strip() != ""))
    return vals if vals else None


def cluster_valid_neighborhoods(
    comp_df: pd.DataFrame,
    k: Optional[int] = None,
    k_values: Optional[List[int]] = None,
    primary_k: Optional[int] = None,
    k_min: int = 3,
    k_max: int = 15,
    scale: bool = False,
    random_state: int = 0,
    max_cells_for_k_selection: int = 100000,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Cluster valid cells.

    Modes:
    1. --k_values is provided:
        fit KMeans for every k
        write CN_k{k} columns, e.g. CN_k6, CN_k8, CN_k10
        copy one selected k result to CN; controlled by --primary_k
    2. --k is provided:
        fit one KMeans and write CN plus CN_k{k}
    3. neither is provided:
        choose best k by CH score from [k_min, k_max]
    """
    out = comp_df.copy()
    prop_cols = [c for c in out.columns if c.startswith("neigh_prop__")]
    if len(prop_cols) == 0:
        raise KeyError("No neighborhood proportion columns found.")

    valid = out["valid_neighborhood"].astype(bool).to_numpy()
    valid_df = out.loc[valid].copy()

    if valid_df.empty:
        raise RuntimeError("No valid cells available for clustering.")

    X = valid_df[prop_cols].to_numpy(dtype=np.float32)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    if scale:
        X_cluster = StandardScaler().fit_transform(X)
    else:
        X_cluster = X

    if k_values is not None and len(k_values) > 0:
        k_values = sorted(set(int(x) for x in k_values))
        k_scores = score_k_values(
            X_cluster,
            k_values=k_values,
            random_state=random_state,
            max_cells_for_k_selection=max_cells_for_k_selection,
        )
        if primary_k is None:
            primary_k = k_values[0]
        if int(primary_k) not in k_values:
            raise ValueError(f"primary_k={primary_k} is not included in k_values={k_values}")
        k_list_to_fit = k_values
        chosen_k = int(primary_k)

    elif k is not None:
        k_list_to_fit = [int(k)]
        k_scores = score_k_values(
            X_cluster,
            k_values=k_list_to_fit,
            random_state=random_state,
            max_cells_for_k_selection=max_cells_for_k_selection,
        )
        chosen_k = int(k)

    else:
        chosen_k, k_scores = choose_k_by_ch_score(
            X_cluster,
            k_min=k_min,
            k_max=k_max,
            random_state=random_state,
            max_cells_for_k_selection=max_cells_for_k_selection,
        )
        k_list_to_fit = [int(chosen_k)]

    out["CN"] = pd.Series([pd.NA] * len(out), dtype="object")
    out["CN_k"] = int(chosen_k)

    center_records = []

    for current_k in k_list_to_fit:
        current_k = int(current_k)
        if current_k < 2 or current_k >= X_cluster.shape[0]:
            log(f"[WARN] skip invalid k={current_k}")
            continue

        log(f"[CLUSTER] fitting KMeans k={current_k} on {X_cluster.shape[0]} valid cells")
        km = KMeans(n_clusters=current_k, n_init=50, random_state=random_state)
        labels = km.fit_predict(X_cluster)

        label_col = f"CN_k{current_k}"
        out[label_col] = pd.Series([pd.NA] * len(out), dtype="object")
        out.loc[valid, label_col] = pd.Series(
            [f"CN{int(i)}" for i in labels],
            index=out.index[valid],
        ).values

        if current_k == int(chosen_k):
            out.loc[valid, "CN"] = out.loc[valid, label_col].values

        centers = pd.DataFrame(km.cluster_centers_, columns=prop_cols)
        centers.insert(0, "CN", [f"CN{i}" for i in range(current_k)])
        centers.insert(0, "k", current_k)
        centers["label_column"] = label_col
        centers["note"] = "centers are in scaled feature space" if scale else "centers are phenotype proportions"

        size_df = out.loc[valid].groupby(label_col).size().reset_index(name="n_cells")
        size_df = size_df.rename(columns={label_col: "CN"})
        centers = centers.merge(size_df, on="CN", how="left")
        center_records.append(centers)

    center_df = pd.concat(center_records, axis=0, ignore_index=True) if center_records else pd.DataFrame()

    return out, k_scores, center_df


def plot_k_selection(k_scores: pd.DataFrame, output_dir: Path):
    if k_scores is None or k_scores.empty:
        return
    if "calinski_harabasz_score" not in k_scores.columns:
        return
    if k_scores["calinski_harabasz_score"].isna().all():
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(6, 4))
    plt.plot(k_scores["k"], k_scores["calinski_harabasz_score"], marker="o")
    plt.xlabel("k")
    plt.ylabel("Calinski-Harabasz score")
    plt.title("K selection for neighborhood clustering")
    plt.tight_layout()
    plt.savefig(output_dir / "k_selection_calinski_harabasz.png", dpi=180)
    plt.close()


def run_pipeline(
    input_table: Path,
    output_dir: Path,
    radius: float,
    min_neighbors: int,
    max_neighbors: Optional[int],
    min_neighbors_quantile: Optional[float],
    max_neighbors_quantile: Optional[float],
    k: Optional[int],
    k_values: Optional[List[int]],
    primary_k: Optional[int],
    k_min: int,
    k_max: int,
    scale: bool,
    sample_col: str,
    cellid_col: str,
    x_col: str,
    y_col: str,
    phenotype_col: str,
    drop_undefined: bool,
    chunk_size: int,
    random_state: int,
    max_cells_for_k_selection: int,
    prefer_output: str,
    workers: int,
):
    output_dir = Path(output_dir)
    plot_dir = output_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    log(f"[READ] {input_table}")
    df = read_table(input_table)
    df = prepare_input_table(
        df,
        sample_col=sample_col,
        cellid_col=cellid_col,
        x_col=x_col,
        y_col=y_col,
        phenotype_col=phenotype_col,
        drop_undefined=drop_undefined,
    )

    log(f"[INFO] cells after preprocessing: {len(df)}")
    log(f"[INFO] samples: {df[sample_col].nunique()}")
    log(f"[INFO] phenotypes: {df[phenotype_col].nunique()}")

    phenotype_order = sorted(df[phenotype_col].dropna().astype(str).unique().tolist())
    log(f"[INFO] phenotype_order: {phenotype_order}")

    comp_df, sample_manifest = compute_neighborhood_composition_by_sample_parallel(
        df=df,
        radius=radius,
        output_dir=output_dir,
        sample_col=sample_col,
        cellid_col=cellid_col,
        x_col=x_col,
        y_col=y_col,
        phenotype_col=phenotype_col,
        phenotype_order=phenotype_order,
        chunk_size=chunk_size,
        workers=workers,
        prefer_output=prefer_output,
    )

    comp_df, filter_params = add_valid_neighborhood_flag(
        comp_df,
        min_neighbors=min_neighbors,
        max_neighbors=max_neighbors,
        min_neighbors_quantile=min_neighbors_quantile,
        max_neighbors_quantile=max_neighbors_quantile,
    )

    summary_df = summarize_neighbor_counts(comp_df, sample_col=sample_col)
    plot_neighbor_count_diagnostics(comp_df, plot_dir, sample_col=sample_col)

    clustered_df, k_scores, cn_centers = cluster_valid_neighborhoods(
        comp_df,
        k=k,
        k_values=k_values,
        primary_k=primary_k,
        k_min=k_min,
        k_max=k_max,
        scale=scale,
        random_state=random_state,
        max_cells_for_k_selection=max_cells_for_k_selection,
    )
    plot_k_selection(k_scores, plot_dir)

    # Large cell-level table: keep as parquet/pickle.
    comp_path = write_table(clustered_df, output_dir / "neighborhood_composition", prefer=prefer_output)

    # Small statistic tables: write as CSV for easier inspection.
    summary_path = write_small_csv(summary_df, output_dir / "neighborhood_count_summary")
    k_path = write_small_csv(k_scores, output_dir / "k_selection_scores")
    center_path = write_small_csv(cn_centers, output_dir / "neighborhood_cluster_centers")

    params = pd.DataFrame([{
        "input_table": str(input_table),
        "radius": radius,
        "min_neighbors": min_neighbors,
        "max_neighbors": max_neighbors,
        "min_neighbors_quantile": min_neighbors_quantile,
        "max_neighbors_quantile": max_neighbors_quantile,
        "effective_min_neighbors": filter_params["effective_min_neighbors"],
        "effective_max_neighbors": filter_params["effective_max_neighbors"],
        "n_cells": filter_params["n_cells"],
        "n_valid_cells": filter_params["n_valid_cells"],
        "valid_fraction": filter_params["valid_fraction"],
        "k": int(clustered_df["CN_k"].dropna().iloc[0]) if "CN_k" in clustered_df.columns and clustered_df["CN_k"].notna().any() else np.nan,
        "k_values": ",".join(map(str, k_values)) if k_values is not None else "",
        "primary_k": primary_k if primary_k is not None else "",
        "scale": scale,
        "drop_undefined": drop_undefined,
        "workers": workers,
        "sample_col": sample_col,
        "cellid_col": cellid_col,
        "x_col": x_col,
        "y_col": y_col,
        "phenotype_col": phenotype_col,
    }])
    params_path = write_small_csv(params, output_dir / "run_params")

    log("[DONE]")
    log(f"[OUTPUT] neighborhood composition + CN: {comp_path}")
    log(f"[OUTPUT] neighbor count summary:      {summary_path}")
    log(f"[OUTPUT] k selection scores:          {k_path}")
    log(f"[OUTPUT] CN centers:                  {center_path}")
    log(f"[OUTPUT] run params:                  {params_path}")
    log(f"[OUTPUT] plots:                       {plot_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute sample-wise radius-neighborhood composition and cluster valid cells."
    )
    parser.add_argument("--input_table", required=True, type=str, help="pkl/parquet/csv table from phenotype_filtered_h5ad_no_csv_v2.py")
    parser.add_argument("--output_dir", required=True, type=str)
    parser.add_argument("--radius", default=300, type=float)

    parser.add_argument("--min_neighbors", default=10, type=int)
    parser.add_argument("--max_neighbors", default=None, type=int)
    parser.add_argument("--min_neighbors_quantile", default=None, type=float, help="Optional lower quantile, e.g. 0.01")
    parser.add_argument("--max_neighbors_quantile", default=None, type=float, help="Optional upper quantile, e.g. 0.99")

    parser.add_argument("--k", default=None, type=int, help="Fixed number of CN clusters. If omitted, choose k by CH score.")
    parser.add_argument("--k_values", nargs="+", default=None, help="Multiple k values, e.g. --k_values 6 8 10 12 or --k_values 6,8,10,12.")
    parser.add_argument("--primary_k", default=None, type=int, help="When --k_values is used, this k is also copied to the generic CN column. Default: first k in k_values.")
    parser.add_argument("--k_min", default=3, type=int)
    parser.add_argument("--k_max", default=15, type=int)
    parser.add_argument("--scale", action="store_true", help="Scale neighborhood proportion columns before KMeans.")

    parser.add_argument("--sample_col", default="sample", type=str)
    parser.add_argument("--cellid_col", default="cellid", type=str)
    parser.add_argument("--x_col", default="spatial_x", type=str)
    parser.add_argument("--y_col", default="spatial_y", type=str)
    parser.add_argument("--phenotype_col", default="celltype", type=str)

    parser.add_argument("--drop_undefined", action="store_true", help="Remove undefined cells before neighborhood counting.")
    parser.add_argument("--chunk_size", default=20000, type=int)
    parser.add_argument("--workers", default=1, type=int, help="Number of samples processed in parallel for neighborhood counting.")
    parser.add_argument("--random_state", default=0, type=int)
    parser.add_argument("--max_cells_for_k_selection", default=100000, type=int)
    parser.add_argument("--prefer_output", default="parquet", choices=["parquet", "pickle"])

    args = parser.parse_args()

    run_pipeline(
        input_table=Path(args.input_table),
        output_dir=Path(args.output_dir),
        radius=args.radius,
        min_neighbors=args.min_neighbors,
        max_neighbors=args.max_neighbors,
        min_neighbors_quantile=args.min_neighbors_quantile,
        max_neighbors_quantile=args.max_neighbors_quantile,
        k=args.k,
        k_values=_parse_k_values(args.k_values),
        primary_k=args.primary_k,
        k_min=args.k_min,
        k_max=args.k_max,
        scale=args.scale,
        sample_col=args.sample_col,
        cellid_col=args.cellid_col,
        x_col=args.x_col,
        y_col=args.y_col,
        phenotype_col=args.phenotype_col,
        drop_undefined=args.drop_undefined,
        chunk_size=args.chunk_size,
        random_state=args.random_state,
        max_cells_for_k_selection=args.max_cells_for_k_selection,
        prefer_output=args.prefer_output,
        workers=args.workers,
    )


if __name__ == "__main__":
    main()
