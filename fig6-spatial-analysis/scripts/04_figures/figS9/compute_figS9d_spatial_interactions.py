#!/usr/bin/env python3
"""Compute B0 FigS9d by sample, pooling all three registered panels."""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree


INPUT = Path("/data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh/all_samples_cell_phenotypes_B0.parquet")
REMOTE_ROOT = Path("/data/home/chenwenhui/mIHC/FigS9_redraw_20260731_v2")
GROUPS = ("Pre_CR", "Pre_NR")
RADIUS = 100.0
PERMUTATIONS = 1000

_FRAME: pd.DataFrame | None = None
_INDICES: dict[str, np.ndarray] | None = None
_CATEGORIES: list[str] | None = None
_PERMUTATIONS = PERMUTATIONS
_BASE_SEED = 20260731


def directed_counts_from_pairs(codes: np.ndarray, pairs: np.ndarray, k: int) -> np.ndarray:
    """Expand each undirected radius edge into both Scimap directions."""
    if pairs.size == 0:
        return np.zeros((k, k), dtype=np.int64)
    a = codes[pairs[:, 0]].astype(np.int64, copy=False)
    b = codes[pairs[:, 1]].astype(np.int64, copy=False)
    undirected = np.bincount(a * k + b, minlength=k * k).reshape(k, k)
    return undirected + undirected.T


def marginal_null(
    row_margins: np.ndarray,
    column_margins: np.ndarray,
    permutations: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample each exact edge-permutation marginal with Hypergeom(N,Cj,Ri)."""
    row_margins = np.asarray(row_margins, dtype=np.int64)
    column_margins = np.asarray(column_margins, dtype=np.int64)
    k = len(row_margins)
    total = int(row_margins.sum())
    mean = np.full((k, k), np.nan, dtype=np.float64)
    std = np.full((k, k), np.nan, dtype=np.float64)
    if total == 0:
        return mean, std
    rng = np.random.default_rng(seed)
    for i in range(k):
        for j in range(k):
            draws = rng.hypergeometric(
                int(column_margins[j]),
                int(total - column_margins[j]),
                int(row_margins[i]),
                size=permutations,
            )
            mean[i, j] = draws.mean()
            std[i, j] = draws.std(ddof=1)
    return mean, std


def _sample_seed(sample: str) -> int:
    payload = f"{_BASE_SEED}\0{sample}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def _compute_sample(sample: str) -> tuple[pd.DataFrame, dict]:
    assert _FRAME is not None and _INDICES is not None and _CATEGORIES is not None
    started = time.perf_counter()
    sub = _FRAME.iloc[_INDICES[sample]]
    groups = sub["Group"].dropna().unique()
    if len(groups) != 1:
        raise ValueError(f"{sample}: expected one Group, got {groups.tolist()}")
    coords = sub[["spatial_x", "spatial_y"]].to_numpy(dtype=np.float64, copy=False)
    lookup = {name: i for i, name in enumerate(_CATEGORIES)}
    codes = np.fromiter((lookup[x] for x in sub["celltype"]), dtype=np.int16, count=len(sub))
    pairs = cKDTree(coords).query_pairs(r=RADIUS, output_type="ndarray")
    observed = directed_counts_from_pairs(codes, pairs, len(_CATEGORIES))
    null_mean, null_std = marginal_null(
        observed.sum(axis=1), observed.sum(axis=0), _PERMUTATIONS, _sample_seed(sample)
    )
    difference = observed - null_mean
    direction = np.where(np.isnan(difference) | (difference >= 0), 1.0, -1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        magnitude = observed / observed.max(axis=1, keepdims=True)
        z = difference / null_std
    z[np.isnan(z)] = 0.0
    pvalue = stats.norm.sf(np.abs(z)) * 2.0
    score = magnitude * direction
    k = len(_CATEGORIES)
    out = pd.DataFrame(
        {
            "sample": sample,
            "Group": str(groups[0]),
            "phenotype": np.repeat(_CATEGORIES, k),
            "neighbour_phenotype": np.tile(_CATEGORIES, k),
            "score": score.ravel(),
            "pvalue": pvalue.ravel(),
            "observed_count": observed.ravel(),
            "null_mean": null_mean.ravel(),
            "null_std": null_std.ravel(),
        }
    )
    metrics = {
        "sample": sample,
        "Group": str(groups[0]),
        "n_cells": int(len(sub)),
        "n_undirected_edges": int(len(pairs)),
        "n_directed_edges": int(observed.sum()),
        "seconds": time.perf_counter() - started,
    }
    return out, metrics


def summarize(per_sample: pd.DataFrame) -> pd.DataFrame:
    x = per_sample.copy()
    x["binary_score"] = np.where(x["score"] > 0, 1.0, -1.0)
    x.loc[x["score"].isna(), "binary_score"] = np.nan
    x["Group"] = pd.Categorical(x["Group"], GROUPS, ordered=True)
    out = (
        x.groupby(["Group", "phenotype", "neighbour_phenotype"], observed=True, sort=False)
        .agg(
            binary_mean=("binary_score", "mean"),
            score_mean=("score", "mean"),
            pvalue_mean=("pvalue", "mean"),
            n_samples=("sample", "nunique"),
        )
        .reset_index()
        .sort_values(["Group", "phenotype", "neighbour_phenotype"], kind="stable")
    )
    out["Group"] = out["Group"].astype(str)
    out["significant"] = out["pvalue_mean"] <= 0.05
    out["masked_binary_mean"] = out["binary_mean"].where(out["significant"])
    return out


def main() -> None:
    import pyarrow.parquet as pq

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output-root", type=Path, default=REMOTE_ROOT)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--permutations", type=int, default=PERMUTATIONS)
    parser.add_argument("--seed", type=int, default=20260731)
    args = parser.parse_args()
    t0 = time.perf_counter()
    print(f"READ {args.input}", flush=True)
    table = pq.read_table(
        args.input,
        columns=["sample", "spatial_x", "spatial_y", "celltype", "Group"],
        filters=[("Group", "in", list(GROUPS))],
    )
    frame = table.to_pandas()
    print(f"LOADED rows={len(frame):,} samples={frame['sample'].nunique()}", flush=True)
    categories = sorted(frame["celltype"].astype(str).unique().tolist())
    indices = {
        str(sample): np.asarray(idx, dtype=np.int64)
        for sample, idx in frame.groupby("sample", sort=False).indices.items()
    }
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = args.output_root / f"output_{stamp}"
    output.mkdir(parents=True, exist_ok=False)
    per_dir = output / "per_sample"
    per_dir.mkdir()

    global _FRAME, _INDICES, _CATEGORIES, _PERMUTATIONS, _BASE_SEED
    _FRAME, _INDICES, _CATEGORIES = frame, indices, categories
    _PERMUTATIONS, _BASE_SEED = args.permutations, args.seed
    if "fork" not in mp.get_all_start_methods():
        raise RuntimeError("fork start method required")
    workers = min(args.workers, len(indices))
    results, metrics = [], []
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("fork")) as pool:
        futures = {pool.submit(_compute_sample, sample): sample for sample in indices}
        for done, future in enumerate(as_completed(futures), 1):
            result, metric = future.result()
            result.to_csv(per_dir / f"{metric['sample']}.csv", index=False)
            results.append(result)
            metrics.append(metric)
            print(
                f"DONE {done}/{len(futures)} {metric['sample']} cells={metric['n_cells']:,} "
                f"edges={metric['n_directed_edges']:,} sec={metric['seconds']:.1f}",
                flush=True,
            )
    combined = pd.concat(results, ignore_index=True).sort_values(
        ["Group", "sample", "phenotype", "neighbour_phenotype"], kind="stable"
    )
    grouped = summarize(combined)
    combined.to_csv(output / "FigS9d_per_sample_interactions.csv", index=False)
    grouped.to_csv(output / "FigS9d_group_summary.csv", index=False)
    for group in GROUPS:
        grouped[grouped["Group"] == group].to_csv(output / f"FigS9d_{group}_summary.csv", index=False)
    metadata = {
        "input": str(args.input),
        "input_size": args.input.stat().st_size,
        "output": str(output),
        "radius_px": RADIUS,
        "permutations": args.permutations,
        "seed": args.seed,
        "workers": workers,
        "sample_key": "sample",
        "panel_split": False,
        "groups": list(GROUPS),
        "n_samples": len(indices),
        "n_rows": len(frame),
        "phenotypes": categories,
        "sample_metrics": sorted(metrics, key=lambda x: x["sample"]),
        "runtime_seconds": time.perf_counter() - t0,
        "null": "1000 exact hypergeometric marginal draws induced by global edge-neighbour label permutation",
        "scimap_compatibility": "directed counts; ddof=1 zscore p; sign(observed-null); focal-row-max magnitude",
    }
    (output / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"COMPLETE output={output} runtime_seconds={metadata['runtime_seconds']:.1f}", flush=True)


if __name__ == "__main__":
    main()
