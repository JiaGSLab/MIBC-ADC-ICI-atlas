#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter merged h5ad by GeoJSON keep-region, compute threshold plots/CSV,
then cluster the filtered cells per panel using cell-level features.

Thresholding:
- GMM 2-component fit
- GMM intersection threshold
- GMM posterior thresholds: 0.50 / 0.75 / 0.90 for positive component
- Cohen's d annotation
- Smoothed elbow thresholds:
    1) chord-distance elbow
    2) Kneedle

Clustering:
- run separately for each panel in the filtered merged h5ad
- marker features only by default
- methods: kmeans, gmm, agglomerative
- output cluster labels are written back to obs
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import os
# Limit inner BLAS/OpenMP threads; sample-level workers handle outer parallelism.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("BLIS_NUM_THREADS", "1")
import re
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from scipy.signal import savgol_filter
from scipy.stats import norm
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler


def log(msg: str) -> None:
    print(msg, flush=True)


def _json_default(obj):
    """
    JSON serializer for numpy/pandas scalar objects.
    """
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        x = float(obj)
        return None if not np.isfinite(x) else x
    if isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    try:
        if pd.isna(obj):
            return None
    except Exception:
        pass
    return str(obj)


def _threshold_df_to_json_records(threshold_df: pd.DataFrame) -> str:
    """
    Store threshold records as a JSON string in .uns.

    This avoids h5ad write failures caused by storing a list of dictionaries
    with mixed strings, floats, numpy scalars and NaN values.
    """
    if threshold_df is None or threshold_df.empty:
        return "[]"
    df = threshold_df.copy()
    df = df.replace([np.inf, -np.inf], np.nan)
    records = df.where(pd.notnull(df), None).to_dict(orient="records")
    return json.dumps(records, default=_json_default)


def layer_token(layer: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(layer))


def thresholds_uns_key(layer: str) -> str:
    return f"{layer}_thresholds_json"


def sanitize_uns_for_h5ad(adata: ad.AnnData, layer: str = "arcsinh_resid_mean") -> ad.AnnData:
    """
    Remove or convert .uns entries known to cause h5ad write failures.
    """
    legacy_key = f"{layer}_thresholds"
    json_key = thresholds_uns_key(layer)
    if legacy_key in adata.uns:
        val = adata.uns[legacy_key]
        if not isinstance(val, (str, bytes)):
            try:
                adata.uns[json_key] = json.dumps(
                    val, default=_json_default
                )
            except Exception:
                adata.uns[json_key] = str(val)
            del adata.uns[legacy_key]
    return adata



# ---------------------------------------------------------------------
# GeoJSON filtering
# ---------------------------------------------------------------------
def _iter_geometries(obj):
    if isinstance(obj, list):
        for item in obj:
            yield from _iter_geometries(item)
        return
    if not isinstance(obj, dict):
        return

    obj_type = obj.get("type")
    if obj_type == "FeatureCollection":
        for feat in obj.get("features", []):
            yield from _iter_geometries(feat)
    elif obj_type == "Feature":
        geom = obj.get("geometry", None)
        if geom is not None:
            yield geom
    elif obj_type in {"Polygon", "MultiPolygon"}:
        yield obj
    elif "geometry" in obj and isinstance(obj["geometry"], dict):
        yield from _iter_geometries(obj["geometry"])


def _load_polygons_from_geojson(geojson_path: Path) -> List[List[np.ndarray]]:
    with open(geojson_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    polygons: List[List[np.ndarray]] = []
    for geom in _iter_geometries(data):
        gtype = geom.get("type")
        coords = geom.get("coordinates", [])

        if gtype == "Polygon":
            rings = [np.asarray(r, dtype=float) for r in coords if len(r) >= 3]
            if rings:
                polygons.append(rings)
        elif gtype == "MultiPolygon":
            for poly in coords:
                rings = [np.asarray(r, dtype=float) for r in poly if len(r) >= 3]
                if rings:
                    polygons.append(rings)

    if not polygons:
        raise ValueError(f"No Polygon/MultiPolygon geometry found in {geojson_path}")
    return polygons


def _points_in_polygon_with_holes(points: np.ndarray, rings: List[np.ndarray]) -> np.ndarray:
    outer = MplPath(rings[0], closed=True)
    inside = outer.contains_points(points)
    if len(rings) > 1:
        for hole in rings[1:]:
            hole_path = MplPath(hole, closed=True)
            inside &= ~hole_path.contains_points(points)
    return inside


def build_keep_mask(points: np.ndarray, polygons: List[List[np.ndarray]]) -> np.ndarray:
    keep = np.zeros(len(points), dtype=bool)
    for rings in polygons:
        keep |= _points_in_polygon_with_holes(points, rings)
    return keep


def filter_one_h5ad(
    input_h5ad: Path,
    geojson_path: Path,
    output_h5ad: Path,
    x_col: str = "centroidX",
    y_col: str = "centroidY",
    overwrite: bool = False,
) -> Tuple[ad.AnnData, Dict]:
    if output_h5ad.exists() and not overwrite:
        adata = ad.read_h5ad(output_h5ad)
        return adata, {
            "sample": input_h5ad.parent.name,
            "status": "skipped_exists",
            "input_h5ad": str(input_h5ad),
            "geojson": str(geojson_path),
            "output_h5ad": str(output_h5ad),
            "n_before": None,
            "n_after": int(adata.n_obs),
            "message": "output exists; reused existing filtered file",
        }

    adata = ad.read_h5ad(input_h5ad)
    if x_col not in adata.obs.columns or y_col not in adata.obs.columns:
        raise KeyError(f"obs must contain '{x_col}' and '{y_col}'")

    points = adata.obs[[x_col, y_col]].to_numpy(dtype=float)
    polygons = _load_polygons_from_geojson(geojson_path)
    keep = build_keep_mask(points, polygons)

    adata.obs["qc_keep_geojson"] = keep
    filtered = adata[keep].copy()

    filtered.uns["qc_geojson_file"] = str(geojson_path)
    filtered.uns["qc_geojson_filter_applied"] = True
    filtered.uns["qc_geojson_n_before"] = int(adata.n_obs)
    filtered.uns["qc_geojson_n_after"] = int(filtered.n_obs)
    filtered.uns["qc_geojson_x_col"] = x_col
    filtered.uns["qc_geojson_y_col"] = y_col
    filtered.uns["valid_region_source"] = "QC_geojson"

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    filtered.write_h5ad(output_h5ad)

    return filtered, {
        "sample": input_h5ad.parent.name,
        "status": "ok",
        "input_h5ad": str(input_h5ad),
        "geojson": str(geojson_path),
        "output_h5ad": str(output_h5ad),
        "n_before": int(adata.n_obs),
        "n_after": int(filtered.n_obs),
        "message": "",
        "valid_region_source": "QC_geojson",
    }


def filter_one_h5ad_by_shared_region(
    input_h5ad: Path,
    output_h5ad: Path,
    shared_region_col: str = "qc_three_panel_overlap",
    overwrite: bool = False,
) -> Tuple[ad.AnnData, Dict]:
    if output_h5ad.exists() and not overwrite:
        adata = ad.read_h5ad(output_h5ad)
        return adata, {
            "sample": input_h5ad.parent.name,
            "status": "skipped_exists",
            "input_h5ad": str(input_h5ad),
            "geojson": "",
            "output_h5ad": str(output_h5ad),
            "n_before": None,
            "n_after": int(adata.n_obs),
            "message": "output exists; reused existing filtered file",
            "valid_region_source": "shared_three_panel_cell_region",
        }

    adata = ad.read_h5ad(input_h5ad)
    if shared_region_col not in adata.obs.columns:
        raise KeyError(f"obs must contain '{shared_region_col}' for shared-region filtering")

    keep = adata.obs[shared_region_col].fillna(False).astype(bool).to_numpy()
    filtered = adata[keep].copy()

    filtered.uns["valid_region_source"] = "shared_three_panel_cell_region"
    filtered.uns["shared_three_panel_region_col"] = shared_region_col
    filtered.uns["shared_three_panel_region_filter_applied"] = True
    filtered.uns["shared_three_panel_region_n_before"] = int(adata.n_obs)
    filtered.uns["shared_three_panel_region_n_after"] = int(filtered.n_obs)

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    filtered.write_h5ad(output_h5ad)

    return filtered, {
        "sample": input_h5ad.parent.name,
        "status": "ok",
        "input_h5ad": str(input_h5ad),
        "geojson": "",
        "output_h5ad": str(output_h5ad),
        "n_before": int(adata.n_obs),
        "n_after": int(filtered.n_obs),
        "message": "",
        "valid_region_source": "shared_three_panel_cell_region",
    }


# ---------------------------------------------------------------------
# Threshold methods with smoothing
# ---------------------------------------------------------------------
def _cohen_d(mu1: float, mu2: float, sd1: float, sd2: float) -> float:
    pooled = np.sqrt((sd1 ** 2 + sd2 ** 2) / 2.0)
    if pooled < 1e-12:
        return np.nan
    return float((mu2 - mu1) / pooled)


def _make_odd(n: int) -> int:
    return n if n % 2 == 1 else n - 1


def smooth_desc_curve(
    y_desc: np.ndarray,
    method: str = "savgol",
    window: int = 51,
    polyorder: int = 3,
    moving_avg_window: int = 31,
) -> np.ndarray:
    y = np.asarray(y_desc, dtype=float)
    if len(y) < 7:
        return y.copy()

    if method == "none":
        return y.copy()

    if method == "savgol":
        w = min(window, len(y) - (1 - len(y) % 2))
        w = _make_odd(max(w, polyorder + 2))
        if w <= polyorder:
            return y.copy()
        try:
            return savgol_filter(y, window_length=w, polyorder=min(polyorder, w - 2), mode="interp")
        except Exception:
            return y.copy()

    if method == "moving_average":
        w = min(moving_avg_window, len(y))
        w = max(w, 3)
        kernel = np.ones(w, dtype=float) / w
        return np.convolve(y, kernel, mode="same")

    return y.copy()


def _elbow_idx_chord_desc(y_desc: np.ndarray) -> int:
    y = np.asarray(y_desc, dtype=float)
    n = len(y)
    if n < 3:
        return 0
    x = np.arange(1, n + 1, dtype=float)

    start = np.array([x[0], y[0]], dtype=float)
    end = np.array([x[-1], y[-1]], dtype=float)
    line_vec = end - start
    normv = np.linalg.norm(line_vec)
    if normv == 0:
        return 0
    line_unit = line_vec / normv

    pts = np.column_stack((x, y))
    vecs = pts - start
    proj_len = np.dot(vecs, line_unit)
    proj_pts = start + proj_len[:, None] * line_unit

    d = np.linalg.norm(pts - proj_pts, axis=1)
    d[(pts[:, 1] >= proj_pts[:, 1]) | ~np.isfinite(d)] = 0.0
    return int(np.argmax(d))


def elbow_threshold_chord(
    values: np.ndarray,
    smooth_method: str = "savgol",
    smooth_window: int = 51,
    smooth_polyorder: int = 3,
    lower_pct: float = 2.0,
    upper_pct: float = 99.9,
) -> Tuple[Optional[float], Optional[int]]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 5:
        return None, None

    vals_desc = np.sort(vals)[::-1]
    lo, hi = np.percentile(vals_desc, [lower_pct, upper_pct])
    mask = (vals_desc >= lo) & (vals_desc <= hi)

    if mask.sum() >= 3:
        idxs = np.flatnonzero(mask)
        y_raw = vals_desc[mask]
        y_smooth = smooth_desc_curve(y_raw, method=smooth_method, window=smooth_window, polyorder=smooth_polyorder)
        local_idx = _elbow_idx_chord_desc(y_smooth)
        elbow_idx = int(idxs[local_idx])
    else:
        y_smooth = smooth_desc_curve(vals_desc, method=smooth_method, window=smooth_window, polyorder=smooth_polyorder)
        elbow_idx = _elbow_idx_chord_desc(y_smooth)

    return float(vals_desc[elbow_idx]), elbow_idx


def elbow_threshold_kneedle(
    values: np.ndarray,
    smooth_method: str = "savgol",
    smooth_window: int = 51,
    smooth_polyorder: int = 3,
) -> Tuple[Optional[float], Optional[int]]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 5:
        return None, None

    y_raw = np.sort(vals)[::-1]
    y = smooth_desc_curve(y_raw, method=smooth_method, window=smooth_window, polyorder=smooth_polyorder)
    n = len(y)
    x = np.linspace(0, 1, n)

    y_min, y_max = y.min(), y.max()
    if np.isclose(y_max, y_min):
        return None, None
    y_norm = (y - y_min) / (y_max - y_min)
    z = 1.0 - y_norm

    score = z - x
    idx = int(np.argmax(score))
    return float(y_raw[idx]), idx


def _gmm_fit_2(values: np.ndarray, random_state: int = 0):
    """
    Fit a 2-component GMM. Components are ordered by mean:
    low-mean = negative-like, high-mean = positive-like.
    x_grid is restricted to the observed data range to avoid artificial
    negative thresholds from padded tails.
    """
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 20 or np.nanstd(vals) < 1e-8 or len(np.unique(vals)) < 5:
        return None, None, None, None

    gm = GaussianMixture(
        n_components=2,
        covariance_type="full",
        random_state=random_state,
        max_iter=1000,
        reg_covar=1e-6,
    )
    gm.fit(vals.reshape(-1, 1))

    raw_means = gm.means_.ravel()
    order = np.argsort(raw_means)
    means = raw_means[order]
    weights = gm.weights_.ravel()[order]

    cov = gm.covariances_
    if cov.ndim == 3:
        vars_ = cov[:, 0, 0][order]
    elif cov.ndim == 2:
        vars_ = cov[:, 0][order]
    else:
        vars_ = cov[order]
    sigmas = np.sqrt(np.maximum(vars_, 1e-12))

    x_min, x_max = np.nanmin(vals), np.nanmax(vals)
    if np.isclose(x_min, x_max):
        return None, None, None, None

    x_grid = np.linspace(x_min, x_max, 4000)
    pdf_neg = weights[0] * norm.pdf(x_grid, means[0], sigmas[0])
    pdf_pos = weights[1] * norm.pdf(x_grid, means[1], sigmas[1])
    post_pos = pdf_pos / np.clip(pdf_neg + pdf_pos, 1e-12, None)

    info = {
        "weights": weights.tolist(),
        "means": means.tolist(),
        "sigmas": sigmas.tolist(),
        "cohen_d": _cohen_d(float(means[0]), float(means[1]), float(sigmas[0]), float(sigmas[1])),
        "observed_min": float(x_min),
        "observed_max": float(x_max),
        "post_pos_at_mean_neg": float(np.interp(means[0], x_grid, post_pos)),
        "post_pos_at_mean_pos": float(np.interp(means[1], x_grid, post_pos)),
    }
    return info, x_grid, pdf_neg, pdf_pos


def _posterior_threshold(
    x_grid: np.ndarray,
    pdf_neg: np.ndarray,
    pdf_pos: np.ndarray,
    q: float,
    mean_neg: Optional[float] = None,
    mean_pos: Optional[float] = None,
    observed_min: Optional[float] = None,
    observed_max: Optional[float] = None,
) -> Optional[float]:
    """
    Return the first upward crossing where P(high-mean component | x) reaches q.

    This is the positivity threshold. It is constrained to:
    - the observed data range;
    - x >= low-component mean;
    - a real upward crossing from < q to >= q.

    If posterior is already >= q at the lower bound, return None because
    there is no meaningful threshold.
    """
    x_grid = np.asarray(x_grid, dtype=float)
    post = pdf_pos / np.clip(pdf_neg + pdf_pos, 1e-12, None)

    lo = observed_min if observed_min is not None else float(np.nanmin(x_grid))
    hi = observed_max if observed_max is not None else float(np.nanmax(x_grid))
    if mean_neg is not None and np.isfinite(mean_neg):
        lo = max(lo, float(mean_neg))

    valid = np.isfinite(post) & np.isfinite(x_grid) & (x_grid >= lo) & (x_grid <= hi)
    idx = np.where(valid)[0]
    if len(idx) < 2:
        return None

    xv = x_grid[idx]
    pv = post[idx]
    if pv[0] >= q:
        return None

    cross = np.where((pv[:-1] < q) & (pv[1:] >= q))[0]
    if len(cross) == 0:
        return None

    i = int(cross[0])
    x0, x1 = xv[i], xv[i + 1]
    p0, p1 = pv[i], pv[i + 1]
    if np.isclose(p1, p0):
        return float(x1)
    t = (q - p0) / (p1 - p0)
    return float(x0 + t * (x1 - x0))


def _intersection_threshold(
    x_grid: np.ndarray,
    pdf_neg: np.ndarray,
    pdf_pos: np.ndarray,
    mean_neg: Optional[float] = None,
    mean_pos: Optional[float] = None,
    observed_min: Optional[float] = None,
    observed_max: Optional[float] = None,
) -> Tuple[Optional[float], int]:
    """
    Weighted Gaussian-curve intersection.

    Only intersections between the two component means are considered valid.
    If several crossings exist, choose the one closest to the midpoint between
    means. If no valid crossing exists, return None and the total crossing count.
    """
    x_grid = np.asarray(x_grid, dtype=float)
    diff = pdf_neg - pdf_pos
    sign_change = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]
    if len(sign_change) == 0:
        return None, 0

    crossings = []
    for i in sign_change:
        x0, x1 = x_grid[i], x_grid[i + 1]
        d0, d1 = diff[i], diff[i + 1]
        if np.isclose(d1, d0):
            xc = x1
        else:
            t = -d0 / (d1 - d0)
            xc = x0 + t * (x1 - x0)
        crossings.append(float(xc))

    lo = observed_min if observed_min is not None else float(np.nanmin(x_grid))
    hi = observed_max if observed_max is not None else float(np.nanmax(x_grid))
    if mean_neg is not None and mean_pos is not None and np.isfinite(mean_neg) and np.isfinite(mean_pos):
        lo = max(lo, min(float(mean_neg), float(mean_pos)))
        hi = min(hi, max(float(mean_neg), float(mean_pos)))
        mid = 0.5 * (float(mean_neg) + float(mean_pos))
    else:
        mid = 0.5 * (lo + hi)

    valid_crossings = [c for c in crossings if lo <= c <= hi]
    if len(valid_crossings) == 0:
        return None, len(crossings)

    chosen = min(valid_crossings, key=lambda c: abs(c - mid))
    return float(chosen), len(crossings)


def _gmm_good_fit_label(
    cohen_d: float,
    weight_neg: Optional[float] = None,
    weight_pos: Optional[float] = None,
    post50_threshold: Optional[float] = None,
    intersection_threshold: Optional[float] = None,
    mean_neg: Optional[float] = None,
    mean_pos: Optional[float] = None,
) -> str:
    """
    Practical GMM quality label for thresholding.

    Cohen's d measures separation, but threshold validity and component weights
    are also required for a useful positivity threshold.
    """
    if not np.isfinite(cohen_d):
        return "undetermined"

    extreme_weight = False
    if weight_neg is not None and weight_pos is not None:
        if min(weight_neg, weight_pos) < 0.01 or max(weight_neg, weight_pos) > 0.99:
            extreme_weight = True

    has_valid_threshold = False
    for t in [post50_threshold, intersection_threshold]:
        if t is None or not np.isfinite(t):
            continue
        if mean_neg is not None and mean_pos is not None and np.isfinite(mean_neg) and np.isfinite(mean_pos):
            if min(mean_neg, mean_pos) <= t <= max(mean_neg, mean_pos):
                has_valid_threshold = True
        else:
            has_valid_threshold = True

    if cohen_d >= 2.0 and (not extreme_weight) and has_valid_threshold:
        return "good"
    if cohen_d >= 1.0 and has_valid_threshold:
        return "moderate"
    return "poor"


def plot_thresholds_for_filtered_adata(
    adata: ad.AnnData,
    sample: str,
    out_dir: Path,
    layer: str = "arcsinh_resid_mean",
    include_feature_kinds: Sequence[str] = ("marker",),
    smooth_method: str = "savgol",
    smooth_window: int = 51,
    smooth_polyorder: int = 3,
) -> pd.DataFrame:
    if layer not in adata.layers:
        raise KeyError(f"Layer '{layer}' not found in h5ad")

    out_dir.mkdir(parents=True, exist_ok=True)
    X = np.asarray(adata.layers[layer], dtype=np.float32)
    features = list(map(str, adata.var_names))
    feature_kind = (
        adata.var["feature_kind"].astype(str).tolist()
        if "feature_kind" in adata.var.columns
        else ["unknown"] * len(features)
    )

    records = []
    for j, feat in enumerate(features):
        if feature_kind[j] not in include_feature_kinds:
            continue

        log(f"    [threshold] {sample} {feat}")
        vals = X[:, j]
        vals = vals[np.isfinite(vals)]
        rec = {
            "sample": sample,
            "feature": feat,
            "feature_kind": feature_kind[j],
            "n_cells": int(len(vals)),
            "smooth_method": smooth_method,
            "smooth_window": int(smooth_window),
            "smooth_polyorder": int(smooth_polyorder),
            "gmm_weight_neg": np.nan,
            "gmm_weight_pos": np.nan,
            "gmm_mean_neg": np.nan,
            "gmm_mean_pos": np.nan,
            "gmm_sigma_neg": np.nan,
            "gmm_sigma_pos": np.nan,
            "cohen_d": np.nan,
            "gmm_fit_quality": "undetermined",
            "gmm_intersection_threshold": np.nan,
            "gmm_n_intersections": np.nan,
            "gmm_post50_threshold": np.nan,
            "gmm_post75_threshold": np.nan,
            "gmm_post90_threshold": np.nan,
            "elbow_chord_threshold": np.nan,
            "elbow_kneedle_threshold": np.nan,
            "suggested_method": "",
            "suggested_threshold": np.nan,
            "status": "ok",
        }

        if len(vals) < 10:
            rec["status"] = "too_few_values"
            records.append(rec)
            continue

        gmm_info, x_grid, pdf_neg, pdf_pos = _gmm_fit_2(vals, random_state=0)
        thr_chord, _ = elbow_threshold_chord(
            vals,
            smooth_method=smooth_method,
            smooth_window=smooth_window,
            smooth_polyorder=smooth_polyorder,
        )
        thr_kneedle, _ = elbow_threshold_kneedle(
            vals,
            smooth_method=smooth_method,
            smooth_window=smooth_window,
            smooth_polyorder=smooth_polyorder,
        )
        rec["elbow_chord_threshold"] = thr_chord if thr_chord is not None else np.nan
        rec["elbow_kneedle_threshold"] = thr_kneedle if thr_kneedle is not None else np.nan

        fig, ax = plt.subplots(figsize=(7.6, 5.6))
        ax.hist(vals, bins=100, density=True, alpha=0.40, label="Histogram")
        annotation_lines = [f"n={len(vals):,}", f"smooth={smooth_method}"]

        if gmm_info is not None and x_grid is not None:
            w_neg, w_pos = gmm_info["weights"]
            mu_neg, mu_pos = gmm_info["means"]
            sd_neg, sd_pos = gmm_info["sigmas"]
            cohen_d = gmm_info["cohen_d"]

            rec["gmm_weight_neg"] = w_neg
            rec["gmm_weight_pos"] = w_pos
            rec["gmm_mean_neg"] = mu_neg
            rec["gmm_mean_pos"] = mu_pos
            rec["gmm_sigma_neg"] = sd_neg
            rec["gmm_sigma_pos"] = sd_pos
            rec["cohen_d"] = cohen_d
            # Temporary label; updated below after threshold validity is known.
            rec["gmm_fit_quality"] = _gmm_good_fit_label(cohen_d)

            ax.plot(x_grid, pdf_neg, lw=2.0, label=f"GMM neg (w={w_neg:.2f}, μ={mu_neg:.2f}, σ={sd_neg:.2f})")
            ax.plot(x_grid, pdf_pos, lw=2.0, label=f"GMM pos (w={w_pos:.2f}, μ={mu_pos:.2f}, σ={sd_pos:.2f})")
            ax.plot(x_grid, pdf_neg + pdf_pos, lw=2.0, alpha=0.85, label="GMM mixture")

            observed_min = gmm_info.get("observed_min", np.nan)
            observed_max = gmm_info.get("observed_max", np.nan)

            thr_inter, n_intersections = _intersection_threshold(
                x_grid, pdf_neg, pdf_pos,
                mean_neg=mu_neg, mean_pos=mu_pos,
                observed_min=observed_min, observed_max=observed_max,
            )
            thr50 = _posterior_threshold(
                x_grid, pdf_neg, pdf_pos, 0.50,
                mean_neg=mu_neg, mean_pos=mu_pos,
                observed_min=observed_min, observed_max=observed_max,
            )
            thr75 = _posterior_threshold(
                x_grid, pdf_neg, pdf_pos, 0.75,
                mean_neg=mu_neg, mean_pos=mu_pos,
                observed_min=observed_min, observed_max=observed_max,
            )
            thr90 = _posterior_threshold(
                x_grid, pdf_neg, pdf_pos, 0.90,
                mean_neg=mu_neg, mean_pos=mu_pos,
                observed_min=observed_min, observed_max=observed_max,
            )

            rec["gmm_intersection_threshold"] = thr_inter if thr_inter is not None else np.nan
            rec["gmm_n_intersections"] = n_intersections
            rec["gmm_post50_threshold"] = thr50 if thr50 is not None else np.nan
            rec["gmm_post75_threshold"] = thr75 if thr75 is not None else np.nan
            rec["gmm_post90_threshold"] = thr90 if thr90 is not None else np.nan

            rec["gmm_fit_quality"] = _gmm_good_fit_label(
                cohen_d,
                weight_neg=w_neg,
                weight_pos=w_pos,
                post50_threshold=thr50,
                intersection_threshold=thr_inter,
                mean_neg=mu_neg,
                mean_pos=mu_pos,
            )

            if thr_inter is not None:
                ax.axvline(thr_inter, linestyle="--", linewidth=1.5, label=f"GMM intersect={thr_inter:.3f}")
            if thr50 is not None:
                ax.axvline(thr50, linestyle="--", linewidth=1.5, label=f"GMM post50={thr50:.3f}")
            if thr75 is not None:
                ax.axvline(thr75, linestyle="--", linewidth=1.5, label=f"GMM post75={thr75:.3f}")
            if thr90 is not None:
                ax.axvline(thr90, linestyle="--", linewidth=1.5, label=f"GMM post90={thr90:.3f}")

            annotation_lines.append(f"Cohen's d={cohen_d:.2f}")
            annotation_lines.append(f"GMM={rec['gmm_fit_quality']}")

        if thr_chord is not None:
            ax.axvline(thr_chord, linestyle=":", linewidth=1.4, label=f"Elbow chord={thr_chord:.3f}")
        if thr_kneedle is not None:
            ax.axvline(thr_kneedle, linestyle=":", linewidth=1.4, label=f"Kneedle={thr_kneedle:.3f}")

        if np.isfinite(rec["cohen_d"]) and rec["cohen_d"] >= 2.0 and np.isfinite(rec["gmm_post50_threshold"]):
            rec["suggested_method"] = "gmm_post50"
            rec["suggested_threshold"] = rec["gmm_post50_threshold"]
        elif np.isfinite(rec["elbow_kneedle_threshold"]):
            rec["suggested_method"] = "elbow_kneedle"
            rec["suggested_threshold"] = rec["elbow_kneedle_threshold"]
        elif np.isfinite(rec["elbow_chord_threshold"]):
            rec["suggested_method"] = "elbow_chord"
            rec["suggested_threshold"] = rec["elbow_chord_threshold"]
        elif np.isfinite(rec["gmm_intersection_threshold"]):
            rec["suggested_method"] = "gmm_intersection"
            rec["suggested_threshold"] = rec["gmm_intersection_threshold"]

        ax.text(
            0.98, 0.98,
            "\n".join(annotation_lines),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.85),
        )
        ax.set_title(f"{sample} | {feat} | {layer}")
        ax.set_xlabel(layer)
        ax.set_ylabel("Density")
        ax.legend(fontsize=7, loc="best")
        plt.tight_layout()
        plt.savefig(out_dir / f"{feat}_{layer_token(layer)}_gmm_elbow.png", dpi=180, bbox_inches="tight")
        plt.close()

        records.append(rec)

    return pd.DataFrame(records)


def write_thresholds_to_h5ad(
    adata: ad.AnnData,
    threshold_df: pd.DataFrame,
    layer: str = "arcsinh_resid_mean",
) -> ad.AnnData:
    legacy_key = f"{layer}_thresholds"
    json_key = thresholds_uns_key(layer)
    if threshold_df.empty:
        adata.uns.pop(legacy_key, None)
        adata.uns[json_key] = "[]"
        return sanitize_uns_for_h5ad(adata, layer=layer)

    df = threshold_df.copy().set_index("feature", drop=False)
    cols_to_copy = [
        "gmm_weight_neg",
        "gmm_weight_pos",
        "gmm_mean_neg",
        "gmm_mean_pos",
        "gmm_sigma_neg",
        "gmm_sigma_pos",
        "cohen_d",
        "gmm_fit_quality",
        "gmm_intersection_threshold",
        "gmm_n_intersections",
        "gmm_post50_threshold",
        "gmm_post75_threshold",
        "gmm_post90_threshold",
        "elbow_chord_threshold",
        "elbow_kneedle_threshold",
        "suggested_method",
        "suggested_threshold",
    ]

    for col in cols_to_copy:
        is_obj = col.endswith("quality") or col.endswith("method")
        values = pd.Series(index=adata.var_names, dtype=object if is_obj else float)
        overlap = adata.var_names.intersection(df.index)
        values.loc[overlap] = df.loc[overlap, col]
        adata.var[col] = values.values

    # Store thresholds in .uns as JSON string, not as list-of-dicts.
    # The threshold table is also available in the per-sample CSV and adata.var columns.
    adata.uns.pop(legacy_key, None)
    adata.uns[json_key] = _threshold_df_to_json_records(threshold_df)
    return sanitize_uns_for_h5ad(adata, layer=layer)


# ---------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------
def sanitize_name(x: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(x))


def get_matrix_for_panel(
    adata: ad.AnnData,
    panel: str,
    layer: str = "arcsinh_resid_mean",
    include_feature_kinds: Sequence[str] = ("marker",),
    min_finite_frac: float = 0.95,
    min_variance: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    panel_mask = adata.obs["panel"].astype(str).to_numpy() == str(panel)
    if panel_mask.sum() == 0:
        return np.zeros((0, 0), dtype=np.float32), panel_mask, []

    X = np.asarray(adata.layers[layer], dtype=np.float32)[panel_mask]
    var = adata.var.copy()

    if "feature_kind" in var.columns:
        feat_mask = var["feature_kind"].astype(str).isin(include_feature_kinds).to_numpy()
    else:
        feat_mask = np.ones(adata.n_vars, dtype=bool)

    X = X[:, feat_mask]
    feat_names = list(var.index[feat_mask].astype(str))

    # remove mostly-NaN or constant features within this panel
    keep_cols = []
    for j in range(X.shape[1]):
        col = X[:, j]
        finite = np.isfinite(col)
        if finite.mean() < min_finite_frac:
            continue
        if np.nanvar(col) < min_variance:
            continue
        keep_cols.append(j)

    if len(keep_cols) == 0:
        return np.zeros((X.shape[0], 0), dtype=np.float32), panel_mask, []

    X = X[:, keep_cols]
    feat_names = [feat_names[j] for j in keep_cols]

    # simple median fill for the few remaining NaNs
    for j in range(X.shape[1]):
        col = X[:, j]
        finite = np.isfinite(col)
        if not finite.all():
            med = np.nanmedian(col[finite]) if finite.any() else 0.0
            col[~finite] = med
            X[:, j] = col

    return X.astype(np.float32), panel_mask, feat_names


def cluster_panel_matrix(
    X: np.ndarray,
    k_values: Sequence[int],
    methods: Sequence[str] = ("kmeans", "gmm", "agglomerative"),
    pca_n_components: int = 10,
    random_state: int = 0,
) -> Dict[str, np.ndarray]:
    """
    Cluster one panel-level matrix.

    Robustness notes:
    - KMeans usually works reliably.
    - GaussianMixture can fail when one component collapses or has singular covariance.
      In that case, we retry with increasing reg_covar and then diag covariance.
      If all GMM attempts fail, GMM for that k is skipped instead of failing the sample.
    - Agglomerative can be slow for large panels; it is optional.
    """
    out = {}
    if X.shape[0] < 3 or X.shape[1] < 1:
        return out

    scaler = StandardScaler(with_mean=True, with_std=True)
    Xs = scaler.fit_transform(X).astype(np.float64)

    n_pca = min(pca_n_components, Xs.shape[1], Xs.shape[0] - 1)
    if n_pca >= 2:
        pca = PCA(n_components=n_pca, random_state=random_state)
        Z = pca.fit_transform(Xs).astype(np.float64)
    else:
        Z = Xs.astype(np.float64)

    for k in k_values:
        if k < 2 or k >= X.shape[0]:
            continue

        if "kmeans" in methods:
            try:
                km = KMeans(n_clusters=k, random_state=random_state, n_init=20)
                out[f"kmeans_k{k}"] = km.fit_predict(Z).astype(int)
            except Exception as e:
                log(f"      [WARN] KMeans k={k} failed: {e}")

        if "gmm" in methods:
            gmm_done = False
            # Try full covariance first, then stronger regularization, then diag covariance.
            gmm_attempts = [
                ("full", 1e-6),
                ("full", 1e-5),
                ("full", 1e-4),
                ("full", 1e-3),
                ("diag", 1e-6),
                ("diag", 1e-5),
                ("diag", 1e-4),
                ("diag", 1e-3),
            ]
            last_err = None
            for covariance_type, reg_covar in gmm_attempts:
                try:
                    gm = GaussianMixture(
                        n_components=k,
                        covariance_type=covariance_type,
                        reg_covar=reg_covar,
                        random_state=random_state,
                        max_iter=1000,
                        n_init=3,
                        init_params="kmeans",
                    )
                    labels = gm.fit_predict(Z).astype(int)
                    key = f"gmm_k{k}_{covariance_type}_reg{reg_covar:g}"
                    out[key] = labels
                    if covariance_type != "full" or reg_covar != 1e-6:
                        log(
                            f"      [WARN] GMM k={k} succeeded with "
                            f"covariance_type={covariance_type}, reg_covar={reg_covar:g}"
                        )
                    gmm_done = True
                    break
                except Exception as e:
                    last_err = e
                    continue

            if not gmm_done:
                log(f"      [WARN] GMM k={k} skipped after all retries. Last error: {last_err}")

        if "agglomerative" in methods:
            try:
                agg = AgglomerativeClustering(n_clusters=k, linkage="ward")
                out[f"agglomerative_k{k}"] = agg.fit_predict(Z).astype(int)
            except Exception as e:
                log(f"      [WARN] Agglomerative k={k} failed: {e}")

    return out


def plot_cluster_heatmap(
    X: np.ndarray,
    feature_names: Sequence[str],
    labels: np.ndarray,
    output_path: Path,
    title: str,
):
    df = pd.DataFrame(X, columns=feature_names)
    df["cluster"] = labels.astype(int)
    med = df.groupby("cluster")[list(feature_names)].median()
    z = med.apply(lambda col: (col - col.mean()) / (col.std() + 1e-8), axis=0)

    fig, ax = plt.subplots(figsize=(max(8, len(feature_names) * 0.6), max(4, z.shape[0] * 0.4)))
    import seaborn as sns
    sns.heatmap(z, cmap="RdBu_r", center=0, linewidths=0.3, linecolor="grey", ax=ax)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_cluster_pca(
    X: np.ndarray,
    labels: np.ndarray,
    output_path: Path,
    title: str,
):
    if X.shape[0] < 3:
        return
    scaler = StandardScaler(with_mean=True, with_std=True)
    Xs = scaler.fit_transform(X)
    n_pca = min(2, Xs.shape[1], Xs.shape[0] - 1)
    if n_pca < 2:
        return
    pca = PCA(n_components=2, random_state=0)
    Z = pca.fit_transform(Xs)

    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(Z[:, 0], Z[:, 1], c=labels, s=4, alpha=0.7)
    ax.set_title(title)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close()


def run_panel_clustering(
    adata: ad.AnnData,
    sample: str,
    output_dir: Path,
    layer: str = "arcsinh_resid_mean",
    include_feature_kinds: Sequence[str] = ("marker",),
    methods: Sequence[str] = ("kmeans", "gmm", "agglomerative"),
    k_values: Sequence[int] = (6, 8, 10),
    pca_n_components: int = 10,
) -> Tuple[ad.AnnData, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    token = layer_token(layer)
    recs = []

    if "panel" not in adata.obs.columns:
        raise KeyError("obs['panel'] not found in filtered h5ad")

    panels = sorted(pd.Series(adata.obs["panel"].astype(str)).unique())

    for panel in panels:
        panel_safe = sanitize_name(panel)
        X, panel_mask, feature_names = get_matrix_for_panel(
            adata,
            panel=panel,
            layer=layer,
            include_feature_kinds=include_feature_kinds,
        )

        if X.shape[0] < 3 or X.shape[1] < 1:
            recs.append({
                "sample": sample,
                "panel": panel,
                "status": "skipped_too_few_features_or_cells",
                "n_cells": int(panel_mask.sum()),
                "n_features": int(X.shape[1]),
            })
            continue

        results = cluster_panel_matrix(
            X,
            k_values=k_values,
            methods=methods,
            pca_n_components=pca_n_components,
            random_state=0,
        )

        for name, labels in results.items():
            col = f"cluster_{panel_safe}_{token}_{name}"
            all_labels = pd.Series(index=adata.obs.index, dtype="object")
            all_labels.loc[adata.obs.index[panel_mask]] = labels.astype(str)
            adata.obs[col] = all_labels.values

            plot_cluster_heatmap(
                X,
                feature_names,
                labels,
                output_dir / f"{sample}_{panel_safe}_{token}_{name}_heatmap.png",
                title=f"{sample} | {panel} | {layer} | {name}",
            )
            plot_cluster_pca(
                X,
                labels,
                output_dir / f"{sample}_{panel_safe}_{token}_{name}_pca.png",
                title=f"{sample} | {panel} | {layer} | {name}",
            )

            cluster_sizes = pd.Series(labels).value_counts().sort_index().to_dict()
            recs.append({
                "sample": sample,
                "panel": panel,
                "status": "ok",
                "method": name,
                "n_cells": int(panel_mask.sum()),
                "n_features": int(X.shape[1]),
                "feature_names": ",".join(feature_names),
                "cluster_sizes": json.dumps(cluster_sizes),
            })

    return adata, pd.DataFrame(recs)



# ---------------------------------------------------------------------
# Parallel sample worker
# ---------------------------------------------------------------------
def count_panels_in_h5ad(input_h5ad: Path, panel_col: str = "panel") -> Tuple[int, List[str]]:
    """Light panel-count check before full processing."""
    adata = None
    try:
        adata = ad.read_h5ad(input_h5ad, backed="r")
        if panel_col in adata.obs.columns:
            panels = sorted(pd.Series(adata.obs[panel_col].astype(str)).dropna().unique().tolist())
        elif "image_ID" in adata.obs.columns:
            panels = sorted(pd.Series(adata.obs["image_ID"].astype(str)).dropna().unique().tolist())
        elif "panels_loaded" in adata.uns:
            panels = sorted([str(x) for x in adata.uns["panels_loaded"]])
        else:
            panels = []
        return len(panels), panels
    finally:
        try:
            if adata is not None and getattr(adata, "file", None) is not None:
                adata.file.close()
        except Exception:
            pass


def process_one_sample(
    sample: str,
    input_root: str,
    geojson_dir: str,
    x_col: str,
    y_col: str,
    overwrite: bool,
    include_feature_kinds: Sequence[str],
    smooth_method: str,
    smooth_window: int,
    smooth_polyorder: int,
    cluster_layer: str,
    cluster_methods: Sequence[str],
    cluster_k: Sequence[int],
    cluster_pca_n: int,
    required_panels: int,
    panel_col: str,
    shared_region_samples: Sequence[str],
    shared_region_col: str,
) -> Dict:
    input_root = Path(input_root)
    geojson_dir = Path(geojson_dir)
    sample_dir = input_root / sample

    input_h5ad = sample_dir / f"{sample}_merged_panels.h5ad"
    output_h5ad = sample_dir / f"{sample}_merged_panels_filtered.h5ad"
    geojson_path = geojson_dir / f"{sample}.geojson"
    token = layer_token(cluster_layer)
    plot_dir = sample_dir / f"filtered_distributions_{token}"
    threshold_csv = sample_dir / f"{sample}_{token}_thresholds.csv"
    cluster_dir = sample_dir / f"filtered_panel_clustering_{token}"
    cluster_csv = sample_dir / f"{sample}_{token}_panel_clustering_summary.csv"

    log(f"\n[SAMPLE] {sample}")

    base_rec = {
        "sample": sample,
        "input_h5ad": str(input_h5ad),
        "geojson": str(geojson_path),
        "output_h5ad": str(output_h5ad),
        "threshold_csv": str(threshold_csv),
        "cluster_csv": str(cluster_csv),
        "n_before": None,
        "n_after": None,
        "n_panels": None,
        "panels": "",
        "valid_region_source": "QC_geojson",
    }

    if not input_h5ad.exists():
        rec = dict(base_rec)
        rec.update(status="skip_missing_h5ad", message="merged h5ad not found")
        log(f"[SKIP] {sample}: merged h5ad not found")
        return rec

    use_shared_region = sample in set(shared_region_samples)
    if not geojson_path.exists() and not use_shared_region:
        rec = dict(base_rec)
        rec.update(status="skip_missing_geojson", message="geojson not found")
        log(f"[SKIP] {sample}: geojson not found")
        return rec

    n_panels, panels = count_panels_in_h5ad(input_h5ad, panel_col=panel_col)
    log(f"[CHECK] {sample}: n_panels={n_panels}, panels={panels}")
    if required_panels > 0 and n_panels != required_panels:
        rec = dict(base_rec)
        rec.update(
            status="skip_panel_count",
            n_panels=n_panels,
            panels=",".join(panels),
            message=f"expected {required_panels} panels, got {n_panels}",
        )
        log(f"[SKIP] {sample}: {rec['message']}")
        return rec

    try:
        log(f"[FILTER] {sample}")
        if use_shared_region:
            filtered, rec = filter_one_h5ad_by_shared_region(
                input_h5ad=input_h5ad,
                output_h5ad=output_h5ad,
                shared_region_col=shared_region_col,
                overwrite=overwrite,
            )
        else:
            filtered, rec = filter_one_h5ad(
                input_h5ad=input_h5ad,
                geojson_path=geojson_path,
                output_h5ad=output_h5ad,
                x_col=x_col,
                y_col=y_col,
                overwrite=overwrite,
            )

        if "panel" in filtered.obs.columns:
            filtered_panels = sorted(pd.Series(filtered.obs["panel"].astype(str)).dropna().unique().tolist())
        elif "image_ID" in filtered.obs.columns:
            filtered_panels = sorted(pd.Series(filtered.obs["image_ID"].astype(str)).dropna().unique().tolist())
        else:
            filtered_panels = []
        n_filtered_panels = len(filtered_panels)
        log(f"[CHECK FILTERED] {sample}: n_panels={n_filtered_panels}, panels={filtered_panels}")
        if required_panels > 0 and n_filtered_panels != required_panels:
            rec.update(
                status="skip_filtered_panel_count",
                threshold_csv=str(threshold_csv),
                cluster_csv=str(cluster_csv),
                n_panels=n_filtered_panels,
                panels=",".join(filtered_panels),
                message=f"filtered data expected {required_panels} panels, got {n_filtered_panels}",
            )
            log(f"[SKIP] {sample}: {rec['message']}")
            return rec

        log(f"[THRESHOLD] {sample}")
        thr_df = plot_thresholds_for_filtered_adata(
            filtered,
            sample=sample,
            out_dir=plot_dir,
            layer=cluster_layer,
            include_feature_kinds=include_feature_kinds,
            smooth_method=smooth_method,
            smooth_window=smooth_window,
            smooth_polyorder=smooth_polyorder,
        )
        thr_df.to_csv(threshold_csv, index=False)
        filtered = write_thresholds_to_h5ad(filtered, thr_df, layer=cluster_layer)

        log(f"[CLUSTER] {sample}")
        filtered, cluster_df = run_panel_clustering(
            filtered,
            sample=sample,
            output_dir=cluster_dir,
            layer=cluster_layer,
            include_feature_kinds=include_feature_kinds,
            methods=cluster_methods,
            k_values=cluster_k,
            pca_n_components=cluster_pca_n,
        )
        cluster_df.to_csv(cluster_csv, index=False)

        filtered.uns["panel_clustering_params"] = {
            "layer": cluster_layer,
            "include_feature_kinds": list(include_feature_kinds),
            "methods": list(cluster_methods),
            "k_values": list(cluster_k),
            "pca_n_components": int(cluster_pca_n),
        }
        filtered.uns["threshold_params"] = {
            "smooth_method": smooth_method,
            "smooth_window": int(smooth_window),
            "smooth_polyorder": int(smooth_polyorder),
            "elbow_methods": ["chord", "kneedle"],
            "piecewise_removed": True,
        }

        log(f"[WRITE] {sample}: {output_h5ad}")
        filtered = sanitize_uns_for_h5ad(filtered, layer=cluster_layer)
        filtered.write_h5ad(output_h5ad)

        rec["threshold_csv"] = str(threshold_csv)
        rec["cluster_csv"] = str(cluster_csv)
        rec["n_panels"] = n_filtered_panels
        rec["panels"] = ",".join(filtered_panels)
        if rec.get("status") == "ok":
            pct = 100.0 * rec["n_after"] / max(rec["n_before"], 1)
            log(f"[OK] {sample}: kept {rec['n_after']:,} / {rec['n_before']:,} cells ({pct:.2f}%)")
        else:
            log(f"[OK] {sample}: {rec.get('message', '')}")
        return rec

    except Exception as e:
        log(f"[ERROR] {sample}: {e}")
        traceback.print_exc()
        rec = dict(base_rec)
        rec.update(
            status="error",
            n_panels=n_panels,
            panels=",".join(panels),
            message=str(e),
        )
        return rec
    finally:
        gc.collect()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Parallel filter merged h5ad by GeoJSON, compute GMM/chord/Kneedle thresholds, and cluster filtered cells per panel."
    )
    parser.add_argument("--input_root", type=str, required=True)
    parser.add_argument("--geojson_dir", type=str, required=True)
    parser.add_argument("--samples", nargs="+", default=None)
    parser.add_argument("--workers", type=int, default=1, help="Number of samples to process in parallel.")
    parser.add_argument("--x_col", type=str, default="centroidX")
    parser.add_argument("--y_col", type=str, default="centroidY")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--summary_csv", type=str, default=None)

    parser.add_argument("--required_panels", type=int, default=3, help="Skip samples that do not have exactly this many panels. Use <=0 to disable.")
    parser.add_argument("--panel_col", type=str, default="panel", help="obs column used to count panels.")
    parser.add_argument("--shared_region_samples", nargs="+", default=[], help="Samples without QC GeoJSON that should be filtered using the shared three-panel region flag.")
    parser.add_argument("--shared_region_col", type=str, default="qc_three_panel_overlap", help="obs column used for no-QC shared-region filtering.")

    parser.add_argument("--include_feature_kinds", nargs="+", default=["marker"])
    parser.add_argument("--smooth_method", type=str, default="savgol", choices=["savgol", "moving_average", "none"])
    parser.add_argument("--smooth_window", type=int, default=51)
    parser.add_argument("--smooth_polyorder", type=int, default=3)

    parser.add_argument("--cluster_layer", type=str, default="arcsinh_resid_mean")
    parser.add_argument("--layer", type=str, default=None, help="Alias for --cluster_layer; used for thresholding and clustering.")
    parser.add_argument("--cluster_methods", nargs="+", default=["kmeans", "gmm"], choices=["kmeans", "gmm", "agglomerative"])
    parser.add_argument("--cluster_k", nargs="+", type=int, default=[6, 8, 10])
    parser.add_argument("--cluster_pca_n", type=int, default=10)

    args = parser.parse_args()
    analysis_layer = args.layer if args.layer else args.cluster_layer

    input_root = Path(args.input_root)
    if args.samples is not None:
        samples = list(args.samples)
    else:
        samples = sorted([p.name for p in input_root.iterdir() if p.is_dir()])

    log(f"Samples ({len(samples)}): {samples}")
    log(f"Workers: {args.workers}")
    log(f"Required panels: {args.required_panels}")
    log("Threshold elbow methods: chord + kneedle only; piecewise removed")
    log(f"Cluster methods: {args.cluster_methods}; k={args.cluster_k}; PCA n={args.cluster_pca_n}")
    log(f"Analysis layer: {analysis_layer}")

    worker_kwargs = dict(
        input_root=str(args.input_root),
        geojson_dir=str(args.geojson_dir),
        x_col=args.x_col,
        y_col=args.y_col,
        overwrite=args.overwrite,
        include_feature_kinds=args.include_feature_kinds,
        smooth_method=args.smooth_method,
        smooth_window=args.smooth_window,
        smooth_polyorder=args.smooth_polyorder,
        cluster_layer=analysis_layer,
        cluster_methods=args.cluster_methods,
        cluster_k=args.cluster_k,
        cluster_pca_n=args.cluster_pca_n,
        required_panels=args.required_panels,
        panel_col=args.panel_col,
        shared_region_samples=args.shared_region_samples,
        shared_region_col=args.shared_region_col,
    )

    results = []
    if args.workers <= 1:
        for sample in samples:
            results.append(process_one_sample(sample=sample, **worker_kwargs))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            future_to_sample = {
                executor.submit(process_one_sample, sample=sample, **worker_kwargs): sample
                for sample in samples
            }
            for future in as_completed(future_to_sample):
                sample = future_to_sample[future]
                try:
                    rec = future.result()
                except Exception as e:
                    log(f"[ERROR] {sample}: worker failed unexpectedly: {e}")
                    traceback.print_exc()
                    rec = {
                        "sample": sample,
                        "status": "worker_error",
                        "message": str(e),
                    }
                results.append(rec)

    df = pd.DataFrame(results)
    summary_csv = Path(args.summary_csv) if args.summary_csv else input_root / "geojson_filter_cluster_summary.csv"
    df.to_csv(summary_csv, index=False)

    log("\nDone.")
    log(f"Summary CSV: {summary_csv}")
    if len(df) > 0 and "status" in df.columns:
        log(df["status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
