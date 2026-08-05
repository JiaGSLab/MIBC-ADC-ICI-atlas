#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Apply final marker positivity calls to filtered h5ad files.

Rule
----
For each antibody / marker:
1. If Cohen's d >= cutoff:
      threshold = gmm_post75_threshold
      threshold_source = "gmm_post75"
2. If Cohen's d < cutoff:
      threshold = mean(elbow_chord_threshold, elbow_kneedle_threshold)
      threshold_source = "elbow_chord_kneedle_mean"

The script writes:
- adata.obs["positive_{antibody}"] as categorical calls:
    positive / negative / not_measured / no_threshold
- adata.layers["positive_call"]:
    1 = positive, 0 = negative, NaN = not measured or no threshold
- adata.var final threshold metadata
- {sample}_final_marker_thresholds.csv
- final_marker_calling_summary.csv

Example
-------
python3 apply_final_marker_positivity_calls.py \
  --input_root /data/home/chenwenhui/mIHC/CellInfo/BLCA_PixieResidual_noSOM \
  --cohen_d_cutoff 1.7 \
  --layer arcsinh_resid_mean \
  --only_markers \
  --overwrite
"""

from __future__ import annotations

import argparse
import json
import re
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd


def log(msg: str):
    print(msg, flush=True)


def sanitize_feature_name(x: str) -> str:
    x = re.sub(r"[^A-Za-z0-9_]+", "_", str(x))
    x = re.sub(r"_+", "_", x).strip("_")
    if re.match(r"^[0-9]", x):
        x = "m_" + x
    return x


def layer_token(layer: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(layer))


def json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        val = float(obj)
        return None if not np.isfinite(val) else val
    if isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    try:
        if pd.isna(obj):
            return None
    except Exception:
        pass
    return str(obj)


def find_samples(input_root: Path, samples: Optional[List[str]] = None) -> List[str]:
    if samples:
        return list(samples)
    return sorted([p.name for p in input_root.iterdir() if p.is_dir()])


def load_threshold_table(threshold_csv: Path, only_markers: bool = True) -> pd.DataFrame:
    if not threshold_csv.exists():
        raise FileNotFoundError(f"Threshold CSV not found: {threshold_csv}")

    df = pd.read_csv(threshold_csv)
    if "feature" not in df.columns:
        raise KeyError(f"'feature' column not found in {threshold_csv}")

    if only_markers and "feature_kind" in df.columns:
        df = df[df["feature_kind"].astype(str) == "marker"].copy()

    if "status" in df.columns:
        df = df[df["status"].astype(str) == "ok"].copy()

    numeric_cols = [
        "cohen_d",
        "gmm_post75_threshold",
        "gmm_post50_threshold",
        "gmm_post90_threshold",
        "gmm_intersection_threshold",
        "elbow_chord_threshold",
        "elbow_kneedle_threshold",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.drop_duplicates(subset=["feature"], keep="first").copy()
    df = df.set_index("feature", drop=False)
    return df


def finite_or_nan(x) -> float:
    try:
        val = float(x)
    except Exception:
        return np.nan
    return val if np.isfinite(val) else np.nan


def choose_threshold_for_feature(row: pd.Series, cohen_d_cutoff: float = 1.7, allow_fallback: bool = True) -> Dict:
    cohen_d = finite_or_nan(row.get("cohen_d", np.nan))
    gmm75 = finite_or_nan(row.get("gmm_post75_threshold", np.nan))
    elbow_chord = finite_or_nan(row.get("elbow_chord_threshold", np.nan))
    elbow_kneedle = finite_or_nan(row.get("elbow_kneedle_threshold", np.nan))

    elbow_vals = [x for x in [elbow_chord, elbow_kneedle] if np.isfinite(x)]
    elbow_mean = float(np.mean(elbow_vals)) if len(elbow_vals) > 0 else np.nan

    method = "threshold_unavailable"
    threshold = np.nan
    fallback_used = False

    if np.isfinite(cohen_d) and cohen_d >= cohen_d_cutoff:
        rule = f"cohen_d >= {cohen_d_cutoff}: use gmm_post75_threshold"
        if np.isfinite(gmm75):
            threshold = gmm75
            method = "gmm_post75"
        elif allow_fallback and np.isfinite(elbow_mean):
            threshold = elbow_mean
            method = "fallback_elbow_chord_kneedle_mean"
            fallback_used = True
        else:
            method = "gmm_post75_unavailable"

    elif np.isfinite(cohen_d) and cohen_d < cohen_d_cutoff:
        rule = f"cohen_d < {cohen_d_cutoff}: use mean(elbow_chord_threshold, elbow_kneedle_threshold)"
        if len(elbow_vals) == 2:
            threshold = elbow_mean
            method = "elbow_chord_kneedle_mean"
        elif allow_fallback and len(elbow_vals) == 1:
            threshold = elbow_vals[0]
            method = "fallback_single_elbow"
            fallback_used = True
        else:
            method = "elbow_unavailable"

    else:
        rule = "cohen_d unavailable"
        if allow_fallback and np.isfinite(elbow_mean):
            threshold = elbow_mean
            method = "fallback_elbow_cohen_d_unavailable"
            fallback_used = True
        elif allow_fallback and np.isfinite(gmm75):
            threshold = gmm75
            method = "fallback_gmm75_cohen_d_unavailable"
            fallback_used = True

    return {
        "cohen_d": cohen_d,
        "gmm_post75_threshold": gmm75,
        "elbow_chord_threshold": elbow_chord,
        "elbow_kneedle_threshold": elbow_kneedle,
        "final_threshold": threshold,
        "final_threshold_source": method,
        "final_threshold_rule": rule,
        "fallback_used": bool(fallback_used),
        "final_threshold_available": bool(np.isfinite(threshold)),
    }


def ensure_var_columns(adata: ad.AnnData):
    defaults = {
        "final_threshold": np.nan,
        "final_threshold_source": "",
        "final_threshold_rule": "",
        "final_threshold_available": False,
        "final_positive_count": 0,
        "final_negative_count": 0,
        "final_measured_count": 0,
        "final_positive_fraction": np.nan,
        "final_fallback_used": False,
    }
    for col, default in defaults.items():
        if col not in adata.var.columns:
            adata.var[col] = default


def apply_calls_to_adata(
    adata: ad.AnnData,
    threshold_df: pd.DataFrame,
    layer: str = "arcsinh_resid_mean",
    cohen_d_cutoff: float = 1.7,
    only_markers: bool = True,
    allow_fallback: bool = True,
    obs_prefix: str = "positive_",
) -> Tuple[ad.AnnData, pd.DataFrame]:
    if layer not in adata.layers:
        raise KeyError(f"Layer not found: {layer}")

    X = np.asarray(adata.layers[layer], dtype=np.float32)
    if X.shape != (adata.n_obs, adata.n_vars):
        raise ValueError(f"Layer {layer} has shape {X.shape}, expected {(adata.n_obs, adata.n_vars)}")

    if only_markers and "feature_kind" in adata.var.columns:
        feature_mask = adata.var["feature_kind"].astype(str).eq("marker").to_numpy()
    else:
        feature_mask = np.ones(adata.n_vars, dtype=bool)

    positive_layer = np.full((adata.n_obs, adata.n_vars), np.nan, dtype=np.float32)
    ensure_var_columns(adata)

    records = []
    categories = ["positive", "negative", "not_measured", "no_threshold"]

    for j, feature in enumerate(map(str, adata.var_names)):
        if not feature_mask[j]:
            continue

        safe_feature = sanitize_feature_name(feature)
        obs_col = f"{obs_prefix}{safe_feature}"

        if feature not in threshold_df.index:
            calls = np.full(adata.n_obs, "no_threshold", dtype=object)
            adata.obs[obs_col] = pd.Categorical(calls, categories=categories)
            adata.var.loc[feature, "final_threshold_available"] = False
            adata.var.loc[feature, "final_threshold_source"] = "missing_threshold_csv_row"
            adata.var.loc[feature, "final_threshold_rule"] = "feature missing from threshold CSV"
            records.append({
                "feature": feature,
                "obs_column": obs_col,
                "status": "missing_threshold_csv_row",
                "final_threshold_available": False,
                "final_threshold": np.nan,
                "final_threshold_source": "missing_threshold_csv_row",
                "final_threshold_rule": "feature missing from threshold CSV",
                "final_measured_count": 0,
                "final_positive_count": 0,
                "final_negative_count": 0,
                "final_positive_fraction": np.nan,
            })
            continue

        decision = choose_threshold_for_feature(
            threshold_df.loc[feature],
            cohen_d_cutoff=cohen_d_cutoff,
            allow_fallback=allow_fallback,
        )
        threshold = decision["final_threshold"]

        vals = X[:, j]
        measured = np.isfinite(vals)
        calls = np.full(adata.n_obs, "not_measured", dtype=object)

        if np.isfinite(threshold):
            pos = measured & (vals >= threshold)
            neg = measured & (vals < threshold)
            calls[pos] = "positive"
            calls[neg] = "negative"
            positive_layer[pos, j] = 1.0
            positive_layer[neg, j] = 0.0
        else:
            calls[measured] = "no_threshold"

        positive_count = int(np.sum(calls == "positive"))
        negative_count = int(np.sum(calls == "negative"))
        measured_count = int(np.sum(measured))
        denom = positive_count + negative_count
        positive_fraction = positive_count / denom if denom > 0 else np.nan

        adata.obs[obs_col] = pd.Categorical(calls, categories=categories)

        adata.var.loc[feature, "final_threshold"] = threshold
        adata.var.loc[feature, "final_threshold_source"] = decision["final_threshold_source"]
        adata.var.loc[feature, "final_threshold_rule"] = decision["final_threshold_rule"]
        adata.var.loc[feature, "final_threshold_available"] = decision["final_threshold_available"]
        adata.var.loc[feature, "final_positive_count"] = positive_count
        adata.var.loc[feature, "final_negative_count"] = negative_count
        adata.var.loc[feature, "final_measured_count"] = measured_count
        adata.var.loc[feature, "final_positive_fraction"] = positive_fraction
        adata.var.loc[feature, "final_fallback_used"] = decision["fallback_used"]

        rec = {
            "feature": feature,
            "obs_column": obs_col,
            "status": "ok" if decision["final_threshold_available"] else "no_threshold",
            "final_measured_count": measured_count,
            "final_positive_count": positive_count,
            "final_negative_count": negative_count,
            "final_positive_fraction": positive_fraction,
        }
        rec.update(decision)
        records.append(rec)

    token = layer_token(layer)
    positive_layer_name = f"positive_call_{token}"
    params_uns_key = f"final_marker_calling_params_{token}"
    thresholds_uns_key = f"final_marker_thresholds_json_{token}"
    adata.layers[positive_layer_name] = positive_layer
    final_df = pd.DataFrame(records)

    adata.uns[params_uns_key] = {
        "layer": layer,
        "cohen_d_cutoff": float(cohen_d_cutoff),
        "rule": (
            f"cohen_d >= {cohen_d_cutoff}: use gmm_post75_threshold; "
            f"cohen_d < {cohen_d_cutoff}: use mean(elbow_chord_threshold, elbow_kneedle_threshold)"
        ),
        "only_markers": bool(only_markers),
        "allow_fallback": bool(allow_fallback),
        "obs_prefix": obs_prefix,
        "positive_layer": positive_layer_name,
    }

    records_json = final_df.replace([np.inf, -np.inf], np.nan).where(pd.notnull(final_df), None).to_dict(orient="records")
    adata.uns[thresholds_uns_key] = json.dumps(records_json, default=json_default)

    return adata, final_df


def process_one_sample(
    sample: str,
    input_root: str,
    layer: str,
    cohen_d_cutoff: float,
    only_markers: bool,
    allow_fallback: bool,
    overwrite: bool,
    obs_prefix: str,
) -> Dict:
    input_root = Path(input_root)
    sample_dir = input_root / sample
    h5ad_path = sample_dir / f"{sample}_merged_panels_filtered.h5ad"
    token = layer_token(layer)
    threshold_csv = sample_dir / f"{sample}_{token}_thresholds.csv"
    output_csv = sample_dir / f"{sample}_final_marker_thresholds_{token}.csv"

    rec = {
        "sample": sample,
        "status": "",
        "h5ad_path": str(h5ad_path),
        "threshold_csv": str(threshold_csv),
        "output_csv": str(output_csv),
        "message": "",
    }

    if not h5ad_path.exists():
        rec.update(status="skip_missing_filtered_h5ad", message="filtered h5ad not found")
        return rec

    if not threshold_csv.exists():
        rec.update(status="skip_missing_threshold_csv", message="threshold CSV not found")
        return rec

    try:
        log(f"[READ] {sample}: {h5ad_path}")
        adata = ad.read_h5ad(h5ad_path)
        threshold_df = load_threshold_table(threshold_csv, only_markers=only_markers)

        log(f"[CALL] {sample}")
        adata, final_df = apply_calls_to_adata(
            adata=adata,
            threshold_df=threshold_df,
            layer=layer,
            cohen_d_cutoff=cohen_d_cutoff,
            only_markers=only_markers,
            allow_fallback=allow_fallback,
            obs_prefix=obs_prefix,
        )

        final_df.insert(0, "sample", sample)
        final_df.to_csv(output_csv, index=False)

        if overwrite:
            out_h5ad = h5ad_path
        else:
            out_h5ad = sample_dir / f"{sample}_merged_panels_filtered_with_calls_{token}.h5ad"

        log(f"[WRITE] {sample}: {out_h5ad}")
        adata.write_h5ad(out_h5ad)

        rec.update(
            status="ok",
            h5ad_path=str(out_h5ad),
            n_obs=int(adata.n_obs),
            n_vars=int(adata.n_vars),
            n_features_called=int(final_df.shape[0]),
        )

    except Exception as e:
        rec.update(status="error", message=str(e))
        log(f"[ERROR] {sample}: {e}")
        traceback.print_exc()

    return rec


def main():
    parser = argparse.ArgumentParser(description="Apply final marker positivity calls to filtered h5ad files.")
    parser.add_argument("--input_root", required=True, type=str)
    parser.add_argument("--samples", nargs="+", default=None)
    parser.add_argument("--layer", default="arcsinh_resid_mean", type=str)
    parser.add_argument("--cohen_d_cutoff", default=1.7, type=float)
    parser.add_argument("--only_markers", action="store_true")
    parser.add_argument("--no_fallback", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--obs_prefix", default="positive_", type=str)
    parser.add_argument("--workers", default=1, type=int)
    parser.add_argument("--summary_csv", default=None, type=str)
    args = parser.parse_args()

    input_root = Path(args.input_root)
    samples = find_samples(input_root, args.samples)

    log(f"[INFO] input_root: {input_root}")
    log(f"[INFO] samples: {len(samples)}")
    log(f"[INFO] layer: {args.layer}")
    log(f"[INFO] cohen_d_cutoff: {args.cohen_d_cutoff}")
    log(f"[INFO] only_markers: {args.only_markers}")
    log(f"[INFO] allow_fallback: {not args.no_fallback}")
    log(f"[INFO] overwrite: {args.overwrite}")
    log(f"[INFO] workers: {args.workers}")

    worker_kwargs = dict(
        input_root=str(input_root),
        layer=args.layer,
        cohen_d_cutoff=args.cohen_d_cutoff,
        only_markers=args.only_markers,
        allow_fallback=(not args.no_fallback),
        overwrite=args.overwrite,
        obs_prefix=args.obs_prefix,
    )

    results = []
    if args.workers <= 1:
        for sample in samples:
            results.append(process_one_sample(sample=sample, **worker_kwargs))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(process_one_sample, sample=sample, **worker_kwargs): sample
                for sample in samples
            }
            for future in as_completed(futures):
                sample = futures[future]
                try:
                    results.append(future.result())
                except Exception as e:
                    log(f"[ERROR] {sample}: worker failed unexpectedly: {e}")
                    results.append({"sample": sample, "status": "worker_error", "message": str(e)})

    summary = pd.DataFrame(results)
    summary_csv = Path(args.summary_csv) if args.summary_csv else input_root / "final_marker_calling_summary.csv"
    summary.to_csv(summary_csv, index=False)

    log("[DONE]")
    log(f"[SUMMARY] {summary_csv}")
    if "status" in summary.columns:
        log(summary["status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
