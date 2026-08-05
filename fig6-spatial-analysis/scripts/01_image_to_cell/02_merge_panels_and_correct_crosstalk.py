#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Residual-only mIF pipeline (no SOM).

This script reuses the current pixie_som_residual_merge.py logic for:
- tile reading
- optional local background subtraction
- cell feature extraction
- pixel-level AF/donor residual correction
- cell-level residual summary
- panel-level h5ad writing
- sample-level merged h5ad writing

It removes:
- pixel SOM
- cell SOM
- SOM-based QC / heatmaps / dendrograms
- X_pixel_cluster_fractions obsm

Design goal:
- keep the same residual-corrected h5ad/layer structure
- keep speed and memory usage lower than the SOM version
"""

from __future__ import annotations

import argparse
import gc
import glob
import importlib.util
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd


def _load_base_module():
    here = Path(__file__).resolve().parent
    base_path = here / "01_extract_registered_cell_features.py"
    if not base_path.exists():
        raise FileNotFoundError(
            f"Required base script not found next to this script: {base_path}"
        )
    spec = importlib.util.spec_from_file_location("pixie_base_mod", str(base_path))
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


B = _load_base_module()


def _cleanup_memory(*objs):
    try:
        del objs
    except Exception:
        pass
    gc.collect()
    if getattr(B, "_GPU", False):
        try:
            B.cp.get_default_memory_pool().free_all_blocks()
            B.cp.get_default_pinned_memory_pool().free_all_blocks()
        except Exception:
            pass


def build_panel_anndata_nosom(
    cell_df: pd.DataFrame,
    resid_cell_df: Optional[pd.DataFrame],
    crosstalk_cell_df: Optional[pd.DataFrame],
    feature_names: List[str],
    cofactors: Dict[str, float],
    panel: str,
    sample: str,
    sample_meta: Dict[str, object],
    dye_map: Optional[Dict[str, str]] = None,
    correction_models: Optional[Dict[str, Dict[str, object]]] = None,
    crosstalk_cofactors: Optional[Dict[str, float]] = None,
    crosstalk_metadata: Optional[Dict[str, object]] = None,
    pipeline_params: Optional[Dict[str, object]] = None,
) -> ad.AnnData:
    raw_layers, arcsinh_layers, x_stat_map = B._build_stat_matrices(
        cell_df, feature_names, cofactors
    )

    extra_layer_dict = {}
    if resid_cell_df is not None:
        resid_raw_layers, resid_arc_layers, _ = B._build_stat_matrices(
            resid_cell_df, feature_names, cofactors
        )
        extra_layer_dict.update(
            {f"resid_{k.replace('raw_', '')}": v for k, v in resid_raw_layers.items()}
        )
        extra_layer_dict.update(
            {
                f"arcsinh_resid_{k.replace('arcsinh_', '')}": v
                for k, v in resid_arc_layers.items()
            }
        )
    if crosstalk_cell_df is not None:
        crosstalk_cofactors = crosstalk_cofactors or cofactors
        crosstalk_raw_layers, crosstalk_arc_layers, _ = B._build_stat_matrices(
            crosstalk_cell_df, feature_names, crosstalk_cofactors
        )
        extra_layer_dict.update(
            {
                f"crosstalk_{k.replace('raw_', '')}": v
                for k, v in crosstalk_raw_layers.items()
            }
        )
        extra_layer_dict.update(
            {
                f"arcsinh_crosstalk_{k.replace('arcsinh_', '')}": v
                for k, v in crosstalk_arc_layers.items()
            }
        )

    obs_cols = [
        "Cell_ID",
        "image_ID",
        "centroidX",
        "centroidY",
        "area_nuclear",
        "area_cell",
        "eccentricity",
        "AF_high",
        "qc_suspect",
        "qc_bad_pixel_dominant",
        "bad_pixel_frac",
        "global_label",
        "label",
    ]
    obs_cols = [c for c in obs_cols if c in cell_df.columns]
    obs = cell_df[obs_cols].copy()
    if "Cell_ID" not in obs.columns:
        obs["Cell_ID"] = cell_df.index.astype(str)
    obs["sample_ID"] = sample
    obs["panel"] = panel
    for k, v in sample_meta.items():
        obs[k] = v
    obs.index = obs["Cell_ID"].astype(str)

    var = pd.DataFrame(index=pd.Index(feature_names, name="feature"))
    var["feature_kind"] = [B._feature_kind(f) for f in feature_names]
    var["panel"] = panel
    var["arcsinh_cofactor"] = [
        cofactors.get(f, B.ARCSINH_COFACTOR_DEFAULT) for f in feature_names
    ]
    var["default_x_stat"] = [x_stat_map[f] for f in feature_names]
    if dye_map is None:
        dye_map = {}
    var["dye"] = [
        dye_map.get(f, f if f in ["AF", "DAPI"] else np.nan) for f in feature_names
    ]
    if correction_models is None:
        correction_models = {}
    var["residual_correction_applied"] = [
        bool(correction_models.get(f, {}).get("applied", False))
        if f not in ["AF", "DAPI"]
        else False
        for f in feature_names
    ]
    var["residual_selected_donors"] = [
        ",".join(correction_models.get(f, {}).get("selected_donors", []))
        if f not in ["AF", "DAPI"]
        else ""
        for f in feature_names
    ]
    var["corr_af_raw"] = [
        float(correction_models.get(f, {}).get("corr_af", np.nan))
        if f not in ["AF", "DAPI"]
        else np.nan
        for f in feature_names
    ]
    var["corr_af_resid"] = [
        float(correction_models.get(f, {}).get("corr_resid_af", np.nan))
        if f not in ["AF", "DAPI"]
        else np.nan
        for f in feature_names
    ]

    x_layer = extra_layer_dict.get(
        "arcsinh_resid_selected", arcsinh_layers["arcsinh_selected"]
    )
    adata = ad.AnnData(X=x_layer, obs=obs, var=var)
    for k, v in raw_layers.items():
        adata.layers[k] = v
    for k, v in arcsinh_layers.items():
        adata.layers[k] = v
    for k, v in extra_layer_dict.items():
        adata.layers[k] = v

    adata.uns["panel"] = panel
    adata.uns["sample_ID"] = sample
    adata.uns["sample_metadata"] = sample_meta
    adata.uns["correction_models"] = correction_models
    adata.uns["dye_map"] = dye_map
    adata.uns["pipeline_params"] = pipeline_params if pipeline_params is not None else {}
    adata.uns["crosstalk_correction_params"] = crosstalk_metadata if crosstalk_metadata is not None else {}
    adata.uns["crosstalk_channel_order"] = list((crosstalk_metadata or {}).get("crosstalk_channel_order", []))
    adata.uns["crosstalk_matrix_version"] = (crosstalk_metadata or {}).get("crosstalk_matrix_version", "not_applied")
    adata.uns["af_correction_source"] = (crosstalk_metadata or {}).get("af_correction_source", "not_applied")
    adata.uns["exposure_normalization_applied"] = bool((crosstalk_metadata or {}).get("exposure_normalization_applied", False))
    adata.uns["crosstalk_arcsinh_cofactors"] = (
        {k: float(v) for k, v in (crosstalk_cofactors or {}).items()}
        if crosstalk_cofactors is not None
        else {}
    )
    adata.uns["note"] = (
        "No-SOM residual pipeline. Pixel-level AF/donor residual correction was applied, "
        "then residual summaries were aggregated to cell level. "
        "If crosstalk metadata is present, an exposure-normalized AF+crosstalk corrected branch is stored in parallel layers. "
        "No pixel/cell SOM results are stored in this file."
    )
    return adata


def process_panel_nosom(
    fluo_image: str,
    af_image: str,
    raw_qptiff_path: Optional[str],
    crosstalk_config_path: Optional[str],
    panel: str,
    sample: str,
    sample_V: str,
    sample_mask_dir: str,
    output_dir: str,
    sample_meta: Dict[str, object],
    bg_method: Optional[str] = None,
    bg_kwargs: Optional[Dict] = None,
    pixels_per_cell_residual: int = 0,
    af_corr_threshold: float = 0.20,
    donor_corr_threshold: float = 0.15,
    max_fit_pixels: int = 2_000_000,
    background_fit_quantile: float = 0.10,
) -> Optional[ad.AnnData]:
    t0 = B.datetime.now()
    B.slide_io.kill_jvm()

    number_part = os.path.splitext(os.path.splitext(os.path.basename(fluo_image))[0])[0].split("-")[-1]
    reader_cls_fluo = B.slide_io.get_slide_reader(fluo_image)
    reader_fluo = reader_cls_fluo(fluo_image)
    reader_cls_af = B.slide_io.get_slide_reader(af_image)
    reader_af = reader_cls_af(af_image)
    ch_names = reader_fluo.metadata.channel_names
    w_fluo, h_fluo = reader_fluo.metadata.slide_dimensions[0]
    w_af, h_af = reader_af.metadata.slide_dimensions[0]
    if w_fluo != w_af or h_fluo != h_af:
        print(f"[SKIP] dimension mismatch: {fluo_image}")
        B.slide_io.kill_jvm()
        return None

    try:
        antibodies = B.ANTIBODY_MAPPING[sample_V][number_part][len(ch_names)]
    except KeyError:
        raise ValueError(
            f"Channel count mismatch for {fluo_image}: channels={len(ch_names)}"
        )

    dye_map = B.infer_antibody_dye_map(antibodies)
    donor_map = B.infer_donor_map(antibodies)
    print(
        f"  Panel={panel}, antibodies={antibodies}, dyes={dye_map}, size=({w_fluo}×{h_fluo})"
    )
    if bg_method:
        print(f"  Background subtraction: {bg_method} {bg_kwargs}")
    print(
        f"  pixels_per_cell_residual={pixels_per_cell_residual} (<=0 means all pixels)"
    )
    if raw_qptiff_path and crosstalk_config_path:
        print(f"  Crosstalk branch: enabled ({raw_qptiff_path})")
    else:
        print("  Crosstalk branch: disabled")

    fluo_vips = reader_fluo.slide2vips(level=0)
    af_vips = reader_af.slide2vips(level=0)

    af_full = B.pyvips_to_numpy(af_vips)
    af_small = B.cv2.resize(
        af_full, dsize=(0, 0), fx=0.2, fy=0.2, interpolation=B.cv2.INTER_AREA
    )
    del af_full
    af_thr = B.threshold_multiotsu(af_small, classes=B.AF_OTSU_CLASSES)[1]
    af_thr = float(np.clip(af_thr, B.AF_THRESHOLD_FLOOR, B.AF_THRESHOLD_CEIL))
    del af_small

    all_cell_dfs: List[pd.DataFrame] = []
    all_pixel_dfs: List[pd.DataFrame] = []
    label_offset = 0
    bg_tracker = (
        B.TileBackgroundTracker(["DAPI", "AF"] + list(antibodies)) if bg_method else None
    )

    for y in range(0, h_fluo, B.TILE_WH - B.OVERLAP):
        for x in range(0, w_fluo, B.TILE_WH - B.OVERLAP):
            cw = min(B.TILE_WH, w_fluo - x)
            ch = min(B.TILE_WH, h_fluo - y)
            if cw <= B.MIN_WH or ch <= B.MIN_WH:
                continue
            mask_file = os.path.join(sample_mask_dir, f"{panel}_mask_nuclear_{y}_{x}.tif")
            if not os.path.exists(mask_file):
                continue

            with B.Image.new_from_file(mask_file, access="sequential") as mask_vips:
                mask_np = B.pyvips_to_numpy(mask_vips)
            fluo_tile = B.pyvips_to_numpy(fluo_vips.extract_area(x, y, cw, ch))
            af_tile = B.pyvips_to_numpy(af_vips.extract_area(x, y, cw, ch))
            cell_df, pixel_df, tile_bg = B.cell_and_pixel_extraction(
                mask_np,
                fluo_tile,
                af_tile,
                antibodies,
                pixel_sample_k=pixels_per_cell_residual,
                bg_method=bg_method,
                bg_kwargs=bg_kwargs,
            )
            if bg_tracker is not None and tile_bg is not None:
                bg_tracker.add_tile(
                    tile_id=f"{panel}_{y}_{x}", x=x, y=y, backgrounds=tile_bg
                )
            if cell_df is not None and len(cell_df) > 0:
                max_lbl = int(cell_df["label"].max())
                cell_df["global_label"] = cell_df["label"] + label_offset
                cell_df["centroidX"] += x
                cell_df["centroidY"] += y
                if pixel_df is not None and len(pixel_df) > 0:
                    pixel_df["cell_label"] = pixel_df["cell_label"].astype(np.int64) + label_offset
                    pixel_df["pixel_x"] += x
                    pixel_df["pixel_y"] += y
                    all_pixel_dfs.append(pixel_df)
                all_cell_dfs.append(cell_df)
                label_offset += max_lbl + 1

            _cleanup_memory(mask_np, fluo_tile, af_tile, cell_df, pixel_df)

    _cleanup_memory(fluo_vips, af_vips)
    if not all_cell_dfs:
        print(f"  [SKIP] no cells for panel {panel}")
        B.slide_io.kill_jvm()
        return None

    panel_out = os.path.join(output_dir, sample)
    plot_dir = os.path.join(panel_out, f"plots_{panel}")
    B._safe_mkdir(plot_dir)

    cell_df_all = pd.concat(all_cell_dfs, ignore_index=True)
    pixel_df_all = pd.concat(all_pixel_dfs, ignore_index=True) if all_pixel_dfs else None
    _cleanup_memory(all_cell_dfs, all_pixel_dfs)

    if bg_tracker is not None and len(bg_tracker.records) > 0:
        bg_df = bg_tracker.to_dataframe()
        bg_df.to_csv(os.path.join(panel_out, f"tile_backgrounds_{panel}.csv"), index=False)
        bg_tracker.plot_backgrounds(os.path.join(plot_dir, "background_distributions.png"))
        bg_tracker.plot_background_heatmap(
            os.path.join(plot_dir, "background_spatial_heatmap.png"), w_fluo, h_fluo
        )
        bg_summary = bg_tracker.get_summary()
        _cleanup_memory(bg_df)
    else:
        bg_summary = None

    cell_df_all.insert(0, "image_ID", panel)
    cell_df_all, _ = B.remove_duplicates_by_distance(cell_df_all)
    kept_labels = set(cell_df_all["global_label"].values)
    if pixel_df_all is not None:
        pixel_df_all = pixel_df_all[pixel_df_all["cell_label"].isin(kept_labels)].reset_index(drop=True)

    cell_df_all["Cell_ID"] = (
        cell_df_all["image_ID"]
        + "_"
        + cell_df_all["centroidY"].round(0).astype(int).astype(str)
        + "_"
        + cell_df_all["centroidX"].round(0).astype(int).astype(str)
    )
    cell_df_all = cell_df_all.set_index("Cell_ID", drop=False)

    if bg_method and "AF_cell_p90" in cell_df_all.columns:
        af_dynamic_thr = float(np.percentile(cell_df_all["AF_cell_p90"].values, 95))
        cell_df_all["AF_high"] = cell_df_all["AF_cell_p90"] > af_dynamic_thr
    else:
        cell_df_all["AF_high"] = (
            cell_df_all.get("AF_cell_p90", pd.Series(False, index=cell_df_all.index)) > af_thr
        )

    cell_df_clean = cell_df_all.copy()
    pixel_df_clean = pixel_df_all
    if pixel_df_clean is None or len(pixel_df_clean) < 100:
        print(f"  [SKIP] too few sampled pixels for panel {panel}")
        B.slide_io.kill_jvm()
        return None

    ab_cols_pure = [c for c in antibodies if c in pixel_df_clean.columns]
    pixel_df_clean, cell_pixel_qc = B.pixel_level_qc(pixel_df_clean, ab_cols_pure)
    cell_df_clean = cell_df_clean.merge(
        cell_pixel_qc[["cell_label", "bad_pixel_frac", "qc_bad_pixel_dominant"]],
        left_on="global_label",
        right_on="cell_label",
        how="left",
    )
    cell_df_clean["bad_pixel_frac"] = cell_df_clean["bad_pixel_frac"].fillna(0)
    cell_df_clean["qc_bad_pixel_dominant"] = cell_df_clean["qc_bad_pixel_dominant"].fillna(False)
    if "cell_label" in cell_df_clean.columns:
        cell_df_clean = cell_df_clean.drop(columns=["cell_label"])

    pixel_df_corr, correction_models = B.fit_apply_pixel_residual_models(
        pixel_df_clean,
        antibodies=ab_cols_pure,
        donor_map=donor_map,
        af_corr_threshold=af_corr_threshold,
        donor_corr_threshold=donor_corr_threshold,
        max_fit_pixels=max_fit_pixels,
        background_fit_quantile=background_fit_quantile,
    )
    B.plot_raw_vs_residual_diagnostics(
        pixel_df_corr, ab_cols_pure, correction_models, plot_dir
    )

    resid_summary = B.summarise_residual_pixels_to_cells(pixel_df_corr, ab_cols_pure)
    resid_cell_df = cell_df_clean.copy()
    if not resid_summary.empty:
        resid_cell_df = resid_cell_df.merge(
            resid_summary,
            left_on="global_label",
            right_on="cell_label",
            how="left",
            suffixes=("", "_resid"),
        )
        if "cell_label" in resid_cell_df.columns:
            resid_cell_df = resid_cell_df.drop(columns=["cell_label"])
        for ab in ab_cols_pure:
            for stat in ["mean", "max", "min", "std", "p50", "p75", "p90"]:
                col = f"intensity_{stat}_{ab}"
                resid_col = f"{col}_resid"
                if resid_col in resid_cell_df.columns:
                    resid_cell_df[col] = resid_cell_df[resid_col].astype(np.float32)
                    resid_cell_df = resid_cell_df.drop(columns=[resid_col])

    cell_df_clean["qc_suspect"] = cell_df_clean["qc_bad_pixel_dominant"] | cell_df_clean["AF_high"]
    resid_cell_df["qc_bad_pixel_dominant"] = cell_df_clean["qc_bad_pixel_dominant"].values
    resid_cell_df["bad_pixel_frac"] = cell_df_clean["bad_pixel_frac"].values
    resid_cell_df["AF_high"] = cell_df_clean["AF_high"].values
    resid_cell_df["qc_suspect"] = cell_df_clean["qc_suspect"].values

    crosstalk_cell_df = None
    crosstalk_cofactors = None
    crosstalk_metadata: Dict[str, object] = {"status": "not_requested"}
    if raw_qptiff_path and crosstalk_config_path:
        try:
            correction_config = B.load_crosstalk_correction_config(crosstalk_config_path)
            pixel_df_crosstalk, crosstalk_metadata = B.apply_crosstalk_correction_to_pixel_df(
                pixel_df=pixel_df_clean,
                antibodies=ab_cols_pure,
                dye_map=dye_map,
                qptiff_path=raw_qptiff_path,
                correction_config=correction_config,
            )
            feature_to_column = {"AF": "crosstalk_AF", "DAPI": "crosstalk_DAPI"}
            feature_to_column.update({ab: f"crosstalk_{ab}" for ab in ab_cols_pure})
            crosstalk_summary = B.summarize_feature_columns_to_cells(
                pixel_df_crosstalk,
                feature_to_column=feature_to_column,
            )
            crosstalk_cell_df = cell_df_clean.copy()
            if not crosstalk_summary.empty:
                crosstalk_cell_df = crosstalk_cell_df.merge(
                    crosstalk_summary,
                    left_on="global_label",
                    right_on="cell_label",
                    how="left",
                    suffixes=("", "_crosstalk"),
                )
                if "cell_label" in crosstalk_cell_df.columns:
                    crosstalk_cell_df = crosstalk_cell_df.drop(columns=["cell_label"])
                for feature in ["AF", "DAPI"] + list(ab_cols_pure):
                    for stat in ["mean", "max", "min", "std", "p50", "p75", "p90"]:
                        col = B._feature_stat_col(feature, stat)
                        corrected_col = f"{col}_crosstalk"
                        if corrected_col in crosstalk_cell_df.columns:
                            crosstalk_cell_df[col] = crosstalk_cell_df[corrected_col].astype(np.float32)
                            crosstalk_cell_df = crosstalk_cell_df.drop(columns=[corrected_col])
            crosstalk_cell_df["qc_bad_pixel_dominant"] = cell_df_clean["qc_bad_pixel_dominant"].values
            crosstalk_cell_df["bad_pixel_frac"] = cell_df_clean["bad_pixel_frac"].values
            crosstalk_cell_df["AF_high"] = cell_df_clean["AF_high"].values
            crosstalk_cell_df["qc_suspect"] = cell_df_clean["qc_suspect"].values
            crosstalk_feature_df = pd.DataFrame(
                {
                    "AF": pixel_df_crosstalk["crosstalk_AF"].to_numpy(dtype=np.float32),
                    "DAPI": pixel_df_crosstalk["crosstalk_DAPI"].to_numpy(dtype=np.float32),
                }
            )
            for ab in ab_cols_pure:
                col = f"crosstalk_{ab}"
                if col in pixel_df_crosstalk.columns:
                    crosstalk_feature_df[ab] = pixel_df_crosstalk[col].to_numpy(dtype=np.float32)
            crosstalk_metadata["status"] = "ok"
            crosstalk_cofactors = B.estimate_cofactors(
                crosstalk_feature_df,
                ["AF", "DAPI"] + list(ab_cols_pure),
                already_bg_subtracted=True,
                min_cofactor=B.CROSSTALK_ARCSINH_MIN_COFACTOR,
                fallback_cofactor=B.CROSSTALK_ARCSINH_MIN_COFACTOR,
            )
            crosstalk_metadata["crosstalk_arcsinh_min_cofactor"] = float(
                B.CROSSTALK_ARCSINH_MIN_COFACTOR
            )
            crosstalk_metadata["crosstalk_cofactors"] = {
                k: float(v) for k, v in crosstalk_cofactors.items()
            }
        except Exception as e:
            crosstalk_cell_df = None
            crosstalk_cofactors = None
            crosstalk_metadata = {
                "status": "failed",
                "error": str(e),
                "qptiff_path": str(raw_qptiff_path),
                "config_path": str(crosstalk_config_path),
            }
            print(f"  [WARNING] crosstalk branch failed for {panel}: {e}")

    panel_feature_names = ["AF", "DAPI"] + list(ab_cols_pure)
    cofactors = B.estimate_cofactors(pixel_df_corr, panel_feature_names, already_bg_subtracted=True)
    raw_layers, arcsinh_layers, _ = B._build_stat_matrices(
        cell_df_clean, panel_feature_names, cofactors
    )
    resid_raw_layers, resid_arcsinh_layers, _ = B._build_stat_matrices(
        resid_cell_df, panel_feature_names, cofactors
    )
    try:
        B.plot_feature_distributions(
            raw_layers,
            arcsinh_layers,
            panel_feature_names,
            plot_dir,
            residual_layers=resid_raw_layers,
            residual_arcsinh_layers=resid_arcsinh_layers,
        )
    except Exception as e:
        print(f"  [WARNING] plot_feature_distributions failed for {panel}: {e}")

    B.plot_af_scatter(resid_cell_df, ab_cols_pure, plot_dir, cofactors)

    pipeline_params = {
        "bg_method": bg_method,
        "bg_kwargs": bg_kwargs,
        "bg_summary": bg_summary,
        "distribution_stats": list(B.DISTRIBUTION_STATS),
        "default_x_stat": B.DEFAULT_X_STAT,
        "dye_map": dye_map,
        "donor_map": donor_map,
        "correction_thresholds": {
            "af_corr_threshold": af_corr_threshold,
            "donor_corr_threshold": donor_corr_threshold,
            "max_fit_pixels": max_fit_pixels,
            "background_fit_quantile": background_fit_quantile,
        },
        "pixels_per_cell_residual": pixels_per_cell_residual,
        "crosstalk_branch_status": crosstalk_metadata.get("status", "not_requested"),
        "nosom": True,
    }

    adata = build_panel_anndata_nosom(
        cell_df_clean,
        resid_cell_df,
        crosstalk_cell_df,
        panel_feature_names,
        cofactors,
        panel,
        sample,
        sample_meta,
        dye_map=dye_map,
        correction_models=correction_models,
        crosstalk_cofactors=crosstalk_cofactors,
        crosstalk_metadata=crosstalk_metadata,
        pipeline_params=pipeline_params,
    )
    h5ad_path = os.path.join(panel_out, f"{panel}.h5ad")
    adata.write_h5ad(h5ad_path)
    print(f"  Saved panel h5ad → {h5ad_path} ({adata.n_obs} cells × {adata.n_vars} features)")
    print(f"  Panel {panel} done in {(B.datetime.now() - t0).total_seconds():.0f}s")

    _cleanup_memory(
        cell_df_all, pixel_df_all, cell_df_clean, resid_cell_df,
        crosstalk_cell_df, pixel_df_clean, pixel_df_corr, resid_summary, raw_layers, arcsinh_layers,
        resid_raw_layers, resid_arcsinh_layers,
    )
    B.slide_io.kill_jvm()
    return adata


def process_single_sample_nosom(
    root_dir: str,
    sample: str,
    fluo_dir: str,
    af_dir: str,
    raw_qptiff_dir: Optional[str],
    crosstalk_config_path: Optional[str],
    output_dir: str,
    registrar_dir: str,
    sample_info_file: str,
    bg_method: Optional[str] = None,
    bg_kwargs: Optional[Dict] = None,
    merge_overlap_threshold: float = 200.0,
    pixels_per_cell_residual: int = 0,
    af_corr_threshold: float = 0.20,
    donor_corr_threshold: float = 0.15,
    max_fit_pixels: int = 2_000_000,
    background_fit_quantile: float = 0.10,
) -> bool:
    print(f"\\n{'=' * 80}\\nSAMPLE: {sample}\\n{'=' * 80}")

    sample_fluo = os.path.join(fluo_dir, sample)
    sample_af = os.path.join(af_dir, sample)
    sample_mask = os.path.join(root_dir, sample, "mask")
    sample_out = os.path.join(output_dir, sample)
    B._safe_mkdir(sample_out)

    fluo_list = sorted(glob.glob(os.path.join(sample_fluo, "*.ome.tiff")))
    af_list = sorted(glob.glob(os.path.join(sample_af, "*.ome.tiff")))
    if not os.path.exists(sample_mask):
        print(f"  [SKIP] missing mask dir: {sample_mask}")
        return False
    if not fluo_list:
        print(f"  [SKIP] no fluo images in {sample_fluo}")
        return False
    if not af_list:
        print(f"  [SKIP] no AF images in {sample_af}")
        return False
    if len(fluo_list) != len(af_list):
        print(f"  [SKIP] fluo ({len(fluo_list)}) vs AF ({len(af_list)}) mismatch")
        return False

    sample_info_df = B.load_sample_info(sample_info_file)
    sample_meta = B.get_sample_metadata(sample, sample_info_df)
    s_ver = B.sample_version(sample)

    panel_h5ads = []
    n_panel_ok = 0
    n_panel_fail = 0
    for fluo_image in fluo_list:
        basename = os.path.splitext(os.path.splitext(os.path.basename(fluo_image))[0])[0]
        panel = re.sub("-", "_", basename)
        af_image = os.path.join(sample_af, f"{basename}.ome.tiff")
        raw_qptiff_path = None
        if raw_qptiff_dir:
            candidate = os.path.join(raw_qptiff_dir, sample, f"{basename}.qptiff")
            if os.path.exists(candidate):
                raw_qptiff_path = candidate
        if not os.path.exists(af_image):
            print(f"  [SKIP] missing AF image: {af_image}")
            n_panel_fail += 1
            continue
        print(f"  Processing panel file: {fluo_image}")
        try:
            process_panel_nosom(
                fluo_image,
                af_image,
                raw_qptiff_path,
                crosstalk_config_path,
                panel,
                sample,
                s_ver,
                sample_mask,
                output_dir,
                sample_meta,
                bg_method=bg_method,
                bg_kwargs=bg_kwargs,
                pixels_per_cell_residual=pixels_per_cell_residual,
                af_corr_threshold=af_corr_threshold,
                donor_corr_threshold=donor_corr_threshold,
                max_fit_pixels=max_fit_pixels,
                background_fit_quantile=background_fit_quantile,
            )
        except Exception as e:
            n_panel_fail += 1
            print(f"  [ERROR] panel {panel} failed in sample {sample}: {e}")
            import traceback
            traceback.print_exc()
            continue

        h5ad_path = os.path.join(sample_out, f"{panel}.h5ad")
        if os.path.exists(h5ad_path):
            panel_h5ads.append(h5ad_path)
            n_panel_ok += 1
        else:
            n_panel_fail += 1
            print(f"  [WARN] panel {panel} finished without output h5ad: {h5ad_path}")

    if not panel_h5ads:
        print(f"  [FAIL] no panel h5ad produced for {sample}")
        return False

    registrar_file = (
        Path(registrar_dir) / sample / sample / "data" / f"{sample}_registrar.pickle"
    )
    if not registrar_file.exists():
        print(f"  [WARN] missing registrar for {sample}; skip sample-level merged h5ad: {registrar_file}")
        print(f"  Sample summary: panels ok={n_panel_ok}, failed={n_panel_fail}, merged_h5ad=SKIPPED")
        return True

    try:
        B.merge_panel_h5ads_to_sample_h5ad(
            sample=sample,
            sample_h5ad_dir=sample_out,
            registrar_dir=registrar_dir,
            sample_meta=sample_meta,
            output_h5ad=os.path.join(sample_out, f"{sample}_merged_panels.h5ad"),
            panel_h5ads=panel_h5ads,
            overlap_threshold=merge_overlap_threshold,
        )
        print(f"  Sample summary: panels ok={n_panel_ok}, failed={n_panel_fail}, merged_h5ad=OK")
    except Exception as e:
        print(f"  [WARN] sample-level merge failed for {sample}: {e}")
        import traceback
        traceback.print_exc()
        print(f"  Sample summary: panels ok={n_panel_ok}, failed={n_panel_fail}, merged_h5ad=FAILED")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Residual-only no-SOM mIF pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--root_dir", type=str, required=True)
    parser.add_argument("--fluo_dir", type=str, required=True)
    parser.add_argument("--af_dir", type=str, required=True)
    parser.add_argument("--raw_qptiff_dir", type=str, default=None)
    parser.add_argument("--crosstalk_config_path", type=str, default=None)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--registrar_dir", type=str, required=True)
    parser.add_argument("--sample_info_file", type=str, required=True)
    parser.add_argument("--samples", nargs="+", default=None)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--ignore_blacklist",
        action="store_true",
        help="Use the provided sample list verbatim instead of dropping legacy BLACKLIST samples.",
    )
    parser.add_argument(
        "--bg_method",
        type=str,
        default=None,
        choices=["percentile", "mode", "rolling_ball", "morphological", "otsu"],
    )
    parser.add_argument("--bg_percentile", type=float, default=B.BG_PERCENTILE_DEFAULT)
    parser.add_argument("--merge_overlap_threshold", type=float, default=200.0)
    parser.add_argument(
        "--pixels_per_cell_residual",
        type=int,
        default=B.PIXELS_PER_CELL_RESIDUAL,
        help="Pixels per cell for residual correction / summary; <=0 means all pixels",
    )
    parser.add_argument("--af_corr_threshold", type=float, default=B.AF_CORR_THRESHOLD)
    parser.add_argument("--donor_corr_threshold", type=float, default=B.DONOR_CORR_THRESHOLD)
    parser.add_argument("--max_fit_pixels", type=int, default=B.MAX_FIT_PIXELS)
    parser.add_argument("--background_fit_quantile", type=float, default=B.BACKGROUND_FIT_QUANTILE)

    args = parser.parse_args()

    samples = args.samples
    if samples is None:
        samples = [d for d in os.listdir(args.root_dir) if os.path.isdir(os.path.join(args.root_dir, d))]
    if not args.ignore_blacklist:
        samples = [s for s in samples if s not in B.BLACKLIST]

    bg_kwargs = (
        {"percentile": args.bg_percentile}
        if args.bg_method == "percentile"
        else ({} if args.bg_method else None)
    )

    print(f"Samples ({len(samples)}): {samples}")
    print(f"Workers: {args.workers}")
    print(f"Residual-only pipeline; no SOM will be run.")
    print(f"pixels_per_cell_residual: {args.pixels_per_cell_residual}")
    print(f"AF/donor corr thresholds: {args.af_corr_threshold}/{args.donor_corr_threshold}")
    print(f"max_fit_pixels: {args.max_fit_pixels}; background_fit_quantile: {args.background_fit_quantile}")
    print(f"raw_qptiff_dir: {args.raw_qptiff_dir}")
    print(f"crosstalk_config_path: {args.crosstalk_config_path}")

    if args.workers <= 1:
        for sample in samples:
            try:
                ok = process_single_sample_nosom(
                    root_dir=args.root_dir,
                    sample=sample,
                    fluo_dir=args.fluo_dir,
                    af_dir=args.af_dir,
                    raw_qptiff_dir=args.raw_qptiff_dir,
                    crosstalk_config_path=args.crosstalk_config_path,
                    output_dir=args.output_dir,
                    registrar_dir=args.registrar_dir,
                    sample_info_file=args.sample_info_file,
                    bg_method=args.bg_method,
                    bg_kwargs=bg_kwargs,
                    merge_overlap_threshold=args.merge_overlap_threshold,
                    pixels_per_cell_residual=args.pixels_per_cell_residual,
                    af_corr_threshold=args.af_corr_threshold,
                    donor_corr_threshold=args.donor_corr_threshold,
                    max_fit_pixels=args.max_fit_pixels,
                    background_fit_quantile=args.background_fit_quantile,
                )
                print(f"[{'OK' if ok else 'FAIL'}] {sample}")
            except Exception as e:
                print(f"[ERROR] {sample}: {e}")
                import traceback
                traceback.print_exc()
            _cleanup_memory()
        return

    # conservative parallel branch
    from concurrent.futures import ProcessPoolExecutor, as_completed

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                process_single_sample_nosom,
                args.root_dir,
                sample,
                args.fluo_dir,
                args.af_dir,
                args.raw_qptiff_dir,
                args.crosstalk_config_path,
                args.output_dir,
                args.registrar_dir,
                args.sample_info_file,
                args.bg_method,
                bg_kwargs,
                args.merge_overlap_threshold,
                args.pixels_per_cell_residual,
                args.af_corr_threshold,
                args.donor_corr_threshold,
                args.max_fit_pixels,
                args.background_fit_quantile,
            ): sample
            for sample in samples
        }
        for future in as_completed(futures):
            sample = futures[future]
            try:
                ok = future.result()
                print(f"[{'OK' if ok else 'FAIL'}] {sample}")
            except Exception as e:
                print(f"[ERROR] {sample}: {e}")
                import traceback
                traceback.print_exc()

    print("\\n" + "=" * 80)
    print("All samples processed.")
    print("=" * 80)


if __name__ == "__main__":
    main()
