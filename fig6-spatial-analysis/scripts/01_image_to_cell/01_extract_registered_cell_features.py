#!/usr/bin/env python3
"""
Integrated Pixie-style mIF pipeline.

What this script does
---------------------
1. Per panel:
   - Read fluorescence OME-TIFF + AF image + tiled nuclear masks
   - Optional per-tile local background subtraction
   - Extract per-cell summary intensities (mean/p50/p75/p90) for DAPI, AF and antibodies
   - Sample per-cell pixels for Pixie-style level-1/level-2 SOM clustering
   - Save one panel-level AnnData (.h5ad)
   - Save diagnostic plots, including raw + arcsinh distributions for all summary stats

2. Per sample:
   - Read the three panel-level h5ad files
   - Warp coordinates using VALIS registrar (same logic as merge_panel_and_phenotyping.py)
   - Merge all cells and all antibody features into one sample-level AnnData (.h5ad)
   - Add sample_info.csv columns into obs
   - Add qc_three_panel_overlap flag (do not drop cells by default)

Design choices
--------------
- No per-panel CSV outputs for cell intensities. Quantile intensities are stored in h5ad layers.
- Panel-level h5ad layers include raw/arcsinh matrices for mean, p50, p75, p90.
- Residual correction and residual summary can use all pixels, while SOM can use a separate per-cell subsample for speed.
- Sample-level h5ad merges those layers across panels using union-of-features.
- Missing features for a given panel are stored as NaN in the merged sample-level h5ad.
- Panel-level obsm['X_pixel_cluster_fractions'] is preserved, but not merged across panels because
  each panel is clustered independently.
"""

from __future__ import annotations

import argparse
import gc
import glob
import json
import os
import re
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import cv2
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import tifffile
from minisom import MiniSom
from pyvips import Image
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.ndimage import uniform_filter
from scipy.optimize import nnls
from numpy.linalg import LinAlgError
from scipy.sparse import issparse
from scipy.spatial import KDTree
from skimage import measure
from skimage.filters import threshold_multiotsu, threshold_otsu
from skimage.morphology import disk, opening
from skimage.segmentation import expand_labels
from valis import registration, slide_io

matplotlib.use("Agg")
warnings.filterwarnings("ignore", category=FutureWarning)

try:
    import cupy as cp

    _GPU = True
except ImportError:
    _GPU = False

# -----------------------------------------------------------------------------
# CONSTANTS
# -----------------------------------------------------------------------------
TILE_WH = 4096
OVERLAP = 50
MIN_WH = 100
EXPAND_DISTANCE = 10

PIXELS_PER_CELL_RESIDUAL = (
    0  # <=0 means use all pixels per cell for residual correction / summary
)
PIXELS_PER_CELL_SOM = (
    64  # per-cell pixel subsampling used only for SOM / clustering speed
)
PIXEL_SOM_GRID = (20, 20)
CELL_SOM_GRID = (15, 15)
SOM_SIGMA = 1.5
SOM_LR = 0.5
PIXEL_SOM_ITERS = 200_000
CELL_SOM_ITERS = 100_000
N_PIXEL_METACLUSTERS = 40
N_CELL_METACLUSTERS = 20
MORPH_WEIGHT = 0.3

ARCSINH_COFACTOR_DEFAULT = 5.0
CROSSTALK_ARCSINH_MIN_COFACTOR = 1e-4
NORM_LO_PCT = 0.5
NORM_HI_PCT = 99.5

AF_OTSU_CLASSES = 4
AF_THRESHOLD_FLOOR = 30
AF_THRESHOLD_CEIL = 50

BG_METHOD_DEFAULT = "percentile"
BG_PERCENTILE_DEFAULT = 5.0

DISTRIBUTION_STATS = ("mean", "p50", "p75", "p90")
DEFAULT_X_STAT = {"DAPI": "p50", "AF": "p90", "marker": "p90"}

# Pixel residual correction
AF_CORR_THRESHOLD = 0.20
DONOR_CORR_THRESHOLD = 0.15
MAX_FIT_PIXELS = 2_000_000
BACKGROUND_FIT_QUANTILE = 0.10
AF_FIT_BINS = 64
AF_SPLINE_KNOTS = 3

# Akoya dye order after DAPI. If a 5-antibody panel is present, 780 is omitted.
DYE_ORDER_6 = ["570", "690", "480", "620", "780", "520"]
DYE_ORDER_5 = ["570", "690", "480", "620", "520"]
OPAL_CHANNEL_ORDER = [
    "Opal 480",
    "Opal 520",
    "Opal 570",
    "Opal 620",
    "Opal 690",
    "Opal 780",
]
AF_CHANNEL_NAME = "Sample AF"
ADJACENT_DYE_MAP = {
    "DAPI": ["480"],
    "480": ["DAPI", "520"],
    "520": ["480", "570"],
    "570": ["520", "620"],
    "620": ["570", "690"],
    "690": ["620", "780"],
    "780": ["690"],
}

ANTIBODY_MAPPING = {
    "v1": {
        "01": {
            7: ["CD31", "CD20", "CD45", "CD21", "PanCK", "CD45RO"],
            6: ["CD31", "CD20", "CD45", "CD21", "CD45RO"],
        },
        "02": {7: ["CD8", "TIM3", "CD4", "PD1", "CD3", "FOXP3"]},
        "03": {7: ["HER2", "CD68", "CD163", "CD11c", "KI67", "HLADR"]},
    },
    "v2": {
        "01": {7: ["HLADR", "CD11c", "CD20", "CD163", "CD21", "CD68"]},
        "02": {7: ["CD31", "KI67", "PanCK", "HER2", "PDL1", "SMA"]},
        "03": {7: ["FOXP3", "CD4", "TIM3", "CD45RO", "PD1", "CD8"]},
    },
}

V1_SAMPLES = {
    "pt9-pre",
    "pt9-post",
    "pt10-pre",
    "pt10-post",
    "pt16-pre",
    "pt17-pre",
    "pt18-pre",
    "pt19-pre",
    "pt20-pre",
    "pt21-pre",
    "pt23-pre",
    "pt23-post",
    "pt24-post",
    "pt25-pre",
    "pt25-post",
    "pt26-pre",
    "pt26-post",
    "pt27-pre",
    "pt27-post",
    "pt28-pre",
    "pt28-post",
    "pt29-pre",
    "pt29-post",
}

BLACKLIST = {
    "pt26-post",
    "pt57-post",
    "pt37-post",
    "pt33-pre",
    "pt20-post",
    "pt18-post",
    "pt51-post",
    "pt45-post",
    "pt16-post",
}

# -----------------------------------------------------------------------------
# GENERAL UTILITIES
# -----------------------------------------------------------------------------


def _safe_mkdir(path: str | Path):
    os.makedirs(path, exist_ok=True)


def pyvips_to_numpy(image: Image) -> np.ndarray:
    fmt2dt = {
        "uchar": np.uint8,
        "char": np.int8,
        "ushort": np.uint16,
        "short": np.int16,
        "uint": np.uint32,
        "int": np.int32,
        "float": np.float32,
        "double": np.float64,
        "complex": np.complex64,
        "dpcomplex": np.complex128,
    }
    return np.ndarray(
        buffer=image.write_to_memory(),
        dtype=fmt2dt[image.format],
        shape=[image.height, image.width, image.bands],
    )


def quartiles(labeled_mask, intensity):
    return np.percentile(intensity[labeled_mask], q=(50, 75, 90))


def normalize_sample_name(sample: str) -> str:
    """Normalize pt9-pre -> pt09-pre; keep suffix unchanged."""
    m = re.match(r"^(pt)(\d+)(-.+)$", sample)
    if not m:
        return sample
    prefix, num, suffix = m.groups()
    return f"{prefix}{int(num):02d}{suffix}"


def sample_version(sample: str) -> str:
    return "v1" if sample in V1_SAMPLES else "v2"


def sample_key_from_info_row(row: pd.Series) -> str:
    sid = str(row["Sample_ID"])
    tp = str(row["Time_point"]).split("-")[0]
    return normalize_sample_name(f"{sid}-{tp}")


def load_sample_info(sample_info_file: str) -> pd.DataFrame:
    df = pd.read_csv(sample_info_file)
    if "Sample_ID" not in df.columns or "Time_point" not in df.columns:
        raise ValueError(
            "sample_info.csv must contain at least 'Sample_ID' and 'Time_point' columns"
        )
    df = df.copy()
    df["_sample_key"] = df.apply(sample_key_from_info_row, axis=1)
    return df


def get_sample_metadata(sample: str, sample_info_df: pd.DataFrame) -> Dict[str, object]:
    key = normalize_sample_name(sample)
    matched = sample_info_df[sample_info_df["_sample_key"] == key]
    if matched.empty:
        return {"sample_name": sample}
    row = matched.iloc[0].to_dict()
    row.pop("_sample_key", None)
    row["sample_name"] = sample
    return row


def _panel_name_for_registrar(panel_name: str) -> str:
    return panel_name.replace("_", "-")


def _to_dense(x) -> np.ndarray:
    if issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=np.float32)


def _cleanup_memory(*objs):
    """Best-effort cleanup for large intermediate objects in worker processes."""
    try:
        del objs
    except Exception:
        pass
    gc.collect()
    if _GPU:
        try:
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()
        except Exception:
            pass


# -----------------------------------------------------------------------------
# BACKGROUND ESTIMATION
# -----------------------------------------------------------------------------


def estimate_background_percentile(
    tile: np.ndarray, mask: Optional[np.ndarray] = None, percentile: float = 5.0
) -> np.ndarray:
    H, W, C = tile.shape
    if mask is not None and mask.sum() > 100:
        cell_mask = mask > 0
    else:
        cell_mask = np.ones((H, W), dtype=bool)
    bg = np.zeros(C, dtype=np.float32)
    for c in range(C):
        vals = tile[cell_mask, c]
        bg[c] = max(float(np.percentile(vals, percentile)), 1.0)
    return bg


def estimate_background_mode(
    tile: np.ndarray, mask: Optional[np.ndarray] = None, n_bins: int = 256
) -> np.ndarray:
    H, W, C = tile.shape
    if mask is not None and mask.sum() > 100:
        cell_mask = mask > 0
    else:
        cell_mask = np.ones((H, W), dtype=bool)
    bg = np.zeros(C, dtype=np.float32)
    for c in range(C):
        vals = tile[cell_mask, c]
        vmin, vmax = vals.min(), vals.max()
        if vmax - vmin < 1:
            bg[c] = float(vmin)
            continue
        hist, bins = np.histogram(vals, bins=n_bins, range=(vmin, vmax))
        i = int(np.argmax(hist))
        bg[c] = max(float((bins[i] + bins[i + 1]) / 2), 1.0)
    return bg


def estimate_background_rolling_ball(tile: np.ndarray, radius: int = 50) -> np.ndarray:
    H, W, C = tile.shape
    bg_img = np.zeros_like(tile, dtype=np.float32)
    for c in range(C):
        smoothed = uniform_filter(
            tile[..., c].astype(np.float32), size=radius * 2 + 1, mode="reflect"
        )
        bg_img[..., c] = np.minimum(smoothed, tile[..., c])
    return bg_img.mean(axis=(0, 1))


def estimate_background_morphological(tile: np.ndarray, radius: int = 30) -> np.ndarray:
    H, W, C = tile.shape
    bg_img = np.zeros_like(tile, dtype=np.float32)
    selem = disk(radius)
    for c in range(C):
        bg_img[..., c] = opening(tile[..., c].astype(np.float32), selem)
    return bg_img.mean(axis=(0, 1))


def estimate_background_otsu(
    tile: np.ndarray, mask: Optional[np.ndarray] = None
) -> np.ndarray:
    H, W, C = tile.shape
    if mask is not None and mask.sum() > 100:
        cell_mask = mask > 0
    else:
        cell_mask = np.ones((H, W), dtype=bool)
    bg = np.zeros(C, dtype=np.float32)
    for c in range(C):
        vals = tile[cell_mask, c]
        try:
            thr = threshold_otsu(vals)
            bg_vals = vals[vals <= thr]
            bg[c] = (
                float(np.median(bg_vals))
                if len(bg_vals) > 0
                else float(np.percentile(vals, 5))
            )
        except Exception:
            bg[c] = float(np.percentile(vals, 5))
    return bg


def subtract_local_background(
    tile: np.ndarray,
    mask: Optional[np.ndarray] = None,
    method: str = BG_METHOD_DEFAULT,
    **kwargs,
) -> Tuple[np.ndarray, np.ndarray]:
    tile = tile.astype(np.float32)
    if method == "percentile":
        bg = estimate_background_percentile(tile, mask, **kwargs)
    elif method == "mode":
        bg = estimate_background_mode(tile, mask, **kwargs)
    elif method == "rolling_ball":
        bg = estimate_background_rolling_ball(tile, **kwargs)
    elif method == "morphological":
        bg = estimate_background_morphological(tile, **kwargs)
    elif method == "otsu":
        bg = estimate_background_otsu(tile, mask, **kwargs)
    else:
        raise ValueError(f"Unknown background method: {method}")
    corrected = np.maximum(tile - bg[None, None, :], 0)
    return corrected, bg


class TileBackgroundTracker:
    def __init__(self, channel_names: List[str]):
        self.channel_names = channel_names
        self.records: List[Dict[str, object]] = []

    def add_tile(self, tile_id: str, x: int, y: int, backgrounds: Dict[str, float]):
        self.records.append(
            {"tile_id": tile_id, "tile_x": x, "tile_y": y, **backgrounds}
        )

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.records)

    def get_summary(self) -> Dict[str, Dict[str, float]]:
        df = self.to_dataframe()
        out: Dict[str, Dict[str, float]] = {}
        for ch in self.channel_names:
            col = f"{ch}_bg"
            if col not in df.columns:
                continue
            vals = df[col].dropna()
            if len(vals) == 0:
                continue
            out[ch] = {
                "mean": float(vals.mean()),
                "std": float(vals.std()),
                "min": float(vals.min()),
                "max": float(vals.max()),
                "median": float(vals.median()),
                "q25": float(vals.quantile(0.25)),
                "q75": float(vals.quantile(0.75)),
            }
        return out

    def flag_outliers(self, z_threshold: float = 3.0) -> pd.DataFrame:
        df = self.to_dataframe()
        flags = pd.DataFrame(index=df.index)
        for ch in self.channel_names:
            col = f"{ch}_bg"
            if col in df.columns:
                vals = df[col]
                z = (vals - vals.mean()) / (vals.std() + 1e-8)
                flags[f"{ch}_outlier"] = np.abs(z) > z_threshold
        flags["any_outlier"] = flags.any(axis=1) if len(flags.columns) else False
        return df.join(flags)

    def plot_backgrounds(self, output_path: str):
        df = self.to_dataframe()
        bg_cols = [f"{ch}_bg" for ch in self.channel_names if f"{ch}_bg" in df.columns]
        if not bg_cols:
            return
        n = len(bg_cols)
        n_cols = min(4, n)
        n_rows = int(np.ceil(n / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows))
        axes = np.ravel([axes]) if n == 1 else np.ravel(axes)
        for i, col in enumerate(bg_cols):
            ax = axes[i]
            vals = df[col].dropna()
            if len(vals):
                sns.histplot(vals, bins=30, kde=True, ax=ax)
                ax.axvline(
                    vals.median(),
                    color="red",
                    linestyle="--",
                    label=f"median={vals.median():.1f}",
                )
                ax.set_title(col.replace("_bg", ""))
                ax.legend(fontsize=8)
        for i in range(n, len(axes)):
            axes[i].set_visible(False)
        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()

    def plot_background_heatmap(self, output_path: str, img_w: int, img_h: int):
        df = self.to_dataframe()
        bg_cols = [f"{ch}_bg" for ch in self.channel_names if f"{ch}_bg" in df.columns]
        if not bg_cols or len(df) < 4:
            return
        fig, axes = plt.subplots(1, len(bg_cols), figsize=(4 * len(bg_cols), 4))
        axes = [axes] if len(bg_cols) == 1 else axes
        for ax, col in zip(axes, bg_cols):
            sc = ax.scatter(
                df["tile_x"], df["tile_y"], c=df[col], cmap="viridis", s=50, alpha=0.8
            )
            ax.set_xlim(0, img_w)
            ax.set_ylim(img_h, 0)
            ax.set_aspect("equal")
            ax.set_title(col.replace("_bg", " background"))
            plt.colorbar(sc, ax=ax, shrink=0.6)
        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()


# -----------------------------------------------------------------------------
# CELL FEATURE EXTRACTION
# -----------------------------------------------------------------------------


def cell_and_pixel_extraction(
    mask_np: np.ndarray,
    fluo_full_np: np.ndarray,
    af_tile_np: np.ndarray,
    antibodies: List[str],
    pixel_sample_k: int = PIXELS_PER_CELL_RESIDUAL,
    bg_method: Optional[str] = None,
    bg_kwargs: Optional[Dict] = None,
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[Dict[str, float]]]:
    mask_nuc = mask_np[:, :, 0]
    mask_cell = expand_labels(mask_nuc, distance=EXPAND_DISTANCE)

    tile_backgrounds = None
    if bg_method is not None:
        if bg_kwargs is None:
            bg_kwargs = {"percentile": BG_PERCENTILE_DEFAULT}
        fluo_corrected, fluo_bg = subtract_local_background(
            fluo_full_np, mask=mask_cell, method=bg_method, **bg_kwargs
        )
        af_corrected, af_bg = subtract_local_background(
            af_tile_np, mask=mask_cell, method=bg_method, **bg_kwargs
        )
        tile_backgrounds = {
            "DAPI_bg": float(fluo_bg[0]) if len(fluo_bg) > 0 else 0.0,
            "AF_bg": float(af_bg[0]) if len(af_bg) > 0 else 0.0,
        }
        for i, ab in enumerate(antibodies):
            if len(fluo_bg) > i + 1:
                tile_backgrounds[f"{ab}_bg"] = float(fluo_bg[i + 1])
        af_2d = af_corrected[:, :, 0]
        dapi_2d = fluo_corrected[..., 0]
        fluo_ab = fluo_corrected[..., 1:]
    else:
        af_2d = af_tile_np[:, :, 0]
        dapi_2d = fluo_full_np[..., 0]
        fluo_ab = fluo_full_np[..., 1:]

    props_nuc_af = measure.regionprops_table(
        mask_nuc,
        intensity_image=af_2d,
        cache=True,
        properties=(
            "label",
            "intensity_mean",
            "intensity_max",
            "intensity_min",
            "intensity_std",
        ),
        extra_properties=(quartiles,),
    )
    df_nuc_af = pd.DataFrame(props_nuc_af)
    df_nuc_af.columns = [
        "label",
        "AF_nuclear_mean",
        "AF_nuclear_max",
        "AF_nuclear_min",
        "AF_nuclear_std",
        "AF_nuclear_p50",
        "AF_nuclear_p75",
        "AF_nuclear_p90",
    ]

    props_cell_af = measure.regionprops_table(
        mask_cell,
        intensity_image=af_2d,
        cache=True,
        properties=(
            "label",
            "intensity_mean",
            "intensity_max",
            "intensity_min",
            "intensity_std",
        ),
        extra_properties=(quartiles,),
    )
    df_cell_af = pd.DataFrame(props_cell_af)
    df_cell_af.columns = [
        "label",
        "AF_cell_mean",
        "AF_cell_max",
        "AF_cell_min",
        "AF_cell_std",
        "AF_cell_p50",
        "AF_cell_p75",
        "AF_cell_p90",
    ]

    props_dapi = measure.regionprops_table(
        mask_nuc,
        intensity_image=dapi_2d,
        cache=True,
        properties=(
            "label",
            "intensity_mean",
            "intensity_max",
            "intensity_min",
            "intensity_std",
        ),
        extra_properties=(quartiles,),
    )
    df_dapi = pd.DataFrame(props_dapi)
    df_dapi.columns = [
        "label",
        "DAPI_mean",
        "DAPI_max",
        "DAPI_min",
        "DAPI_std",
        "DAPI_p50",
        "DAPI_p75",
        "DAPI_p90",
    ]

    result_df = None
    for c, antibody in enumerate(antibodies):
        props_nuc = measure.regionprops_table(
            mask_nuc,
            intensity_image=fluo_ab[..., c],
            cache=True,
            properties=(
                "label",
                "centroid",
                "area",
                "eccentricity",
                "intensity_mean",
                "intensity_max",
                "intensity_min",
                "intensity_std",
            ),
            extra_properties=(quartiles,),
        )
        tmp_nuc = pd.DataFrame(props_nuc)
        tmp_nuc.columns = [
            "label",
            "centroidY",
            "centroidX",
            "area_nuclear",
            "eccentricity",
            f"nuclear_mean_{antibody}",
            f"nuclear_max_{antibody}",
            f"nuclear_min_{antibody}",
            f"nuclear_std_{antibody}",
            f"nuclear_p50_{antibody}",
            f"nuclear_p75_{antibody}",
            f"nuclear_p90_{antibody}",
        ]
        props_cell = measure.regionprops_table(
            mask_cell,
            intensity_image=fluo_ab[..., c],
            properties=(
                "label",
                "area",
                "intensity_mean",
                "intensity_max",
                "intensity_min",
                "intensity_std",
            ),
            extra_properties=(quartiles,),
        )
        tmp_cell = pd.DataFrame(props_cell)
        tmp_cell.columns = [
            "label",
            "area_cell",
            f"cell_mean_{antibody}",
            f"cell_max_{antibody}",
            f"cell_min_{antibody}",
            f"cell_std_{antibody}",
            f"cell_p50_{antibody}",
            f"cell_p75_{antibody}",
            f"cell_p90_{antibody}",
        ]
        tmp = tmp_nuc.merge(tmp_cell, on="label", how="inner")
        src = "nuclear" if antibody in ("KI67", "FOXP3") else "cell"
        for stat in DISTRIBUTION_STATS:
            tmp[f"intensity_{stat}_{antibody}"] = tmp[f"{src}_{stat}_{antibody}"]
        if result_df is None:
            result_df = tmp
        else:
            result_df = result_df.merge(
                tmp,
                on=[
                    "label",
                    "centroidY",
                    "centroidX",
                    "area_nuclear",
                    "area_cell",
                    "eccentricity",
                ],
                how="inner",
            )

    result_df = df_dapi.merge(result_df, on="label", how="inner")
    result_df = df_nuc_af.merge(result_df, on="label", how="inner")
    result_df = df_cell_af.merge(result_df, on="label", how="inner")

    channel_names = ["AF", "DAPI"] + list(antibodies)
    all_ch = np.concatenate([af_2d[..., None], dapi_2d[..., None], fluo_ab], axis=-1)
    ab_channel_idx = list(range(2, 2 + len(antibodies)))
    pixel_records = []
    rng = np.random.default_rng(42)
    n_zero_filtered = 0

    for region in measure.regionprops(mask_cell):
        coords = region.coords
        n = len(coords)

        if pixel_sample_k is None or pixel_sample_k <= 0:
            if n < 2:
                continue
            sc = coords
        else:
            k = min(n, pixel_sample_k)
            if k < 2:
                continue
            idx = rng.choice(n, k, replace=False)
            sc = coords[idx]

        intensities = all_ch[sc[:, 0], sc[:, 1], :]
        ab_intensities = intensities[:, ab_channel_idx]
        nonzero_mask = ab_intensities.sum(axis=1) > 0

        if nonzero_mask.sum() == 0:
            n_zero_filtered += int(len(sc))
            continue

        n_zero_filtered += int((~nonzero_mask).sum())

        sc = sc[nonzero_mask]
        intensities = intensities[nonzero_mask]
        labels = np.full((len(sc), 1), region.label, dtype=np.int32)
        yx = sc.astype(np.float32)
        row = np.hstack([labels, yx, intensities])
        pixel_records.append(row)

    if pixel_records:
        pixel_arr = np.vstack(pixel_records)
        pixel_cols = ["cell_label", "pixel_y", "pixel_x"] + channel_names
        pixel_df = pd.DataFrame(pixel_arr, columns=pixel_cols)
        pixel_df["cell_label"] = pixel_df["cell_label"].astype(np.int32)
        if tile_backgrounds is not None:
            tile_backgrounds["n_zero_pixels_filtered"] = n_zero_filtered
    else:
        pixel_df = None

    return result_df, pixel_df, tile_backgrounds


# -----------------------------------------------------------------------------
# PIXEL QC + SOM
# -----------------------------------------------------------------------------


def pixel_level_qc(
    pixel_df: pd.DataFrame,
    antibody_cols: List[str],
    af_col: str = "AF",
    af_percentile_threshold: float = 95.0,
    pan_positive_threshold: float = 0.8,
    bad_pixel_cell_threshold: float = 0.5,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    pixel_df = pixel_df.copy()
    if af_col in pixel_df.columns:
        af_thr = np.percentile(pixel_df[af_col].values, af_percentile_threshold)
        high_af = pixel_df[af_col] > af_thr
    else:
        high_af = pd.Series(False, index=pixel_df.index)
    ab_cols_present = [c for c in antibody_cols if c in pixel_df.columns]
    if ab_cols_present:
        ab_medians = pixel_df[ab_cols_present].median()
        above_median = pixel_df[ab_cols_present] > ab_medians
        frac_above = above_median.sum(axis=1) / len(ab_cols_present)
        pan_positive = frac_above >= pan_positive_threshold
    else:
        pan_positive = pd.Series(False, index=pixel_df.index)
    pixel_df["pixel_qc_high_af"] = high_af
    pixel_df["pixel_qc_pan_positive"] = pan_positive
    pixel_df["pixel_qc_bad"] = high_af | pan_positive

    cell_stats = pixel_df.groupby("cell_label").agg(
        n_pixels=("pixel_qc_bad", "size"),
        n_bad_pixels=("pixel_qc_bad", "sum"),
        n_high_af=("pixel_qc_high_af", "sum"),
        n_pan_positive=("pixel_qc_pan_positive", "sum"),
    )
    cell_stats["bad_pixel_frac"] = cell_stats["n_bad_pixels"] / cell_stats["n_pixels"]
    cell_stats["qc_bad_pixel_dominant"] = (
        cell_stats["bad_pixel_frac"] >= bad_pixel_cell_threshold
    )
    return pixel_df, cell_stats.reset_index()


def infer_antibody_dye_map(antibodies: Sequence[str]) -> Dict[str, str]:
    if len(antibodies) == 6:
        dyes = DYE_ORDER_6
    elif len(antibodies) == 5:
        dyes = DYE_ORDER_5
    else:
        dyes = [f"CH{i + 1}" for i in range(len(antibodies))]
    return {ab: dyes[i] for i, ab in enumerate(antibodies)}


def infer_donor_map(antibodies: Sequence[str]) -> Dict[str, List[str]]:
    dye_map = infer_antibody_dye_map(antibodies)
    dye_to_feature = {"DAPI": "DAPI"}
    for ab, dye in dye_map.items():
        dye_to_feature[dye] = ab
    donor_map: Dict[str, List[str]] = {}
    for ab, dye in dye_map.items():
        donors = []
        for neighbor_dye in ADJACENT_DYE_MAP.get(dye, []):
            feat = dye_to_feature.get(neighbor_dye)
            if feat is not None and feat != ab:
                donors.append(feat)
        donor_map[ab] = donors
    return donor_map


def canonical_observed_channel_name(label: str) -> str:
    if label is None:
        return ""
    text = str(label).strip()
    token = re.sub(r"[^a-z0-9]+", "", text.lower())
    if token == "dapi":
        return "DAPI"
    if token in {"af", "sampleaf"}:
        return AF_CHANNEL_NAME
    match = re.search(r"(480|520|570|620|690|780)", text)
    if match:
        return f"Opal {match.group(1)}"
    return text


def _parse_qptiff_page_metadata(description: str) -> Dict[str, object]:
    clean = (description or "").replace("\x00", " ")
    name_match = re.search(r"<Name>(.*?)</Name>", clean, flags=re.IGNORECASE | re.DOTALL)
    exposure_match = re.search(r"<ExposureTime>(.*?)</ExposureTime>", clean, flags=re.IGNORECASE | re.DOTALL)
    image_type_match = re.search(r"<ImageType>(.*?)</ImageType>", clean, flags=re.IGNORECASE | re.DOTALL)
    try:
        exposure_value = float(exposure_match.group(1).strip()) if exposure_match else float("nan")
    except ValueError:
        exposure_value = float("nan")
    return {
        "name": name_match.group(1).strip() if name_match else "",
        "exposure_time_ms": exposure_value,
        "image_type": image_type_match.group(1).strip() if image_type_match else "",
    }


def load_qptiff_exposure_map(qptiff_path: str) -> Dict[str, float]:
    exposure_map: Dict[str, float] = {}
    with tifffile.TiffFile(qptiff_path) as tif:
        for page in tif.pages:
            meta = _parse_qptiff_page_metadata(page.description or "")
            if meta["image_type"] and str(meta["image_type"]) != "FullResolution":
                continue
            channel_name = canonical_observed_channel_name(str(meta["name"]))
            exposure_time_ms = float(meta["exposure_time_ms"])
            if not channel_name:
                continue
            exposure_map[channel_name] = exposure_time_ms
    return exposure_map


def load_crosstalk_correction_config(config_path: str) -> Dict[str, object]:
    with open(config_path, "r", encoding="utf-8") as handle:
        config = json.load(handle)
    required = ["opal_observed_channel_order", "opal_inverse_matrix", "af_zero_intercept_median"]
    missing = [key for key in required if key not in config]
    if missing:
        raise KeyError(f"Missing keys in crosstalk config {config_path}: {missing}")
    return config


def summarize_feature_columns_to_cells(
    pixel_df: pd.DataFrame,
    feature_to_column: Dict[str, str],
) -> pd.DataFrame:
    if pixel_df.empty:
        return pd.DataFrame(columns=["cell_label"])
    out = pd.DataFrame(index=np.sort(pixel_df["cell_label"].unique()))
    out.index.name = "cell_label"
    grouped = pixel_df.groupby("cell_label", observed=True)
    for feature, col in feature_to_column.items():
        if col not in pixel_df.columns:
            continue
        g = grouped[col]
        for stat_name, series in {
            "mean": g.mean(),
            "max": g.max(),
            "min": g.min(),
            "std": g.std().fillna(0),
            "p50": g.quantile(0.50),
            "p75": g.quantile(0.75),
            "p90": g.quantile(0.90),
        }.items():
            out[_feature_stat_col(feature, stat_name)] = series.astype(np.float32)
    out = out.reset_index()
    return out


def build_feature_to_observed_channel(
    antibodies: Sequence[str],
    dye_map: Dict[str, str],
) -> Dict[str, str]:
    mapping = {"DAPI": "DAPI", "AF": AF_CHANNEL_NAME}
    for antibody in antibodies:
        dye = str(dye_map.get(antibody, "")).strip()
        if dye:
            mapping[antibody] = canonical_observed_channel_name(dye)
    return mapping


def apply_crosstalk_correction_to_pixel_df(
    pixel_df: pd.DataFrame,
    antibodies: Sequence[str],
    dye_map: Dict[str, str],
    qptiff_path: str,
    correction_config: Dict[str, object],
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    if pixel_df.empty:
        raise ValueError("pixel_df is empty")
    exposure_map = load_qptiff_exposure_map(qptiff_path)
    feature_to_observed = build_feature_to_observed_channel(antibodies, dye_map)
    opal_channels = [canonical_observed_channel_name(ch) for ch in correction_config["opal_observed_channel_order"]]
    inverse = np.asarray(correction_config["opal_inverse_matrix"], dtype=np.float32)
    if inverse.shape[0] != inverse.shape[1] or inverse.shape[0] != len(opal_channels):
        raise ValueError("opal_inverse_matrix shape does not match opal_observed_channel_order")

    required_channels = {"DAPI", AF_CHANNEL_NAME}
    required_channels.update(ch for ch in feature_to_observed.values() if ch.startswith("Opal "))
    missing_exposure = sorted(
        ch for ch in required_channels
        if ch not in exposure_map or not np.isfinite(exposure_map[ch]) or exposure_map[ch] <= 0
    )
    if missing_exposure:
        raise ValueError(f"Missing or invalid exposure metadata for channels: {missing_exposure}")

    corrected = pixel_df.copy()
    dapi_vals = corrected["DAPI"].to_numpy(dtype=np.float32) / float(exposure_map["DAPI"])
    af_vals = corrected["AF"].to_numpy(dtype=np.float32) / float(exposure_map[AF_CHANNEL_NAME])
    observed = np.zeros((len(corrected), len(opal_channels)), dtype=np.float32)
    observed_to_feature = {observed_name: feature for feature, observed_name in feature_to_observed.items()}
    present_feature_channels: Dict[str, str] = {}
    for j, observed_name in enumerate(opal_channels):
        feature = observed_to_feature.get(observed_name)
        if feature is None or feature not in corrected.columns:
            continue
        present_feature_channels[feature] = observed_name
        observed[:, j] = corrected[feature].to_numpy(dtype=np.float32) / float(exposure_map[observed_name])

    af_coefficients = correction_config.get("af_zero_intercept_median", {})
    af_corrected_observed = observed.copy()
    for j, observed_name in enumerate(opal_channels):
        af_coeff = float(af_coefficients.get(observed_name, 0.0))
        af_corrected_observed[:, j] = np.clip(
            af_corrected_observed[:, j] - af_coeff * af_vals,
            0.0,
            None,
        )

    source_signal = af_corrected_observed @ inverse
    negative_fraction = float(np.mean(source_signal < 0))
    source_signal = np.clip(source_signal, 0.0, None).astype(np.float32)

    corrected["crosstalk_DAPI"] = dapi_vals.astype(np.float32)
    corrected["crosstalk_AF"] = af_vals.astype(np.float32)
    for j, observed_name in enumerate(opal_channels):
        feature = observed_to_feature.get(observed_name)
        if feature is None or feature not in corrected.columns:
            continue
        corrected[f"crosstalk_{feature}"] = source_signal[:, j]

    metadata = {
        "status": "ok",
        "qptiff_path": str(qptiff_path),
        "exposure_normalization_applied": True,
        "crosstalk_channel_order": opal_channels,
        "crosstalk_matrix_version": correction_config.get("version", "unknown"),
        "af_correction_source": "target_level_zero_intercept_median",
        "negative_fraction_before_clip": negative_fraction,
        "present_feature_channels": present_feature_channels,
        "exposure_map": {k: float(v) for k, v in exposure_map.items() if np.isfinite(v)},
    }
    return corrected, metadata


def _safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10:
        return 0.0
    x = x[m]
    y = y[m]
    if np.std(x) < 1e-8 or np.std(y) < 1e-8:
        return 0.0
    return float(np.corrcoef(np.arcsinh(x / 5.0), np.arcsinh(y / 5.0))[0, 1])


def _subsample_df(df: pd.DataFrame, max_rows: int, seed: int = 42) -> pd.DataFrame:
    if len(df) <= max_rows:
        return df
    return df.sample(n=max_rows, random_state=seed).reset_index(drop=True)


def subsample_pixels_per_cell(
    pixel_df: pd.DataFrame,
    pixels_per_cell: int,
    seed: int = 42,
) -> pd.DataFrame:
    """Subsample pixels *after* residual correction for SOM only.
    This keeps full-resolution pixels for residual summary while reducing downstream SOM cost.
    If pixels_per_cell <= 0, return the input unchanged.
    """
    if (
        pixel_df is None
        or len(pixel_df) == 0
        or pixels_per_cell is None
        or pixels_per_cell <= 0
    ):
        return pixel_df

    rng = np.random.default_rng(seed)
    out_parts = []
    for _, g in pixel_df.groupby("cell_label", sort=False, observed=True):
        if len(g) <= pixels_per_cell:
            out_parts.append(g)
        else:
            idx = rng.choice(len(g), size=pixels_per_cell, replace=False)
            out_parts.append(g.iloc[idx])
    if not out_parts:
        return pixel_df.iloc[:0].copy()
    return pd.concat(out_parts, axis=0, ignore_index=True)


def _select_background_pixels_for_target(
    fit_df: pd.DataFrame,
    target_col: str,
    af_col: str = "AF",
    n_bins: int = AF_FIT_BINS,
    quantile: float = BACKGROUND_FIT_QUANTILE,
) -> pd.DataFrame:
    fit_df = fit_df[
        [af_col, target_col]
        + [c for c in fit_df.columns if c not in [af_col, target_col]]
    ].copy()
    fit_df = fit_df[np.isfinite(fit_df[af_col]) & np.isfinite(fit_df[target_col])]
    if fit_df.empty:
        return fit_df
    try:
        bins = pd.qcut(
            fit_df[af_col], q=min(n_bins, fit_df[af_col].nunique()), duplicates="drop"
        )
    except ValueError:
        return fit_df
    fit_df["_af_bin"] = bins
    kept = []
    for _, g in fit_df.groupby("_af_bin", observed=True):
        thr = g[target_col].quantile(quantile)
        kept.append(g[g[target_col] <= thr])
    out = pd.concat(kept, axis=0, ignore_index=True) if kept else fit_df.iloc[:0].copy()
    return out.drop(columns=["_af_bin"], errors="ignore")



def _af_design_matrix(
    af: np.ndarray,
    donors: Optional[np.ndarray] = None,
    n_knots: int = AF_SPLINE_KNOTS,
    fixed_knots: Optional[Sequence[float]] = None,
) -> Tuple[np.ndarray, List[str], List[float]]:
    """
    Build AF spline + donor design matrix.

    Important:
    - During fitting, knots are learned from the fitting subset.
    - During application, the same knots can be reused on all pixels.
    """
    af = np.asarray(af, dtype=np.float32).reshape(-1)
    X_cols = [np.ones_like(af), af]
    names = ["intercept", "AF"]
    knot_vals: List[float] = []

    if fixed_knots is not None:
        knots = np.asarray(list(fixed_knots), dtype=np.float32).reshape(-1)
    elif n_knots > 0 and np.nanmax(af) > np.nanmin(af):
        qs = np.linspace(0.2, 0.8, n_knots)
        knots = np.unique(np.quantile(af, qs)).astype(np.float32)
    else:
        knots = np.asarray([], dtype=np.float32)

    for k in knots:
        X_cols.append(np.clip(af - float(k), 0, None))
        names.append(f"AF_relu_gt_{float(k):.4f}")
        knot_vals.append(float(k))

    if donors is not None and np.size(donors) > 0:
        donors = np.asarray(donors, dtype=np.float32)
        if donors.ndim == 1:
            donors = donors.reshape(-1, 1)
        for j in range(donors.shape[1]):
            X_cols.append(donors[:, j])

    X = np.column_stack(X_cols).astype(np.float32)
    return X, names, knot_vals


def _drop_degenerate_columns(
    X: np.ndarray,
    names: List[str],
    var_eps: float = 1e-10,
    corr_eps: float = 0.9999,
) -> Tuple[np.ndarray, List[str]]:
    """Remove near-constant and near-duplicate columns while always keeping intercept."""
    X = np.asarray(X, dtype=np.float64)
    keep = np.ones(X.shape[1], dtype=bool)

    for j in range(1, X.shape[1]):
        if np.nanvar(X[:, j]) < var_eps:
            keep[j] = False

    X1 = X[:, keep]
    names1 = [n for n, k in zip(names, keep) if k]

    keep2 = np.ones(X1.shape[1], dtype=bool)
    for j in range(1, X1.shape[1]):
        if not keep2[j]:
            continue
        xj = X1[:, j]
        sj = np.std(xj)
        if sj < 1e-12:
            keep2[j] = False
            continue
        for i in range(1, j):
            if not keep2[i]:
                continue
            xi = X1[:, i]
            si = np.std(xi)
            if si < 1e-12:
                continue
            r = np.corrcoef(xi, xj)[0, 1]
            if np.isfinite(r) and abs(r) >= corr_eps:
                keep2[j] = False
                break

    X2 = X1[:, keep2]
    names2 = [n for n, k in zip(names1, keep2) if k]
    return X2.astype(np.float32), names2


def _safe_nnls_with_fallback(
    X_fit: np.ndarray,
    y_fit: np.ndarray,
    coef_names: List[str],
    fit_bg: pd.DataFrame,
    af_col: str,
    selected_donors: List[str],
    af_spline_knots: int,
) -> Tuple[Optional[np.ndarray], Optional[List[str]], str]:
    """
    Try complex model first; if singular/degenerate, progressively simplify:
      1) AF + spline + donors
      2) AF + donors (no spline)
      3) AF only
      4) fail -> caller falls back to raw
    Returns: coef, final_coef_names, fit_mode
    """
    y_fit = np.asarray(y_fit, dtype=np.float32)

    def run_model(X, names):
        Xc, namesc = _drop_degenerate_columns(X, names)
        if Xc.shape[1] < 1 or len(namesc) < 1:
            raise LinAlgError("All columns dropped as degenerate")
        coef, _ = nnls(Xc, y_fit)
        return coef, namesc

    try:
        coef, final_names = run_model(X_fit, coef_names)
        return coef, final_names, "full"
    except Exception:
        pass

    donor_fit_simple = (
        fit_bg[selected_donors].to_numpy(dtype=np.float32) if selected_donors else None
    )
    X_simple, names_simple, _ = _af_design_matrix(
        fit_bg[af_col].to_numpy(dtype=np.float32),
        donor_fit_simple,
        n_knots=0,
        fixed_knots=None,
    )
    names_simple = names_simple + [f"donor_{d}" for d in selected_donors]
    try:
        coef, final_names = run_model(X_simple, names_simple)
        return coef, final_names, "af_plus_donor_no_spline"
    except Exception:
        pass

    X_af = np.column_stack(
        [
            np.ones(len(fit_bg), dtype=np.float32),
            fit_bg[af_col].to_numpy(dtype=np.float32),
        ]
    )
    names_af = ["intercept", "AF"]
    try:
        coef, final_names = run_model(X_af, names_af)
        return coef, final_names, "af_only"
    except Exception:
        pass

    return None, None, "failed"


def _build_application_matrix(
    corrected: pd.DataFrame,
    af_col: str,
    selected_donors: List[str],
    fit_mode: str,
    af_spline_knots: int,
    knots: Sequence[float],
) -> Tuple[np.ndarray, List[str]]:
    """Build X_all for the chosen fit mode."""
    if fit_mode == "full":
        donor_all = (
            corrected[selected_donors].to_numpy(dtype=np.float32)
            if selected_donors
            else None
        )
        X_all, names_all, _ = _af_design_matrix(
            corrected[af_col].to_numpy(dtype=np.float32),
            donor_all,
            n_knots=af_spline_knots,
            fixed_knots=knots,
        )
        names_all = names_all + [f"donor_{d}" for d in selected_donors]
        return X_all, names_all

    if fit_mode == "af_plus_donor_no_spline":
        donor_all = (
            corrected[selected_donors].to_numpy(dtype=np.float32)
            if selected_donors
            else None
        )
        X_all, names_all, _ = _af_design_matrix(
            corrected[af_col].to_numpy(dtype=np.float32),
            donor_all,
            n_knots=0,
            fixed_knots=None,
        )
        names_all = names_all + [f"donor_{d}" for d in selected_donors]
        return X_all, names_all

    if fit_mode == "af_only":
        X_all = np.column_stack(
            [
                np.ones(len(corrected), dtype=np.float32),
                corrected[af_col].to_numpy(dtype=np.float32),
            ]
        )
        names_all = ["intercept", "AF"]
        return X_all, names_all

    raise ValueError(f"Unknown fit_mode: {fit_mode}")


def fit_apply_pixel_residual_models(
    pixel_df: pd.DataFrame,
    antibodies: Sequence[str],
    donor_map: Optional[Dict[str, List[str]]] = None,
    af_col: str = "AF",
    dapi_col: str = "DAPI",
    af_corr_threshold: float = AF_CORR_THRESHOLD,
    donor_corr_threshold: float = DONOR_CORR_THRESHOLD,
    max_fit_pixels: int = MAX_FIT_PIXELS,
    background_fit_quantile: float = BACKGROUND_FIT_QUANTILE,
    n_af_bins: int = AF_FIT_BINS,
    af_spline_knots: int = AF_SPLINE_KNOTS,
    seed: int = 42,
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, object]]]:
    if donor_map is None:
        donor_map = infer_donor_map(antibodies)

    corrected = pixel_df.copy()
    fit_pool = _subsample_df(
        pixel_df[[af_col, dapi_col] + list(antibodies)].copy(),
        max_fit_pixels,
        seed=seed,
    )
    model_info: Dict[str, Dict[str, object]] = {}

    for target in antibodies:
        target_vals = fit_pool[target].to_numpy(dtype=np.float32)
        af_vals = fit_pool[af_col].to_numpy(dtype=np.float32)
        corr_af = _safe_corr(target_vals, af_vals)
        donor_candidates = [
            d for d in donor_map.get(target, []) if d in fit_pool.columns
        ]
        donor_corrs = {
            d: _safe_corr(target_vals, fit_pool[d].to_numpy(dtype=np.float32))
            for d in donor_candidates
        }
        selected_donors = [
            d for d, r in donor_corrs.items() if abs(r) >= donor_corr_threshold
        ]
        apply_correction = abs(corr_af) >= af_corr_threshold or len(selected_donors) > 0

        if not apply_correction:
            corrected[f"resid_{target}"] = corrected[target].astype(np.float32)
            model_info[target] = {
                "applied": False,
                "reason": "low correlation",
                "corr_af": corr_af,
                "donor_corrs": donor_corrs,
                "selected_donors": [],
            }
            continue

        fit_cols = [af_col, target] + selected_donors
        fit_df = fit_pool[fit_cols].copy()
        fit_bg = _select_background_pixels_for_target(
            fit_df,
            target_col=target,
            af_col=af_col,
            n_bins=n_af_bins,
            quantile=background_fit_quantile,
        )
        if fit_bg.empty or len(fit_bg) < 100:
            corrected[f"resid_{target}"] = corrected[target].astype(np.float32)
            model_info[target] = {
                "applied": False,
                "reason": "insufficient background pixels",
                "corr_af": corr_af,
                "donor_corrs": donor_corrs,
                "selected_donors": selected_donors,
            }
            continue

        donor_fit = (
            fit_bg[selected_donors].to_numpy(dtype=np.float32)
            if selected_donors
            else None
        )
        X_fit, coef_names, knots = _af_design_matrix(
            fit_bg[af_col].to_numpy(dtype=np.float32),
            donor_fit,
            n_knots=af_spline_knots,
            fixed_knots=None,
        )
        coef_names = coef_names + [f"donor_{d}" for d in selected_donors]
        y_fit = fit_bg[target].to_numpy(dtype=np.float32)

        try:
            coef, final_coef_names, fit_mode = _safe_nnls_with_fallback(
                X_fit, y_fit, coef_names, fit_bg, af_col, selected_donors, af_spline_knots
            )
        except Exception:
            coef, final_coef_names, fit_mode = None, None, "failed"

        if coef is None or final_coef_names is None or fit_mode == "failed":
            corrected[f"resid_{target}"] = corrected[target].astype(np.float32)
            model_info[target] = {
                "applied": False,
                "reason": "nnls_singular_failed_all_fallbacks",
                "corr_af": corr_af,
                "donor_corrs": donor_corrs,
                "selected_donors": selected_donors,
            }
            continue

        try:
            X_all, names_all = _build_application_matrix(
                corrected=corrected,
                af_col=af_col,
                selected_donors=selected_donors,
                fit_mode=fit_mode,
                af_spline_knots=af_spline_knots,
                knots=knots,
            )
            name_to_idx = {n: i for i, n in enumerate(names_all)}
            used_idx = [name_to_idx[n] for n in final_coef_names if n in name_to_idx]
            if len(used_idx) != len(final_coef_names):
                raise ValueError("application_matrix_name_mismatch")
            X_all_used = X_all[:, used_idx]
            pred_all = (X_all_used @ coef).astype(np.float32)
            resid_all = np.clip(
                corrected[target].to_numpy(dtype=np.float32) - pred_all, 0, None
            )
        except Exception:
            corrected[f"resid_{target}"] = corrected[target].astype(np.float32)
            model_info[target] = {
                "applied": False,
                "reason": "application_matrix_build_failed",
                "corr_af": corr_af,
                "donor_corrs": donor_corrs,
                "selected_donors": selected_donors,
            }
            continue

        corrected[f"resid_{target}"] = resid_all
        corr_resid_af = _safe_corr(
            resid_all, corrected[af_col].to_numpy(dtype=np.float32)
        )
        model_info[target] = {
            "applied": True,
            "fit_mode": fit_mode,
            "corr_af": corr_af,
            "corr_resid_af": corr_resid_af,
            "donor_corrs": donor_corrs,
            "selected_donors": selected_donors,
            "coef_names": final_coef_names,
            "coef": [float(x) for x in coef],
            "af_knots": knots,
            "n_fit_pool": int(len(fit_pool)),
            "n_fit_background": int(len(fit_bg)),
            "background_fit_quantile": float(background_fit_quantile),
            "n_af_bins": int(n_af_bins),
        }
    return corrected, model_info


def summarise_residual_pixels_to_cells(
    pixel_df: pd.DataFrame, antibodies: Sequence[str]
) -> pd.DataFrame:
    if pixel_df.empty:
        return pd.DataFrame(columns=["cell_label"])
    out = pd.DataFrame(index=np.sort(pixel_df["cell_label"].unique()))
    out.index.name = "cell_label"
    grouped = pixel_df.groupby("cell_label", observed=True)
    for ab in antibodies:
        col = f"resid_{ab}"
        if col not in pixel_df.columns:
            continue
        g = grouped[col]
        out[f"intensity_mean_{ab}"] = g.mean().astype(np.float32)
        out[f"intensity_max_{ab}"] = g.max().astype(np.float32)
        out[f"intensity_min_{ab}"] = g.min().astype(np.float32)
        out[f"intensity_std_{ab}"] = g.std().fillna(0).astype(np.float32)
        out[f"intensity_p50_{ab}"] = g.quantile(0.50).astype(np.float32)
        out[f"intensity_p75_{ab}"] = g.quantile(0.75).astype(np.float32)
        out[f"intensity_p90_{ab}"] = g.quantile(0.90).astype(np.float32)
    out = out.reset_index()
    return out


def plot_raw_vs_residual_diagnostics(
    pixel_df: pd.DataFrame,
    antibodies: Sequence[str],
    model_info: Dict[str, Dict[str, object]],
    output_dir: str,
    max_points: int = 200_000,
    seed: int = 42,
):
    diag_dir = os.path.join(output_dir, "residual_diagnostics")
    _safe_mkdir(diag_dir)
    view = _subsample_df(pixel_df, max_points, seed=seed)
    for ab in antibodies:
        resid_col = f"resid_{ab}"
        if resid_col not in view.columns:
            continue
        info = model_info.get(ab, {})
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        axes[0].scatter(
            np.arcsinh(view["AF"] / 5.0), np.arcsinh(view[ab] / 5.0), s=1, alpha=0.15
        )
        axes[0].set_xlabel("AF (arcsinh)")
        axes[0].set_ylabel(f"{ab} raw (arcsinh)")
        axes[0].set_title(f"raw vs AF\nr={info.get('corr_af', 0):.3f}")

        axes[1].scatter(
            np.arcsinh(view["AF"] / 5.0),
            np.arcsinh(view[resid_col] / 5.0),
            s=1,
            alpha=0.15,
        )
        axes[1].set_xlabel("AF (arcsinh)")
        axes[1].set_ylabel(f"{ab} residual (arcsinh)")
        axes[1].set_title(f"residual vs AF\nr={info.get('corr_resid_af', 0):.3f}")

        axes[2].hist(np.arcsinh(view[ab] / 5.0), bins=100, alpha=0.5, label="raw")
        axes[2].hist(
            np.arcsinh(view[resid_col] / 5.0), bins=100, alpha=0.5, label="residual"
        )
        axes[2].set_title("distribution (arcsinh)")
        axes[2].legend()
        plt.suptitle(f"{ab} residual correction")
        plt.tight_layout()
        plt.savefig(
            os.path.join(diag_dir, f"{ab}_raw_vs_residual.png"),
            dpi=180,
            bbox_inches="tight",
        )
        plt.close()


def remove_duplicates_by_distance(
    df: pd.DataFrame, distance_threshold: float = 10.0
) -> Tuple[pd.DataFrame, np.ndarray]:
    coords = df[["centroidY", "centroidX"]].to_numpy()
    n = len(df)
    tree = KDTree(coords)
    pairs = tree.query_pairs(distance_threshold)
    if not pairs:
        return df.copy().reset_index(drop=True), np.arange(n)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            if ri < rj:
                parent[rj] = ri
            else:
                parent[ri] = rj

    for i, j in pairs:
        union(i, j)
    clusters: Dict[int, List[int]] = {}
    for i in range(n):
        clusters.setdefault(find(i), []).append(i)
    keep = sorted(min(v) for v in clusters.values())
    return df.iloc[keep].reset_index(drop=True), np.array(keep)


def estimate_cofactors(
    pixel_df: pd.DataFrame,
    feature_cols: List[str],
    already_bg_subtracted: bool = False,
    min_cofactor: Optional[float] = None,
    fallback_cofactor: Optional[float] = None,
) -> Dict[str, float]:
    cofactors = {}
    bg_subtracted_floor = (
        0.5 if min_cofactor is None else max(float(min_cofactor), 1e-12)
    )
    bg_subtracted_fallback = (
        ARCSINH_COFACTOR_DEFAULT / 2
        if fallback_cofactor is None
        else max(float(fallback_cofactor), 1e-12)
    )
    for col in feature_cols:
        vals = pixel_df[col].values
        if already_bg_subtracted:
            pos = vals[vals > 0]
            if len(pos) > 100:
                med = np.median(pos)
                mad = np.median(np.abs(pos - med))
                pseudo_sd = mad * 1.4826
                cofactors[col] = max(float(pseudo_sd), bg_subtracted_floor)
            else:
                cofactors[col] = bg_subtracted_fallback
        else:
            bg = vals[vals < np.percentile(vals, 50)]
            if len(bg) < 10:
                cofactors[col] = ARCSINH_COFACTOR_DEFAULT
                continue
            med = np.median(bg)
            mad = np.median(np.abs(bg - med))
            pseudo_sd = mad * 1.4826
            cofactors[col] = max(float(pseudo_sd), 1.0)
    return cofactors


def arcsinh_transform(X: np.ndarray, cofactors: np.ndarray) -> np.ndarray:
    return np.arcsinh(X / cofactors[None, :])


def percentile_normalize(
    X: np.ndarray, lo_pct: float = NORM_LO_PCT, hi_pct: float = NORM_HI_PCT
) -> np.ndarray:
    X_out = X.copy()
    for c in range(X.shape[1]):
        lo = np.percentile(X_out[:, c], lo_pct)
        hi = np.percentile(X_out[:, c], hi_pct)
        X_out[:, c] = np.clip((X_out[:, c] - lo) / (hi - lo + 1e-8), 0, 1)
    return X_out


def preprocess_pixel_matrix(
    pixel_df: pd.DataFrame,
    feature_cols: List[str],
    cofactors: Optional[Dict[str, float]] = None,
    already_bg_subtracted: bool = False,
) -> np.ndarray:
    X = pixel_df[feature_cols].values.astype(np.float64)
    if cofactors is None:
        cofactors = estimate_cofactors(pixel_df, feature_cols, already_bg_subtracted)
    cof_arr = np.array([cofactors[c] for c in feature_cols], dtype=np.float64)
    X = arcsinh_transform(X, cof_arr)
    X = percentile_normalize(X)
    return X.astype(np.float32)


def _assign_winners_cpu(
    flat_weights: np.ndarray, X: np.ndarray, batch: int = 100_000
) -> np.ndarray:
    n = len(X)
    winners = np.empty(n, dtype=np.int32)
    for s in range(0, n, batch):
        e = min(s + batch, n)
        d = np.sum((X[s:e, None, :] - flat_weights[None, :, :]) ** 2, axis=2)
        winners[s:e] = np.argmin(d, axis=1)
    return winners


def _assign_winners_gpu(
    flat_weights: np.ndarray, X: np.ndarray, batch: int = 500_000
) -> np.ndarray:
    w_gpu = cp.asarray(flat_weights)
    n = len(X)
    winners = np.empty(n, dtype=np.int32)
    for s in range(0, n, batch):
        e = min(s + batch, n)
        x_gpu = cp.asarray(X[s:e])
        d = cp.sum((x_gpu[:, None, :] - w_gpu[None, :, :]) ** 2, axis=2)
        winners[s:e] = cp.asnumpy(cp.argmin(d, axis=1))
    return winners


def assign_winners(som: MiniSom, X: np.ndarray) -> np.ndarray:
    weights = som.get_weights()
    flat_w = weights.reshape(-1, weights.shape[-1]).astype(np.float32)
    X32 = X.astype(np.float32)
    if _GPU:
        return _assign_winners_gpu(flat_w, X32)
    return _assign_winners_cpu(flat_w, X32)


def train_pixel_som(
    X: np.ndarray,
    grid: Tuple[int, int] = PIXEL_SOM_GRID,
    sigma: float = SOM_SIGMA,
    lr: float = SOM_LR,
    n_iter: int = PIXEL_SOM_ITERS,
    seed: int = 42,
    max_train_samples: int = 2_000_000,
) -> Tuple[MiniSom, np.ndarray]:
    som = MiniSom(
        grid[0],
        grid[1],
        X.shape[1],
        sigma=sigma,
        learning_rate=lr,
        random_seed=seed,
        neighborhood_function="gaussian",
    )
    if len(X) > 500:
        som.pca_weights_init(X[: min(len(X), 500_000)])
    else:
        som.random_weights_init(X)
    if len(X) > max_train_samples:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(X), max_train_samples, replace=False)
        X_train = X[idx]
    else:
        X_train = X
    som.train(X_train, n_iter, verbose=False)
    cluster_ids = assign_winners(som, X)
    return som, cluster_ids


def build_cell_cluster_fractions(
    cell_labels: np.ndarray, pixel_cluster_ids: np.ndarray, n_clusters: int
) -> Tuple[np.ndarray, np.ndarray]:
    unique_labels = np.unique(cell_labels)
    label_to_idx = {l: i for i, l in enumerate(unique_labels)}
    counts = np.zeros((len(unique_labels), n_clusters), dtype=np.float32)
    cell_idx = np.array([label_to_idx[l] for l in cell_labels], dtype=np.int32)
    np.add.at(counts, (cell_idx, pixel_cluster_ids), 1)
    row_sums = counts.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    return unique_labels, counts / row_sums


def train_cell_som(
    fractions: np.ndarray,
    grid: Tuple[int, int] = CELL_SOM_GRID,
    sigma: float = SOM_SIGMA,
    lr: float = SOM_LR,
    n_iter: int = CELL_SOM_ITERS,
    seed: int = 42,
) -> Tuple[MiniSom, np.ndarray]:
    som = MiniSom(
        grid[0],
        grid[1],
        fractions.shape[1],
        sigma=sigma,
        learning_rate=lr,
        random_seed=seed,
        neighborhood_function="gaussian",
    )
    if len(fractions) > 500:
        som.pca_weights_init(fractions[: min(len(fractions), 200_000)])
    else:
        som.random_weights_init(fractions)
    som.train(fractions, n_iter, verbose=False)
    cluster_ids = assign_winners(som, fractions)
    return som, cluster_ids


def metacluster_som(
    som: MiniSom, n_meta: int, method: str = "ward"
) -> Tuple[np.ndarray, np.ndarray]:
    weights = som.get_weights()
    flat_w = weights.reshape(-1, weights.shape[-1])
    Z = linkage(flat_w, method=method)
    node_meta = fcluster(Z, t=n_meta, criterion="maxclust") - 1
    return node_meta, Z


def map_to_metaclusters(cell_som_ids: np.ndarray, node_meta: np.ndarray) -> np.ndarray:
    return node_meta[cell_som_ids]


# -----------------------------------------------------------------------------
# DISTRIBUTION + QC PLOTS
# -----------------------------------------------------------------------------


def _feature_stat_col(feature: str, stat: str) -> str:
    if feature == "DAPI":
        return f"DAPI_{stat}"
    if feature == "AF":
        return f"AF_cell_{stat}"
    return f"intensity_{stat}_{feature}"


def _feature_kind(feature: str) -> str:
    if feature == "AF":
        return "AF"
    if feature == "DAPI":
        return "DAPI"
    return "marker"


def _build_stat_matrices(
    cell_df: pd.DataFrame,
    feature_names: List[str],
    cofactors: Dict[str, float],
    stats_to_store: Sequence[str] = DISTRIBUTION_STATS,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, str]]:
    raw_layers: Dict[str, np.ndarray] = {}
    arcsinh_layers: Dict[str, np.ndarray] = {}
    x_stat_map: Dict[str, str] = {}
    for stat in stats_to_store:
        mat = np.full((len(cell_df), len(feature_names)), np.nan, dtype=np.float32)
        for j, feat in enumerate(feature_names):
            col = _feature_stat_col(feat, stat)
            if col in cell_df.columns:
                mat[:, j] = cell_df[col].to_numpy(dtype=np.float32)
        raw_layers[f"raw_{stat}"] = mat
        cof_arr = np.array(
            [cofactors.get(f, ARCSINH_COFACTOR_DEFAULT) for f in feature_names],
            dtype=np.float32,
        )
        arcsinh_layers[f"arcsinh_{stat}"] = np.arcsinh(mat / cof_arr[None, :]).astype(
            np.float32
        )
    selected_raw = np.full((len(cell_df), len(feature_names)), np.nan, dtype=np.float32)
    for j, feat in enumerate(feature_names):
        kind = _feature_kind(feat)
        selected_stat = DEFAULT_X_STAT[kind]
        x_stat_map[feat] = selected_stat
        selected_raw[:, j] = raw_layers[f"raw_{selected_stat}"][:, j]
    cof_arr = np.array(
        [cofactors.get(f, ARCSINH_COFACTOR_DEFAULT) for f in feature_names],
        dtype=np.float32,
    )
    raw_layers["raw_selected"] = selected_raw
    arcsinh_layers["arcsinh_selected"] = np.arcsinh(
        selected_raw / cof_arr[None, :]
    ).astype(np.float32)
    return raw_layers, arcsinh_layers, x_stat_map


def plot_feature_distributions_single(
    raw_layers: Dict[str, np.ndarray],
    arcsinh_layers: Dict[str, np.ndarray],
    feature_names: List[str],
    output_dir: str,
    prefix: str,
    stats_to_plot: Sequence[str] = DISTRIBUTION_STATS,
):
    dist_dir = os.path.join(output_dir, f"distributions_{prefix}")
    _safe_mkdir(dist_dir)
    for j, feat in enumerate(feature_names):
        fig, axes = plt.subplots(
            2, len(stats_to_plot), figsize=(4 * len(stats_to_plot), 6)
        )
        for k, stat in enumerate(stats_to_plot):
            raw_vals = raw_layers[f"raw_{stat}"][:, j]
            arc_vals = arcsinh_layers[f"arcsinh_{stat}"][:, j]
            raw_vals = raw_vals[np.isfinite(raw_vals)]
            arc_vals = arc_vals[np.isfinite(arc_vals)]
            if len(raw_vals):
                axes[0, k].hist(raw_vals, bins=100, alpha=0.8)
            axes[0, k].set_title(f"{prefix} {stat}")
            axes[0, k].set_xlabel(feat)
            if len(arc_vals):
                axes[1, k].hist(arc_vals, bins=100, alpha=0.8)
            axes[1, k].set_title(f"arcsinh {prefix} {stat}")
            axes[1, k].set_xlabel(feat)
        plt.suptitle(f"{feat}: {prefix} raw vs arcsinh distributions")
        plt.tight_layout()
        plt.savefig(
            os.path.join(dist_dir, f"{feat}_{prefix}_raw_vs_arcsinh.png"),
            dpi=180,
            bbox_inches="tight",
        )
        plt.close()


def plot_feature_distributions(
    raw_layers: Dict[str, np.ndarray],
    arcsinh_layers: Dict[str, np.ndarray],
    feature_names: List[str],
    output_dir: str,
    stats_to_plot: Sequence[str] = DISTRIBUTION_STATS,
    residual_layers: Optional[Dict[str, np.ndarray]] = None,
    residual_arcsinh_layers: Optional[Dict[str, np.ndarray]] = None,
):
    plot_feature_distributions_single(
        raw_layers,
        arcsinh_layers,
        feature_names,
        output_dir,
        prefix="raw",
        stats_to_plot=stats_to_plot,
    )
    if residual_layers is not None and residual_arcsinh_layers is not None:
        plot_feature_distributions_single(
            residual_layers,
            residual_arcsinh_layers,
            feature_names,
            output_dir,
            prefix="resid",
            stats_to_plot=stats_to_plot,
        )


def plot_pixel_metacluster_heatmap(
    pixel_df: pd.DataFrame,
    feature_cols: List[str],
    pixel_meta_col: str,
    output_path: str,
    cofactors: Dict[str, float],
):
    X = pixel_df[feature_cols].values.astype(np.float64)
    cof_arr = np.array([cofactors[c] for c in feature_cols], dtype=np.float64)
    X_t = np.arcsinh(X / cof_arr[None, :])
    tmp = pd.DataFrame(X_t, columns=feature_cols)
    tmp["metacluster"] = pixel_df[pixel_meta_col].values
    medians = tmp.groupby("metacluster")[feature_cols].median()
    z = medians.apply(lambda col: (col - col.mean()) / (col.std() + 1e-8))
    cluster_counts = pixel_df[pixel_meta_col].value_counts().sort_index()
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(max(10, len(feature_cols) * 1.2), max(8, len(medians) * 0.4)),
        gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05},
    )
    sns.heatmap(
        z,
        cmap="RdBu_r",
        center=0,
        linewidths=0.3,
        linecolor="grey",
        annot=False,
        ax=axes[0],
        cbar_kws={"shrink": 0.6},
    )
    axes[0].set_title("Pixel metacluster × feature")
    y_pos = np.arange(len(cluster_counts))
    axes[1].barh(y_pos, cluster_counts.values, color="steelblue", height=0.8)
    axes[1].invert_yaxis()
    axes[1].set_yticks(y_pos)
    axes[1].set_yticklabels([])
    axes[1].set_xlabel("# pixels")
    axes[1].set_title("Cluster size")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_cell_metacluster_heatmap(
    cell_df: pd.DataFrame,
    feature_cols: List[str],
    meta_col: str,
    output_path: str,
    cofactors: Dict[str, float],
    flag_pan_positive: bool = True,
    pan_positive_threshold: float = 1.0,
) -> set:
    int_cols = []
    for feat in feature_cols:
        col = _feature_stat_col(feat, DEFAULT_X_STAT[_feature_kind(feat)])
        if col in cell_df.columns:
            int_cols.append((feat, col))
    if not int_cols:
        return set()
    feat_names, col_names = zip(*int_cols)
    X = cell_df[list(col_names)].values.astype(np.float64)
    cof_arr = np.array(
        [cofactors.get(f, ARCSINH_COFACTOR_DEFAULT) for f in feat_names],
        dtype=np.float64,
    )
    X_t = np.arcsinh(X / cof_arr[None, :])
    tmp = pd.DataFrame(X_t, columns=feat_names)
    tmp["metacluster"] = cell_df[meta_col].values
    medians = tmp.groupby("metacluster")[list(feat_names)].median()
    z_initial = medians.apply(lambda col: (col - col.mean()) / (col.std() + 1e-8))
    pan_positive_clusters = set()
    if flag_pan_positive and len(feat_names) >= 3:
        for mc in z_initial.index:
            n_high = (z_initial.loc[mc] > pan_positive_threshold).sum()
            if n_high >= len(feat_names) - 1:
                pan_positive_clusters.add(mc)
        if pan_positive_clusters:
            clean_clusters = [
                mc for mc in medians.index if mc not in pan_positive_clusters
            ]
            if len(clean_clusters) >= 2:
                medians_clean = medians.loc[clean_clusters]
                clean_mean = medians_clean.mean()
                clean_std = medians_clean.std() + 1e-8
                z = (medians - clean_mean) / clean_std
            else:
                z = z_initial
        else:
            z = z_initial
    else:
        z = z_initial
    cluster_counts = cell_df[meta_col].value_counts().sort_index()
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(max(10, len(feat_names) * 1.2), max(8, len(medians) * 0.4)),
        gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05},
    )
    sns.heatmap(
        z,
        cmap="RdBu_r",
        center=0,
        linewidths=0.3,
        linecolor="grey",
        annot=False,
        ax=axes[0],
        cbar_kws={"shrink": 0.6},
    )
    axes[0].set_title("Cell metacluster × feature")
    if pan_positive_clusters:
        ytick_labels = [
            f"{mc} *" if mc in pan_positive_clusters else str(mc) for mc in z.index
        ]
        axes[0].set_yticklabels(ytick_labels)
    y_pos = np.arange(len(cluster_counts))
    colors = [
        "salmon" if mc in pan_positive_clusters else "steelblue"
        for mc in cluster_counts.index
    ]
    axes[1].barh(y_pos, cluster_counts.values, color=colors, height=0.8)
    axes[1].invert_yaxis()
    axes[1].set_yticks(y_pos)
    axes[1].set_yticklabels([])
    axes[1].set_xlabel("# cells")
    axes[1].set_title("Cluster size")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close()
    return pan_positive_clusters


def plot_cell_composition_heatmap(
    cell_fractions: np.ndarray,
    cell_meta: np.ndarray,
    pixel_node_meta: np.ndarray,
    output_path: str,
):
    n_pixel_meta = int(pixel_node_meta.max()) + 1
    meta_frac = np.zeros((int(cell_meta.max()) + 1, n_pixel_meta), dtype=np.float32)
    for pm in range(n_pixel_meta):
        nodes = np.where(pixel_node_meta == pm)[0]
        if len(nodes) == 0:
            continue
        for cm in range(meta_frac.shape[0]):
            mask = cell_meta == cm
            if mask.sum() == 0:
                continue
            meta_frac[cm, pm] = cell_fractions[mask][:, nodes].sum(axis=1).mean()
    row_sums = meta_frac.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    meta_frac /= row_sums
    fig, ax = plt.subplots(
        figsize=(max(10, n_pixel_meta * 0.35), max(6, meta_frac.shape[0] * 0.4))
    )
    sns.heatmap(
        meta_frac, cmap="YlOrRd", linewidths=0.3, linecolor="grey", annot=False, ax=ax
    )
    ax.set_xlabel("Pixel metacluster")
    ax.set_ylabel("Cell metacluster")
    ax.set_title("Cell metacluster composition")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_qc_heatmap(
    cell_df: pd.DataFrame, morph_cols: List[str], meta_col: str, output_path: str
):
    cols = [c for c in morph_cols if c in cell_df.columns]
    if not cols:
        return
    medians = cell_df[[meta_col] + cols].groupby(meta_col)[cols].median()
    z = medians.apply(lambda col: (col - col.mean()) / (col.std() + 1e-8))
    fig, ax = plt.subplots(
        figsize=(max(6, len(cols) * 1.2), max(6, len(medians) * 0.4))
    )
    sns.heatmap(
        z,
        cmap="coolwarm",
        center=0,
        linewidths=0.3,
        linecolor="grey",
        annot=False,
        ax=ax,
    )
    ax.set_title("Morphology / QC summary")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_dendrogram(Z: np.ndarray, n_clusters: int, output_path: str):
    fig, ax = plt.subplots(figsize=(10, 6))
    dendrogram(Z, truncate_mode="lastp", p=min(50, n_clusters * 2), ax=ax)
    ax.set_title(f"Hierarchical clustering (cut={n_clusters})")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def plot_af_scatter(
    cell_df: pd.DataFrame,
    antibody_cols: List[str],
    plot_dir: str,
    cofactors: Dict[str, float],
):
    af_col = "AF_cell_p90" if "AF_cell_p90" in cell_df.columns else None
    if af_col is None:
        return
    for ab in antibody_cols:
        int_col = (
            f"intensity_p90_{ab}" if f"intensity_p90_{ab}" in cell_df.columns else None
        )
        if int_col is None:
            continue
        x = np.arcsinh(
            cell_df[af_col].values / cofactors.get("AF", ARCSINH_COFACTOR_DEFAULT)
        )
        y = np.arcsinh(
            cell_df[int_col].values / cofactors.get(ab, ARCSINH_COFACTOR_DEFAULT)
        )
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(x, y, s=1, alpha=0.3)
        ax.set_xlabel("AF (arcsinh)")
        ax.set_ylabel(f"{ab} (arcsinh p90)")
        ax.set_title(f"AF vs {ab}")
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f"scatter_AF_vs_{ab}.png"), dpi=150)
        plt.close()


def plot_dapi_per_metacluster(cell_df: pd.DataFrame, meta_col: str, output_path: str):
    if "DAPI_mean" not in cell_df.columns:
        return
    fig, ax = plt.subplots(figsize=(max(10, cell_df[meta_col].nunique() * 0.5), 5))
    cell_df.boxplot(column="DAPI_mean", by=meta_col, ax=ax)
    ax.set_title("DAPI per cell metacluster")
    plt.suptitle("")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


# -----------------------------------------------------------------------------
# BUILD PANEL H5AD
# -----------------------------------------------------------------------------


def build_panel_anndata(
    cell_df: pd.DataFrame,
    resid_cell_df: Optional[pd.DataFrame],
    feature_names: List[str],
    cell_fractions: np.ndarray,
    cell_meta: np.ndarray,
    pixel_node_meta: np.ndarray,
    cofactors: Dict[str, float],
    som_params: Dict,
    panel: str,
    sample: str,
    sample_meta: Dict[str, object],
    dye_map: Optional[Dict[str, str]] = None,
    correction_models: Optional[Dict[str, Dict[str, object]]] = None,
) -> ad.AnnData:
    raw_layers, arcsinh_layers, x_stat_map = _build_stat_matrices(
        cell_df, feature_names, cofactors
    )
    extra_layer_dict = {}
    if resid_cell_df is not None:
        resid_raw_layers, resid_arc_layers, _ = _build_stat_matrices(
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
        "qc_low_dapi_cluster",
        "qc_high_af_cluster",
        "qc_tiny_nucleus",
        "qc_huge_nucleus",
        "qc_pan_positive",
        "qc_bad_pixel_dominant",
        "bad_pixel_frac",
        "cell_metacluster",
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
    if "cell_metacluster" in obs.columns:
        obs["cell_metacluster"] = obs["cell_metacluster"].astype("category")
    obs.index = obs["Cell_ID"].astype(str)

    var = pd.DataFrame(index=pd.Index(feature_names, name="feature"))
    var["feature_kind"] = [_feature_kind(f) for f in feature_names]
    var["panel"] = panel
    var["arcsinh_cofactor"] = [
        cofactors.get(f, ARCSINH_COFACTOR_DEFAULT) for f in feature_names
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
    adata.obsm["X_pixel_cluster_fractions"] = cell_fractions
    adata.uns["som_params"] = som_params
    adata.uns["panel"] = panel
    adata.uns["sample_ID"] = sample
    adata.uns["pixel_node_metacluster"] = pixel_node_meta
    adata.uns["sample_metadata"] = sample_meta
    adata.uns["correction_models"] = correction_models
    adata.uns["dye_map"] = dye_map
    return adata


# -----------------------------------------------------------------------------
# PANEL PROCESSING
# -----------------------------------------------------------------------------


def process_panel(
    fluo_image: str,
    af_image: str,
    panel: str,
    sample: str,
    sample_V: str,
    sample_mask_dir: str,
    output_dir: str,
    sample_meta: Dict[str, object],
    pixel_som_grid: Tuple[int, int] = PIXEL_SOM_GRID,
    cell_som_grid: Tuple[int, int] = CELL_SOM_GRID,
    n_pixel_meta: int = N_PIXEL_METACLUSTERS,
    n_cell_meta: int = N_CELL_METACLUSTERS,
    bg_method: Optional[str] = None,
    bg_kwargs: Optional[Dict] = None,
    pixels_per_cell_residual: int = PIXELS_PER_CELL_RESIDUAL,
    pixels_per_cell_som: int = PIXELS_PER_CELL_SOM,
    af_corr_threshold: float = AF_CORR_THRESHOLD,
    donor_corr_threshold: float = DONOR_CORR_THRESHOLD,
    max_fit_pixels: int = MAX_FIT_PIXELS,
    background_fit_quantile: float = BACKGROUND_FIT_QUANTILE,
) -> Optional[ad.AnnData]:
    t0 = datetime.now()
    slide_io.kill_jvm()

    number_part = os.path.splitext(os.path.splitext(os.path.basename(fluo_image))[0])[
        0
    ].split("-")[-1]
    reader_cls_fluo = slide_io.get_slide_reader(fluo_image)
    reader_fluo = reader_cls_fluo(fluo_image)
    reader_cls_af = slide_io.get_slide_reader(af_image)
    reader_af = reader_cls_af(af_image)
    ch_names = reader_fluo.metadata.channel_names
    w_fluo, h_fluo = reader_fluo.metadata.slide_dimensions[0]
    w_af, h_af = reader_af.metadata.slide_dimensions[0]
    if w_fluo != w_af or h_fluo != h_af:
        print(f"[SKIP] dimension mismatch: {fluo_image}")
        slide_io.kill_jvm()
        return None
    try:
        antibodies = ANTIBODY_MAPPING[sample_V][number_part][len(ch_names)]
    except KeyError:
        raise ValueError(
            f"Channel count mismatch for {fluo_image}: channels={len(ch_names)}"
        )

    dye_map = infer_antibody_dye_map(antibodies)
    donor_map = infer_donor_map(antibodies)
    print(
        f"  Panel={panel}, antibodies={antibodies}, dyes={dye_map}, size=({w_fluo}×{h_fluo})"
    )
    if bg_method:
        print(f"  Background subtraction: {bg_method} {bg_kwargs}")
    print(
        f"  pixels_per_cell_residual={pixels_per_cell_residual} (<=0 means all pixels for residual)"
    )
    print(
        f"  pixels_per_cell_som={pixels_per_cell_som} (<=0 means all corrected pixels enter SOM)"
    )

    fluo_vips = reader_fluo.slide2vips(level=0)
    af_vips = reader_af.slide2vips(level=0)

    af_full = pyvips_to_numpy(af_vips)
    af_small = cv2.resize(
        af_full, dsize=(0, 0), fx=0.2, fy=0.2, interpolation=cv2.INTER_AREA
    )
    del af_full
    af_thr = threshold_multiotsu(af_small, classes=AF_OTSU_CLASSES)[1]
    af_thr = float(np.clip(af_thr, AF_THRESHOLD_FLOOR, AF_THRESHOLD_CEIL))
    del af_small

    all_cell_dfs: List[pd.DataFrame] = []
    all_pixel_dfs: List[pd.DataFrame] = []
    label_offset = 0
    bg_tracker = (
        TileBackgroundTracker(["DAPI", "AF"] + list(antibodies)) if bg_method else None
    )

    for y in range(0, h_fluo, TILE_WH - OVERLAP):
        for x in range(0, w_fluo, TILE_WH - OVERLAP):
            cw = min(TILE_WH, w_fluo - x)
            ch = min(TILE_WH, h_fluo - y)
            if cw <= MIN_WH or ch <= MIN_WH:
                continue
            mask_file = os.path.join(
                sample_mask_dir, f"{panel}_mask_nuclear_{y}_{x}.tif"
            )
            if not os.path.exists(mask_file):
                continue
            with Image.new_from_file(mask_file, access="sequential") as mask_vips:
                mask_np = pyvips_to_numpy(mask_vips)
            fluo_tile = pyvips_to_numpy(fluo_vips.extract_area(x, y, cw, ch))
            af_tile = pyvips_to_numpy(af_vips.extract_area(x, y, cw, ch))
            cell_df, pixel_df, tile_bg = cell_and_pixel_extraction(
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
            if cell_df is None or len(cell_df) == 0:
                continue
            max_lbl = int(cell_df["label"].max())
            cell_df["global_label"] = cell_df["label"] + label_offset
            cell_df["centroidX"] += x
            cell_df["centroidY"] += y
            if pixel_df is not None and len(pixel_df) > 0:
                pixel_df["cell_label"] = (
                    pixel_df["cell_label"].astype(np.int64) + label_offset
                )
                pixel_df["pixel_x"] += x
                pixel_df["pixel_y"] += y
                all_pixel_dfs.append(pixel_df)
            all_cell_dfs.append(cell_df)
            label_offset += max_lbl + 1
            _cleanup_memory(mask_np, fluo_tile, af_tile)

    _cleanup_memory(fluo_vips, af_vips)
    if not all_cell_dfs:
        print(f"  [SKIP] no cells for panel {panel}")
        slide_io.kill_jvm()
        return None

    panel_out = os.path.join(output_dir, sample)
    plot_dir = os.path.join(panel_out, f"plots_{panel}")
    _safe_mkdir(plot_dir)

    cell_df_all = pd.concat(all_cell_dfs, ignore_index=True)
    pixel_df_all = (
        pd.concat(all_pixel_dfs, ignore_index=True) if all_pixel_dfs else None
    )
    _cleanup_memory(all_cell_dfs, all_pixel_dfs)

    if bg_tracker is not None and len(bg_tracker.records) > 0:
        bg_df = bg_tracker.to_dataframe()
        bg_df.to_csv(
            os.path.join(panel_out, f"tile_backgrounds_{panel}.csv"), index=False
        )
        bg_tracker.plot_backgrounds(
            os.path.join(plot_dir, "background_distributions.png")
        )
        bg_tracker.plot_background_heatmap(
            os.path.join(plot_dir, "background_spatial_heatmap.png"), w_fluo, h_fluo
        )
        outlier_df = bg_tracker.flag_outliers()
        if outlier_df["any_outlier"].sum() > 0:
            outlier_df[outlier_df["any_outlier"]].to_csv(
                os.path.join(panel_out, f"tile_background_outliers_{panel}.csv"),
                index=False,
            )
        bg_summary = bg_tracker.get_summary()
        _cleanup_memory(bg_df, outlier_df)
    else:
        bg_summary = None

    cell_df_all.insert(0, "image_ID", panel)
    cell_df_all, _ = remove_duplicates_by_distance(cell_df_all)
    kept_labels = set(cell_df_all["global_label"].values)
    if pixel_df_all is not None:
        pixel_df_all = pixel_df_all[
            pixel_df_all["cell_label"].isin(kept_labels)
        ].reset_index(drop=True)

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
            cell_df_all.get("AF_cell_p90", pd.Series(False, index=cell_df_all.index))
            > af_thr
        )

    cell_df_clean = cell_df_all.copy()
    pixel_df_clean = pixel_df_all
    if pixel_df_clean is None or len(pixel_df_clean) < 100:
        print(f"  [SKIP] too few sampled pixels for panel {panel}")
        slide_io.kill_jvm()
        return None

    ab_cols_pure = [c for c in antibodies if c in pixel_df_clean.columns]
    pixel_df_clean, cell_pixel_qc = pixel_level_qc(pixel_df_clean, ab_cols_pure)
    cell_df_clean = cell_df_clean.merge(
        cell_pixel_qc[["cell_label", "bad_pixel_frac", "qc_bad_pixel_dominant"]],
        left_on="global_label",
        right_on="cell_label",
        how="left",
    )
    cell_df_clean["bad_pixel_frac"] = cell_df_clean["bad_pixel_frac"].fillna(0)
    cell_df_clean["qc_bad_pixel_dominant"] = cell_df_clean[
        "qc_bad_pixel_dominant"
    ].fillna(False)
    if "cell_label" in cell_df_clean.columns:
        cell_df_clean = cell_df_clean.drop(columns=["cell_label"])

    # Pixel-level AF / donor residual correction
    pixel_df_corr, correction_models = fit_apply_pixel_residual_models(
        pixel_df_clean,
        antibodies=ab_cols_pure,
        donor_map=donor_map,
        af_corr_threshold=af_corr_threshold,
        donor_corr_threshold=donor_corr_threshold,
        max_fit_pixels=max_fit_pixels,
        background_fit_quantile=background_fit_quantile,
    )
    plot_raw_vs_residual_diagnostics(
        pixel_df_corr, ab_cols_pure, correction_models, plot_dir
    )

    # Residual summaries from corrected pixels
    resid_summary = summarise_residual_pixels_to_cells(pixel_df_corr, ab_cols_pure)
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

    _cleanup_memory(pixel_df_clean, resid_summary)

    # Use residual-corrected marker pixels for SOM; keep AF/DAPI unchanged
    pixel_df_som = subsample_pixels_per_cell(
        pixel_df_corr,
        pixels_per_cell=pixels_per_cell_som,
        seed=42,
    ).copy()
    print(f"  Pixels: residual={len(pixel_df_corr):,}, SOM={len(pixel_df_som):,}")
    for ab in ab_cols_pure:
        resid_col = f"resid_{ab}"
        if resid_col in pixel_df_som.columns:
            pixel_df_som[ab] = pixel_df_som[resid_col].astype(np.float32)
    _cleanup_memory(pixel_df_corr)

    som_feature_cols = ["DAPI", "AF"] + ab_cols_pure
    som_feature_cols = [c for c in som_feature_cols if c in pixel_df_som.columns]
    cofactors = estimate_cofactors(
        pixel_df_som, som_feature_cols, already_bg_subtracted=True
    )
    X_pixel = preprocess_pixel_matrix(
        pixel_df_som, som_feature_cols, cofactors, already_bg_subtracted=True
    )

    pixel_som, pixel_cluster_ids = train_pixel_som(
        X_pixel, grid=pixel_som_grid, n_iter=PIXEL_SOM_ITERS
    )
    pixel_df_som["pixel_cluster"] = pixel_cluster_ids
    n_pixel_nodes = pixel_som_grid[0] * pixel_som_grid[1]
    pixel_node_meta, pixel_Z = metacluster_som(pixel_som, n_meta=n_pixel_meta)
    pixel_df_som["pixel_metacluster"] = pixel_node_meta[pixel_cluster_ids]

    cell_labels_arr = pixel_df_som["cell_label"].values.astype(np.int64)
    unique_labels, cell_fractions = build_cell_cluster_fractions(
        cell_labels_arr, pixel_cluster_ids, n_pixel_nodes
    )
    label_to_frac_idx = {l: i for i, l in enumerate(unique_labels)}
    cell_df_clean = cell_df_clean[
        cell_df_clean["global_label"].isin(set(unique_labels))
    ].copy()
    resid_cell_df = resid_cell_df[
        resid_cell_df["global_label"].isin(set(unique_labels))
    ].copy()
    frac_idx = np.array(
        [label_to_frac_idx[l] for l in cell_df_clean["global_label"].values]
    )
    cell_fractions_aligned = cell_fractions[frac_idx]

    morph_cols = [
        "area_nuclear",
        "area_cell",
        "eccentricity",
        "DAPI_mean",
        "AF_cell_p90",
    ]
    morph_cols = [c for c in morph_cols if c in cell_df_clean.columns]
    if morph_cols:
        morph_raw = cell_df_clean[morph_cols].values.astype(np.float32)
        morph_norm = np.empty_like(morph_raw)
        for i in range(morph_raw.shape[1]):
            lo = np.percentile(morph_raw[:, i], 1)
            hi = np.percentile(morph_raw[:, i], 99)
            morph_norm[:, i] = np.clip((morph_raw[:, i] - lo) / (hi - lo + 1e-8), 0, 1)
        morph_norm *= MORPH_WEIGHT
        cell_features = np.hstack([cell_fractions_aligned, morph_norm])
    else:
        cell_features = cell_fractions_aligned

    _cleanup_memory(cell_fractions, frac_idx, label_to_frac_idx, unique_labels, X_pixel)
    cell_som, cell_cluster_ids = train_cell_som(
        cell_features, grid=cell_som_grid, n_iter=CELL_SOM_ITERS
    )
    cell_node_meta, cell_Z = metacluster_som(cell_som, n_meta=n_cell_meta)
    cell_meta = map_to_metaclusters(cell_cluster_ids, cell_node_meta)
    cell_df_clean["cell_metacluster"] = cell_meta
    resid_cell_df["cell_metacluster"] = cell_meta

    if "DAPI_mean" in cell_df_clean.columns:
        mc_dapi = cell_df_clean.groupby("cell_metacluster")["DAPI_mean"].median()
        global_p10 = np.percentile(cell_df_clean["DAPI_mean"].values, 10)
        low_dapi_mcs = set(mc_dapi[mc_dapi < global_p10].index)
        cell_df_clean["qc_low_dapi_cluster"] = cell_df_clean["cell_metacluster"].isin(
            low_dapi_mcs
        )
    else:
        cell_df_clean["qc_low_dapi_cluster"] = False
    if "AF_cell_p90" in cell_df_clean.columns:
        mc_af = cell_df_clean.groupby("cell_metacluster")["AF_cell_p90"].median()
        af_qc_thr = (
            np.percentile(cell_df_clean["AF_cell_p90"].values, 90)
            if bg_method
            else af_thr
        )
        cell_df_clean["qc_high_af_cluster"] = cell_df_clean["cell_metacluster"].isin(
            set(mc_af[mc_af > af_qc_thr].index)
        )
    else:
        cell_df_clean["qc_high_af_cluster"] = False
    if "area_nuclear" in cell_df_clean.columns:
        mc_area = cell_df_clean.groupby("cell_metacluster")["area_nuclear"].median()
        tiny_thr = np.percentile(cell_df_clean["area_nuclear"].values, 3)
        huge_thr = np.percentile(cell_df_clean["area_nuclear"].values, 99)
        cell_df_clean["qc_tiny_nucleus"] = cell_df_clean["cell_metacluster"].isin(
            set(mc_area[mc_area < tiny_thr].index)
        )
        cell_df_clean["qc_huge_nucleus"] = cell_df_clean["cell_metacluster"].isin(
            set(mc_area[mc_area > huge_thr].index)
        )
    else:
        cell_df_clean["qc_tiny_nucleus"] = False
        cell_df_clean["qc_huge_nucleus"] = False

    cell_df_clean["qc_suspect"] = (
        cell_df_clean["qc_low_dapi_cluster"]
        | cell_df_clean["qc_high_af_cluster"]
        | cell_df_clean["qc_tiny_nucleus"]
        | cell_df_clean["qc_huge_nucleus"]
        | cell_df_clean["qc_bad_pixel_dominant"]
    )
    for qc_col in [
        "qc_low_dapi_cluster",
        "qc_high_af_cluster",
        "qc_tiny_nucleus",
        "qc_huge_nucleus",
        "qc_suspect",
    ]:
        resid_cell_df[qc_col] = cell_df_clean[qc_col].values

    panel_feature_names = ["AF", "DAPI"] + list(ab_cols_pure)
    raw_layers, arcsinh_layers, _ = _build_stat_matrices(
        cell_df_clean, panel_feature_names, cofactors
    )
    resid_raw_layers, resid_arcsinh_layers, _ = _build_stat_matrices(
        resid_cell_df, panel_feature_names, cofactors
    )
    try:
        plot_feature_distributions(
            raw_layers,
            arcsinh_layers,
            panel_feature_names,
            plot_dir,
            residual_layers=resid_raw_layers,
            residual_arcsinh_layers=resid_arcsinh_layers,
        )
    except Exception as e:
        print(f"  [WARNING] plot_feature_distributions failed for {panel}: {e}")

    plot_pixel_metacluster_heatmap(
        pixel_df_som,
        som_feature_cols,
        "pixel_metacluster",
        os.path.join(plot_dir, "heatmap_pixel_metacluster.png"),
        cofactors,
    )
    plot_cell_metacluster_heatmap(
        resid_cell_df,
        som_feature_cols,
        "cell_metacluster",
        os.path.join(plot_dir, "heatmap_cell_metacluster_all.png"),
        cofactors,
    )
    pan_positive_clusters = plot_cell_metacluster_heatmap(
        resid_cell_df,
        ab_cols_pure,
        "cell_metacluster",
        os.path.join(plot_dir, "heatmap_cell_metacluster_markers.png"),
        cofactors,
    )
    if pan_positive_clusters:
        cell_df_clean["qc_pan_positive"] = cell_df_clean["cell_metacluster"].isin(
            pan_positive_clusters
        )
        cell_df_clean["qc_suspect"] = (
            cell_df_clean["qc_suspect"] | cell_df_clean["qc_pan_positive"]
        )
    else:
        cell_df_clean["qc_pan_positive"] = False
    resid_cell_df["qc_pan_positive"] = cell_df_clean["qc_pan_positive"].values
    resid_cell_df["qc_suspect"] = cell_df_clean["qc_suspect"].values

    plot_cell_composition_heatmap(
        cell_fractions_aligned,
        cell_meta,
        pixel_node_meta,
        os.path.join(plot_dir, "heatmap_cell_vs_pixel_composition.png"),
    )
    plot_qc_heatmap(
        cell_df_clean,
        morph_cols,
        "cell_metacluster",
        os.path.join(plot_dir, "heatmap_qc_morphology.png"),
    )
    plot_dendrogram(
        pixel_Z, n_pixel_meta, os.path.join(plot_dir, "dendrogram_pixel_meta.png")
    )
    plot_dendrogram(
        cell_Z, n_cell_meta, os.path.join(plot_dir, "dendrogram_cell_meta.png")
    )
    plot_af_scatter(resid_cell_df, ab_cols_pure, plot_dir, cofactors)
    plot_dapi_per_metacluster(
        cell_df_clean,
        "cell_metacluster",
        os.path.join(plot_dir, "boxplot_dapi_per_metacluster.png"),
    )

    som_params = {
        "pixel_som_grid": list(pixel_som_grid),
        "cell_som_grid": list(cell_som_grid),
        "pixel_som_iters": PIXEL_SOM_ITERS,
        "cell_som_iters": CELL_SOM_ITERS,
        "n_pixel_metaclusters": n_pixel_meta,
        "n_cell_metaclusters": n_cell_meta,
        "arcsinh_cofactors": {k: float(v) for k, v in cofactors.items()},
        "af_threshold": float(af_thr),
        "pixels_per_cell_residual": pixels_per_cell_residual,
        "pixels_per_cell_som": pixels_per_cell_som,
        "som_feature_cols": som_feature_cols,
        "morph_cols_in_cell_som": morph_cols,
        "morph_weight": MORPH_WEIGHT if morph_cols else 0,
        "bg_method": bg_method,
        "bg_kwargs": bg_kwargs,
        "bg_summary": bg_summary,
        "distribution_stats": list(DISTRIBUTION_STATS),
        "default_x_stat": DEFAULT_X_STAT,
        "dye_map": dye_map,
        "donor_map": donor_map,
        "correction_thresholds": {
            "af_corr_threshold": af_corr_threshold,
            "donor_corr_threshold": donor_corr_threshold,
            "max_fit_pixels": max_fit_pixels,
            "background_fit_quantile": background_fit_quantile,
        },
    }

    adata = build_panel_anndata(
        cell_df_clean,
        resid_cell_df,
        panel_feature_names,
        cell_fractions_aligned,
        cell_meta,
        pixel_node_meta,
        cofactors,
        som_params,
        panel,
        sample,
        sample_meta,
        dye_map=dye_map,
        correction_models=correction_models,
    )
    h5ad_path = os.path.join(panel_out, f"{panel}.h5ad")
    adata.write_h5ad(h5ad_path)
    print(
        f"  Saved panel h5ad → {h5ad_path} ({adata.n_obs} cells × {adata.n_vars} features)"
    )
    print(f"  Panel {panel} done in {(datetime.now() - t0).total_seconds():.0f}s")

    _cleanup_memory(
        cell_df_all, pixel_df_all, cell_df_clean, resid_cell_df, pixel_df_som,
        raw_layers, arcsinh_layers, resid_raw_layers, resid_arcsinh_layers,
        cell_features, cell_fractions_aligned, pixel_cluster_ids, cell_cluster_ids,
        pixel_node_meta, pixel_Z, cell_Z
    )
    slide_io.kill_jvm()
    return adata


# -----------------------------------------------------------------------------
# MERGE PANEL H5ADs TO SAMPLE H5AD
# -----------------------------------------------------------------------------


def _warp_obs_coords(
    obs: pd.DataFrame, slide_obj, x_col: str = "centroidX", y_col: str = "centroidY"
) -> pd.DataFrame:
    obs = obs.copy()
    raw_coords = obs[[x_col, y_col]].to_numpy(dtype=np.float32)
    warped_xy = slide_obj.warp_xy(raw_coords)
    obs["orig_centroidX"] = obs[x_col].astype(np.float32)
    obs["orig_centroidY"] = obs[y_col].astype(np.float32)
    obs[x_col] = np.rint(warped_xy[:, 0]).astype(np.int32)
    obs[y_col] = np.rint(warped_xy[:, 1]).astype(np.int32)
    return obs


def compute_three_panel_overlap_flag(
    obs: pd.DataFrame,
    panel_col: str = "image_ID",
    x_col: str = "centroidX",
    y_col: str = "centroidY",
    threshold: float = 200.0,
    n_required_panels: int = 3,
) -> np.ndarray:
    coords = obs[[x_col, y_col]].to_numpy(dtype=np.float32)
    panel_ids = obs[panel_col].astype(str).to_numpy()
    tree = KDTree(coords)
    keep = np.zeros(len(obs), dtype=bool)
    for i, coord in enumerate(coords):
        nbr = tree.query_ball_point(coord, threshold)
        nbr = [j for j in nbr if j != i]
        nbr_panels = set(panel_ids[nbr])
        if len(nbr_panels) >= n_required_panels:
            keep[i] = True
    return keep



def merge_panel_h5ads_to_sample_h5ad(
    sample: str,
    sample_h5ad_dir: str,
    registrar_dir: str,
    sample_meta: Dict[str, object],
    output_h5ad: Optional[str] = None,
    panel_h5ads: Optional[Sequence[str]] = None,
    overlap_threshold: float = 200.0,
) -> ad.AnnData:
    sample_h5ad_dir = Path(sample_h5ad_dir)
    if panel_h5ads is None:
        panel_h5ads = sorted(
            str(p)
            for p in sample_h5ad_dir.glob("*.h5ad")
            if not p.name.endswith("_merged_panels.h5ad")
        )
    if len(panel_h5ads) == 0:
        raise FileNotFoundError(f"No panel h5ad found under: {sample_h5ad_dir}")

    registrar_file = (
        Path(registrar_dir) / sample / sample / "data" / f"{sample}_registrar.pickle"
    )
    if not registrar_file.exists():
        raise FileNotFoundError(f"Registrar file not found: {registrar_file}")
    registrar_obj = registration.load_registrar(str(registrar_file))

    panel_meta = []
    union_vars: List[str] = []
    seen_vars = set()
    union_layers: List[str] = []
    seen_layers = set()
    panel_feature_map: Dict[str, List[str]] = {}
    panel_files_map: Dict[str, str] = {}
    panel_cofactor_map: Dict[str, Dict[str, float]] = {}
    obs_parts = []
    n_total = 0

    # First pass: gather metadata, warped obs, union vars/layers
    for fp in panel_h5ads:
        adata = ad.read_h5ad(fp)
        panel = str(adata.uns.get("panel", Path(fp).stem))
        panel_reg = _panel_name_for_registrar(panel)
        slide_obj = registrar_obj.get_slide(panel_reg)

        obs = adata.obs.copy()
        if "Cell_ID" not in obs.columns:
            obs["Cell_ID"] = adata.obs_names.astype(str)
        obs["orig_Cell_ID"] = obs["Cell_ID"].astype(str)
        obs["sample_ID"] = sample
        obs["panel"] = panel
        obs["panel_registrar_name"] = panel_reg
        obs["source_h5ad"] = str(fp)
        for k, v in sample_meta.items():
            obs[k] = v
        obs = _warp_obs_coords(obs, slide_obj)
        obs.index = obs["orig_Cell_ID"].astype(str)

        panel_meta.append(
            {
                "file": str(fp),
                "panel": panel,
                "var_names": list(map(str, adata.var_names)),
                "layers": list(adata.layers.keys()),
                "n_obs": int(adata.n_obs),
            }
        )
        obs_parts.append(obs)
        n_total += int(adata.n_obs)

        panel_feature_map[panel] = list(map(str, adata.var_names))
        panel_files_map[panel] = str(fp)
        if "arcsinh_cofactor" in adata.var.columns:
            panel_cofactor_map[panel] = adata.var["arcsinh_cofactor"].astype(float).to_dict()
        elif "som_params" in adata.uns and "arcsinh_cofactors" in adata.uns["som_params"]:
            panel_cofactor_map[panel] = {
                str(k): float(v)
                for k, v in adata.uns["som_params"]["arcsinh_cofactors"].items()
            }
        else:
            panel_cofactor_map[panel] = {}

        for v in map(str, adata.var_names):
            if v not in seen_vars:
                seen_vars.add(v)
                union_vars.append(v)
        for layer_name in adata.layers.keys():
            if layer_name not in seen_layers:
                seen_layers.add(layer_name)
                union_layers.append(layer_name)

        _cleanup_memory(adata, obs)

    var_to_idx = {v: i for i, v in enumerate(union_vars)}
    merged_layers = {
        layer_name: np.full((n_total, len(union_vars)), np.nan, dtype=np.float32)
        for layer_name in union_layers
    }
    merged_X = np.full((n_total, len(union_vars)), np.nan, dtype=np.float32)

    # Second pass: fill matrices one panel at a time
    row_start = 0
    for meta in panel_meta:
        adata = ad.read_h5ad(meta["file"])
        n = adata.n_obs
        cols = [var_to_idx[str(v)] for v in adata.var_names]
        merged_X[row_start : row_start + n, cols] = _to_dense(adata.X)
        for layer_name in union_layers:
            if layer_name in adata.layers:
                merged_layers[layer_name][row_start : row_start + n, cols] = _to_dense(
                    adata.layers[layer_name]
                )
        row_start += n
        _cleanup_memory(adata)

    obs_merged = pd.concat(obs_parts, axis=0, join="outer")
    obs_merged = obs_merged.iloc[:n_total].copy()
    obs_merged["qc_three_panel_overlap"] = compute_three_panel_overlap_flag(
        obs_merged, threshold=overlap_threshold
    )

    var = pd.DataFrame(index=pd.Index(union_vars, name="feature"))
    var["present_in_panels"] = [
        ",".join([p for p, feats in panel_feature_map.items() if feat in feats])
        for feat in union_vars
    ]
    var["feature_kind"] = [_feature_kind(f) for f in union_vars]
    var["n_panels_present"] = (
        var["present_in_panels"].str.split(",").map(lambda x: len([i for i in x if i]))
    )

    merged = ad.AnnData(X=merged_X, obs=obs_merged, var=var)
    for layer_name, arr in merged_layers.items():
        merged.layers[layer_name] = arr
    merged.uns["sample_ID"] = sample
    merged.uns["sample_metadata"] = sample_meta
    merged.uns["panels_loaded"] = [m["panel"] for m in panel_meta]
    merged.uns["panel_files"] = panel_files_map
    merged.uns["panel_feature_map"] = panel_feature_map
    merged.uns["panel_arcsinh_cofactors"] = panel_cofactor_map
    merged.uns["merge_params"] = {
        "overlap_threshold": overlap_threshold,
        "note": (
            "Per-panel h5ad files were coordinate-warped with VALIS and merged by union of features. "
            "qc_three_panel_overlap marks cells with nearby cells from all three panels after warping. "
            "Per-panel pixel-fraction obsm matrices are not merged because panel SOMs are trained independently."
        ),
    }
    merged.obs_names = merged.obs["orig_Cell_ID"].astype(str)
    merged.obs_names_make_unique()

    if output_h5ad is None:
        output_h5ad = sample_h5ad_dir / f"{sample}_merged_panels.h5ad"
    merged.write_h5ad(output_h5ad)
    print(
        f"  Saved sample-level merged h5ad → {output_h5ad} ({merged.n_obs} cells × {merged.n_vars} features)"
    )
    _cleanup_memory(obs_parts, panel_meta, merged_layers, merged_X, obs_merged, var)
    return merged


# -----------------------------------------------------------------------------
# SAMPLE PROCESSING
# -----------------------------------------------------------------------------


def process_single_sample(
    root_dir: str,
    sample: str,
    fluo_dir: str,
    af_dir: str,
    output_dir: str,
    registrar_dir: str,
    sample_info_file: str,
    pixel_som_grid: Tuple[int, int] = PIXEL_SOM_GRID,
    cell_som_grid: Tuple[int, int] = CELL_SOM_GRID,
    n_pixel_meta: int = N_PIXEL_METACLUSTERS,
    n_cell_meta: int = N_CELL_METACLUSTERS,
    bg_method: Optional[str] = None,
    bg_kwargs: Optional[Dict] = None,
    merge_overlap_threshold: float = 200.0,
    pixels_per_cell_residual: int = PIXELS_PER_CELL_RESIDUAL,
    pixels_per_cell_som: int = PIXELS_PER_CELL_SOM,
    af_corr_threshold: float = AF_CORR_THRESHOLD,
    donor_corr_threshold: float = DONOR_CORR_THRESHOLD,
    max_fit_pixels: int = MAX_FIT_PIXELS,
    background_fit_quantile: float = BACKGROUND_FIT_QUANTILE,
) -> bool:
    print(f"\n{'=' * 80}\nSAMPLE: {sample}\n{'=' * 80}")

    sample_fluo = os.path.join(fluo_dir, sample)
    sample_af = os.path.join(af_dir, sample)
    sample_mask = os.path.join(root_dir, sample, "mask")
    sample_out = os.path.join(output_dir, sample)
    _safe_mkdir(sample_out)

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

    sample_info_df = load_sample_info(sample_info_file)
    sample_meta = get_sample_metadata(sample, sample_info_df)
    s_ver = sample_version(sample)

    panel_h5ads = []
    n_panel_ok = 0
    n_panel_fail = 0
    for fluo_image in fluo_list:
        basename = os.path.splitext(os.path.splitext(os.path.basename(fluo_image))[0])[0]
        panel = re.sub("-", "_", basename)
        af_image = os.path.join(sample_af, f"{basename}.ome.tiff")
        if not os.path.exists(af_image):
            print(f"  [SKIP] missing AF image: {af_image}")
            n_panel_fail += 1
            continue
        print(f"  Processing panel file: {fluo_image}")
        try:
            process_panel(
                fluo_image,
                af_image,
                panel,
                sample,
                s_ver,
                sample_mask,
                output_dir,
                sample_meta,
                pixel_som_grid=pixel_som_grid,
                cell_som_grid=cell_som_grid,
                n_pixel_meta=n_pixel_meta,
                n_cell_meta=n_cell_meta,
                bg_method=bg_method,
                bg_kwargs=bg_kwargs,
                pixels_per_cell_residual=pixels_per_cell_residual,
                pixels_per_cell_som=pixels_per_cell_som,
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
        print(
            f"  [WARN] missing registrar for {sample}; skip sample-level merged h5ad: {registrar_file}"
        )
        print(
            f"  Sample summary: panels ok={n_panel_ok}, failed={n_panel_fail}, merged_h5ad=SKIPPED"
        )
        return True

    try:
        merge_panel_h5ads_to_sample_h5ad(
            sample=sample,
            sample_h5ad_dir=sample_out,
            registrar_dir=registrar_dir,
            sample_meta=sample_meta,
            output_h5ad=os.path.join(sample_out, f"{sample}_merged_panels.h5ad"),
            panel_h5ads=panel_h5ads,
            overlap_threshold=merge_overlap_threshold,
        )
        print(
            f"  Sample summary: panels ok={n_panel_ok}, failed={n_panel_fail}, merged_h5ad=OK"
        )
    except Exception as e:
        print(f"  [WARN] sample-level merge failed for {sample}: {e}")
        import traceback
        traceback.print_exc()
        print(
            f"  Sample summary: panels ok={n_panel_ok}, failed={n_panel_fail}, merged_h5ad=FAILED"
        )

    return True


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------


def main(
    root_dir: str,
    fluo_dir: str,
    af_dir: str,
    output_dir: str,
    registrar_dir: str,
    sample_info_file: str,
    samples: Optional[List[str]] = None,
    workers: int = 4,
    pixel_som_grid: Tuple[int, int] = PIXEL_SOM_GRID,
    cell_som_grid: Tuple[int, int] = CELL_SOM_GRID,
    n_pixel_meta: int = N_PIXEL_METACLUSTERS,
    n_cell_meta: int = N_CELL_METACLUSTERS,
    bg_method: Optional[str] = None,
    bg_percentile: float = BG_PERCENTILE_DEFAULT,
    merge_overlap_threshold: float = 200.0,
    pixels_per_cell_residual: int = PIXELS_PER_CELL_RESIDUAL,
    pixels_per_cell_som: int = PIXELS_PER_CELL_SOM,
    af_corr_threshold: float = AF_CORR_THRESHOLD,
    donor_corr_threshold: float = DONOR_CORR_THRESHOLD,
    max_fit_pixels: int = MAX_FIT_PIXELS,
    background_fit_quantile: float = BACKGROUND_FIT_QUANTILE,
):
    if samples is None:
        samples = [
            d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))
        ]
    samples = [s for s in samples if s not in BLACKLIST]
    bg_kwargs = (
        {"percentile": bg_percentile}
        if bg_method == "percentile"
        else ({} if bg_method else None)
    )

    print(f"Samples ({len(samples)}): {samples}")
    print(f"Workers: {workers}")
    print(f"GPU: {'YES' if _GPU else 'NO'}")
    print(f"Pixel SOM grid: {pixel_som_grid}; Cell SOM grid: {cell_som_grid}")
    print(f"Pixel meta: {n_pixel_meta}; Cell meta: {n_cell_meta}")
    print(f"Background subtraction: {bg_method if bg_method else 'DISABLED'}")
    print(f"Sample info file: {sample_info_file}")
    print(f"Registrar dir: {registrar_dir}")
    print(f"pixels_per_cell_residual: {pixels_per_cell_residual}")
    print(f"pixels_per_cell_som: {pixels_per_cell_som}")
    print(f"AF/donor corr thresholds: {af_corr_threshold}/{donor_corr_threshold}")
    print(
        f"max_fit_pixels: {max_fit_pixels}; background_fit_quantile: {background_fit_quantile}"
    )

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                process_single_sample,
                root_dir,
                sample,
                fluo_dir,
                af_dir,
                output_dir,
                registrar_dir,
                sample_info_file,
                pixel_som_grid,
                cell_som_grid,
                n_pixel_meta,
                n_cell_meta,
                bg_method,
                bg_kwargs,
                merge_overlap_threshold,
                pixels_per_cell_residual,
                pixels_per_cell_som,
                af_corr_threshold,
                donor_corr_threshold,
                max_fit_pixels,
                background_fit_quantile,
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

    print("\n" + "=" * 80)
    print("All samples processed.")
    print("=" * 80)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Integrated Pixie-style pipeline with per-panel h5ad and merged sample-level h5ad",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--root_dir",
        type=str,
        required=True,
        help="Root directory containing sample mask folders",
    )
    parser.add_argument(
        "--fluo_dir",
        type=str,
        required=True,
        help="Directory containing fluorescence images by sample",
    )
    parser.add_argument(
        "--af_dir",
        type=str,
        required=True,
        help="Directory containing AF images by sample",
    )
    parser.add_argument(
        "--output_dir", type=str, required=True, help="Output directory"
    )
    parser.add_argument(
        "--registrar_dir",
        type=str,
        required=True,
        help="VALIS registrar root directory",
    )
    parser.add_argument(
        "--sample_info_file", type=str, required=True, help="Path to sample_info.csv"
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help="Optional subset of sample directories to process",
    )
    parser.add_argument(
        "--workers", type=int, default=4, help="Parallel worker processes"
    )

    parser.add_argument("--pixel_som_x", type=int, default=PIXEL_SOM_GRID[0])
    parser.add_argument("--pixel_som_y", type=int, default=PIXEL_SOM_GRID[1])
    parser.add_argument("--cell_som_x", type=int, default=CELL_SOM_GRID[0])
    parser.add_argument("--cell_som_y", type=int, default=CELL_SOM_GRID[1])
    parser.add_argument("--n_pixel_meta", type=int, default=N_PIXEL_METACLUSTERS)
    parser.add_argument("--n_cell_meta", type=int, default=N_CELL_METACLUSTERS)

    parser.add_argument(
        "--bg_method",
        type=str,
        default=None,
        choices=["percentile", "mode", "rolling_ball", "morphological", "otsu"],
        help="Local background estimation method; if omitted, no local background subtraction is applied",
    )
    parser.add_argument(
        "--bg_percentile",
        type=float,
        default=BG_PERCENTILE_DEFAULT,
        help="Percentile for bg_method=percentile",
    )
    parser.add_argument(
        "--merge_overlap_threshold",
        type=float,
        default=200.0,
        help="Radius used to compute qc_three_panel_overlap in merged h5ad",
    )
    parser.add_argument(
        "--pixels_per_cell_residual",
        type=int,
        default=PIXELS_PER_CELL_RESIDUAL,
        help="Pixels kept per cell for residual correction / summary; <=0 means use all pixels",
    )
    parser.add_argument(
        "--pixels_per_cell_som",
        type=int,
        default=PIXELS_PER_CELL_SOM,
        help="Pixels kept per cell for SOM only; <=0 means use all corrected pixels",
    )
    parser.add_argument(
        "--af_corr_threshold",
        type=float,
        default=AF_CORR_THRESHOLD,
        help="Min |corr(target, AF)| to trigger AF residual correction",
    )
    parser.add_argument(
        "--donor_corr_threshold",
        type=float,
        default=DONOR_CORR_THRESHOLD,
        help="Min |corr(target, donor)| to include donor in correction model",
    )
    parser.add_argument(
        "--max_fit_pixels",
        type=int,
        default=MAX_FIT_PIXELS,
        help="Max pixels used to fit residual models per panel",
    )
    parser.add_argument(
        "--background_fit_quantile",
        type=float,
        default=BACKGROUND_FIT_QUANTILE,
        help="Within each AF bin, fit using target pixels below this quantile",
    )

    args = parser.parse_args()
    main(
        root_dir=args.root_dir,
        fluo_dir=args.fluo_dir,
        af_dir=args.af_dir,
        output_dir=args.output_dir,
        registrar_dir=args.registrar_dir,
        sample_info_file=args.sample_info_file,
        samples=args.samples,
        workers=args.workers,
        pixel_som_grid=(args.pixel_som_x, args.pixel_som_y),
        cell_som_grid=(args.cell_som_x, args.cell_som_y),
        n_pixel_meta=args.n_pixel_meta,
        n_cell_meta=args.n_cell_meta,
        bg_method=args.bg_method,
        bg_percentile=args.bg_percentile,
        merge_overlap_threshold=args.merge_overlap_threshold,
        pixels_per_cell_residual=args.pixels_per_cell_residual,
        pixels_per_cell_som=args.pixels_per_cell_som,
        af_corr_threshold=args.af_corr_threshold,
        donor_corr_threshold=args.donor_corr_threshold,
        max_fit_pixels=args.max_fit_pixels,
        background_fit_quantile=args.background_fit_quantile,
    )
