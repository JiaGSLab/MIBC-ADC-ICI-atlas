#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phenotype filtered mIHC h5ad files using final marker positivity calls.

Main outputs:
- writes phenotype columns into each filtered h5ad
- writes phenotype columns into each filtered h5ad
- does not write CSV files by default
- optionally writes a fast combined cell table as parquet/pickle for downstream spatial analyses

Key decisions implemented:
- Skip samples without real PanCK measurement by default when --require_pancK is used.
- HER2+ is a separate tumor phenotype and has priority over PanCK+ tumor.
- PanCK+ is sufficient but not necessary for tumor.
- CD31 and SMA are not treated as conflicting lineages; CD31+SMA+ is Vascular.
- HLA-DR/HLADR alone does not define DC/APC; it is subtype/status only.
- CD3-only is not used to define T cells.
- Broad multi-positive cells are classified as undefined.
- In limited conflicts, GMM-derived positivity has priority; ties use intensity margin.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd


PANEL_GROUPS: Dict[str, List[str]] = {
    "v1_panel01": ["CD31", "CD20", "CD45", "CD21", "PanCK", "CD45RO"],
    "v1_panel01_no_PanCK": ["CD31", "CD20", "CD45", "CD21", "CD45RO"],
    "v1_panel02_T": ["CD8", "TIM3", "CD4", "PD1", "CD3", "FOXP3"],
    "v1_panel03_HER2_myeloid": ["HER2", "CD68", "CD163", "CD11c", "KI67", "HLADR"],
    "v2_panel01_myeloid_B": ["HLADR", "CD11c", "CD20", "CD163", "CD21", "CD68"],
    "v2_panel02_epi_stroma": ["CD31", "KI67", "PanCK", "HER2", "PDL1", "SMA"],
    "v2_panel03_T": ["FOXP3", "CD4", "TIM3", "CD45RO", "PD1", "CD8"],
}


def log(msg: str) -> None:
    print(msg, flush=True)


def norm_marker(x: str) -> str:
    s = str(x)
    if s.upper() in {"HLADR", "HLA_DR", "HLA-DR"}:
        return "HLADR"
    if s.upper() in {"PANK", "PANCK"}:
        return "PanCK"
    return s


def sanitize_marker(x: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_]+", "_", str(x))
    s = re.sub(r"_+", "_", s).strip("_")
    if re.match(r"^[0-9]", s):
        s = "m_" + s
    return s


def find_samples(input_root: Path, samples: Optional[List[str]]) -> List[str]:
    if samples:
        return list(samples)
    return sorted([p.name for p in input_root.iterdir() if p.is_dir()])


def add_cellid_and_spatial_columns(adata: ad.AnnData, sample: str) -> ad.AnnData:
    """
    Add globally unique cellid and canonical spatial coordinate columns.

    spatial_x/spatial_y:
        Prefer X_centroid/Y_centroid when available because these are expected
        to be coordinates after registration/reference mapping.
        Fall back to centroidX/centroidY when reference coordinates are absent.

    cellid:
        sample + "__" + existing Cell_ID/cellid/obs_name + "__x{spatial_x}__y{spatial_y}"

    If duplicates still exist, append "__dup{n}" to make cellid unique.
    """
    # 1) canonical spatial coordinates first
    if {"X_centroid", "Y_centroid"}.issubset(adata.obs.columns):
        adata.obs["spatial_x"] = pd.to_numeric(adata.obs["X_centroid"], errors="coerce")
        adata.obs["spatial_y"] = pd.to_numeric(adata.obs["Y_centroid"], errors="coerce")
        adata.uns["spatial_coordinate_source"] = "X_centroid/Y_centroid"
    elif {"centroidX", "centroidY"}.issubset(adata.obs.columns):
        adata.obs["spatial_x"] = pd.to_numeric(adata.obs["centroidX"], errors="coerce")
        adata.obs["spatial_y"] = pd.to_numeric(adata.obs["centroidY"], errors="coerce")
        adata.uns["spatial_coordinate_source"] = "centroidX/centroidY"
    else:
        adata.obs["spatial_x"] = np.nan
        adata.obs["spatial_y"] = np.nan
        adata.uns["spatial_coordinate_source"] = "not_available"

    # 2) base cell ID
    if "Cell_ID" in adata.obs.columns:
        base = adata.obs["Cell_ID"].astype(str)
    elif "cellid" in adata.obs.columns:
        base = adata.obs["cellid"].astype(str)
    else:
        base = pd.Series(adata.obs.index.astype(str), index=adata.obs.index)

    # Remove existing sample prefix to avoid repeated prefixes on rerun.
    prefix = f"{sample}__"
    base = base.str.replace(rf"^{re.escape(prefix)}", "", regex=True)

    # 3) coordinate suffix
    sx = pd.to_numeric(adata.obs["spatial_x"], errors="coerce")
    sy = pd.to_numeric(adata.obs["spatial_y"], errors="coerce")

    def _coord_to_str(v):
        if pd.isna(v):
            return "NA"
        # Coordinates are pixel-like; use rounded integer coordinates for stable IDs.
        return str(int(round(float(v))))

    x_str = sx.map(_coord_to_str).astype(str)
    y_str = sy.map(_coord_to_str).astype(str)

    adata.obs["cellid_base"] = base.values
    adata.obs["cellid"] = (
        prefix
        + base.astype(str)
        + "__x"
        + x_str
        + "__y"
        + y_str
    ).values

    # 4) guarantee uniqueness even if segmentation/registration produced duplicated coordinates.
    duplicated = adata.obs["cellid"].duplicated(keep=False)
    if duplicated.any():
        dup_rank = adata.obs.groupby("cellid", observed=False).cumcount()
        adata.obs.loc[duplicated, "cellid"] = (
            adata.obs.loc[duplicated, "cellid"].astype(str)
            + "__dup"
            + dup_rank.loc[duplicated].astype(str)
        )
        adata.uns["cellid_duplicate_suffix_added"] = True
    else:
        adata.uns["cellid_duplicate_suffix_added"] = False

    adata.uns["cellid_rule"] = (
        "cellid = sample + '__' + Cell_ID/cellid/obs_name + "
        "'__x{round(spatial_x)}__y{round(spatial_y)}'; duplicates get '__dup{n}'"
    )

    return adata

def json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        x = float(obj)
        return None if not np.isfinite(x) else x
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    try:
        if pd.isna(obj):
            return None
    except Exception:
        pass
    return str(obj)


def write_run_summary_json(records: List[Dict], out_path: Path):
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2, ensure_ascii=False, default=json_default)


def extract_obs_table_from_h5ad(h5ad_path: Path, sample: str, include_positive_cols: bool = False) -> pd.DataFrame:
    """
    Extract a lightweight obs table from a phenotyped h5ad.

    This is used only when --write_combined_table is requested.
    """
    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.copy()
    try:
        adata.file.close()
    except Exception:
        pass

    keep = [
        "cellid", "sample", "sample_ID", "panel", "image_ID", "imageid",
        "Cell_ID", "spatial_x", "spatial_y",
        "centroidX", "centroidY", "X_centroid", "Y_centroid",
        "celltype", "subcelltype", "phenotype_panel_group", "phenotype_rule",
        "phenotype_conflict", "phenotype_conflict_markers",
        "phenotype_removed_markers", "phenotype_positive_markers",
        "phenotype_n_positive_markers", "phenotype_n_positive_lineage_families",
        "phenotype_lineage_families", "phenotype_broad_multipositive",
    ]
    cols = [c for c in keep if c in obs.columns]
    if include_positive_cols:
        cols += [c for c in obs.columns if c.startswith("positive_")]

    df = obs[cols].copy()
    if "sample" not in df.columns:
        df.insert(0, "sample", sample)
    df.insert(0, "obs_name", df.index.astype(str))
    return df.reset_index(drop=True)


def write_combined_table_from_h5ads(run_records: List[Dict], output_root: Path, include_positive_cols: bool = False):
    """
    Write one combined cell table without CSV.

    Prefer parquet for speed/size. If pyarrow/fastparquet is unavailable, fall back to pickle.
    """
    ok = [r for r in run_records if r.get("status") == "ok"]
    dfs = []
    for r in ok:
        h5ad = Path(r["output_h5ad"])
        if h5ad.exists():
            dfs.append(extract_obs_table_from_h5ad(h5ad, r["sample"], include_positive_cols=include_positive_cols))
    if not dfs:
        return None

    combined = pd.concat(dfs, ignore_index=True)
    parquet_path = output_root / "all_samples_cell_phenotypes.parquet"
    try:
        combined.to_parquet(parquet_path, index=False)
        return str(parquet_path)
    except Exception:
        pkl_path = output_root / "all_samples_cell_phenotypes.pkl"
        combined.to_pickle(pkl_path)
        return str(pkl_path)


def get_marker_index(adata: ad.AnnData, marker: str) -> Optional[int]:
    marker = norm_marker(marker)
    var_names = list(map(str, adata.var_names))
    if marker in var_names:
        return var_names.index(marker)
    if marker == "HLADR":
        for alt in ["HLA-DR", "HLA_DR"]:
            if alt in var_names:
                return var_names.index(alt)
    return None


def infer_positive_col(adata: ad.AnnData, marker: str) -> Optional[str]:
    candidates = [
        f"positive_{sanitize_marker(marker)}",
        f"positive_{marker}",
        f"{marker}_positive",
    ]
    for c in candidates:
        if c in adata.obs.columns:
            return c
    return None


def is_gmm_source(source: object) -> bool:
    return str(source).startswith("gmm") if source is not None else False


def marker_source(adata: ad.AnnData, marker: str) -> str:
    j = get_marker_index(adata, marker)
    if j is None or "final_threshold_source" not in adata.var.columns:
        return ""
    return str(adata.var.iloc[j]["final_threshold_source"])


def marker_threshold(adata: ad.AnnData, marker: str) -> float:
    j = get_marker_index(adata, marker)
    if j is None or "final_threshold" not in adata.var.columns:
        return np.nan
    try:
        x = float(adata.var.iloc[j]["final_threshold"])
        return x if np.isfinite(x) else np.nan
    except Exception:
        return np.nan


def sample_has_real_pancK(adata: ad.AnnData, layer: str) -> bool:
    j = get_marker_index(adata, "PanCK")
    if j is None or layer not in adata.layers:
        return False
    vals = np.asarray(adata.layers[layer])[:, j]
    finite = np.isfinite(vals)
    if finite.sum() == 0:
        return False
    return bool(np.nanstd(vals[finite]) > 1e-8 and np.nanmax(vals[finite]) > 0)


def panel_measured_markers(adata: ad.AnnData, panel: str, layer: str, finite_frac_min: float = 0.5) -> List[str]:
    if "panel" not in adata.obs.columns:
        return []
    panel_mask = adata.obs["panel"].astype(str).to_numpy() == str(panel)
    if panel_mask.sum() == 0:
        return []
    X = np.asarray(adata.layers[layer])
    measured = []
    for j, marker in enumerate(map(str, adata.var_names)):
        if "feature_kind" in adata.var.columns and str(adata.var.iloc[j]["feature_kind"]) != "marker":
            continue
        vals = X[panel_mask, j]
        if np.isfinite(vals).mean() >= finite_frac_min:
            measured.append(norm_marker(marker))
    return sorted(set(measured))


def detect_panel_group(markers: Sequence[str]) -> str:
    s = set(map(norm_marker, markers))
    for group, members in PANEL_GROUPS.items():
        if s == set(members):
            return group
    best_group, best_score = "unknown", 0.0
    for group, members in PANEL_GROUPS.items():
        members_set = set(members)
        score = len(s & members_set) / max(len(s | members_set), 1)
        if score > best_score:
            best_group, best_score = group, score
    return best_group if best_score >= 0.80 else "unknown"


def get_positive_matrix_for_panel(adata: ad.AnnData, panel_indices: np.ndarray, markers: Sequence[str], layer: str):
    n = len(panel_indices)
    pos, margin = {}, {}
    source_map, threshold_map = {}, {}
    X = np.asarray(adata.layers[layer])

    for marker in markers:
        marker = norm_marker(marker)
        col = infer_positive_col(adata, marker)
        if col is not None:
            vals = adata.obs.iloc[panel_indices][col].astype(str).to_numpy()
            pos[marker] = vals == "positive"
        else:
            pos[marker] = np.zeros(n, dtype=bool)

        j = get_marker_index(adata, marker)
        thr = marker_threshold(adata, marker)
        threshold_map[marker] = thr
        source_map[marker] = marker_source(adata, marker)

        if j is not None and np.isfinite(thr):
            intensity = X[panel_indices, j].astype(float)
            margin[marker] = intensity - thr
        else:
            margin[marker] = np.full(n, -np.inf, dtype=float)

    idx = adata.obs.index[panel_indices]
    return pd.DataFrame(pos, index=idx), pd.DataFrame(margin, index=idx), source_map, threshold_map


def positive_markers(pos_row: pd.Series) -> List[str]:
    return [m for m in pos_row.index if bool(pos_row[m])]


def join_markers(markers: Sequence[str]) -> str:
    return ";".join([f"{m}+" for m in markers]) if markers else "undefined"


def lineages_present(pos_row: pd.Series, family_map: Dict[str, List[str]]) -> List[str]:
    out = []
    for fam, markers in family_map.items():
        if any((m in pos_row.index and bool(pos_row[m])) for m in markers):
            out.append(fam)
    return out


def family_score(pos_row: pd.Series, margin_row: pd.Series, source_map: Dict[str, str], markers: Sequence[str]):
    best_marker, best_gmm, best_margin = "", 0, -np.inf
    for m in markers:
        if m in pos_row.index and bool(pos_row[m]):
            g = 1 if is_gmm_source(source_map.get(m, "")) else 0
            mg = float(margin_row.get(m, -np.inf))
            if (g > best_gmm) or (g == best_gmm and mg > best_margin):
                best_marker, best_gmm, best_margin = m, g, mg
    return best_gmm, best_margin, best_marker


def choose_family(pos_row: pd.Series, margin_row: pd.Series, source_map: Dict[str, str], family_map: Dict[str, List[str]], families: Sequence[str]) -> str:
    scores = []
    for fam in families:
        g, mg, bm = family_score(pos_row, margin_row, source_map, family_map[fam])
        scores.append((fam, g, mg, bm))
    scores = sorted(scores, key=lambda x: (x[1], x[2]), reverse=True)
    return scores[0][0] if scores else ""


def broad_threshold(n_markers: int, custom: int = -1) -> int:
    if custom and custom > 0:
        return custom
    return 5 if n_markers >= 6 else 4


def make_result(celltype, final_pos, group, rule, conflict="none", removed=None, lineage_families=None, broad=False, raw_pos=None):
    raw_pos = list(raw_pos) if raw_pos is not None else list(final_pos)
    final_pos = list(final_pos)
    removed = list(removed) if removed is not None else []
    lineage_families = list(lineage_families) if lineage_families is not None else []
    return {
        "celltype": celltype,
        "subcelltype": join_markers(final_pos),
        "phenotype_panel_group": group,
        "phenotype_rule": rule,
        "phenotype_conflict": conflict,
        "phenotype_conflict_markers": ";".join(raw_pos),
        "phenotype_removed_markers": ";".join(removed),
        "phenotype_positive_markers": ";".join(final_pos),
        "phenotype_n_positive_markers": len(raw_pos),
        "phenotype_n_positive_lineage_families": len(lineage_families),
        "phenotype_lineage_families": ";".join(lineage_families),
        "phenotype_broad_multipositive": bool(broad),
    }


def phenotype_group1(pos_row, margin_row, source_map, broad_thr, max_families):
    group = "v1_panel01"
    raw = positive_markers(pos_row)
    family_map = {"tumor": ["PanCK"], "vascular": ["CD31"], "B": ["CD20", "CD21"]}
    fams = lineages_present(pos_row, family_map)
    if len(raw) >= broad_thr:
        return make_result("undefined", raw, group, "broad multipositive", "broad_multipositive_artifact", lineage_families=fams, broad=True, raw_pos=raw)
    if len(fams) >= max_families:
        return make_result("undefined", raw, group, ">=3 lineage families positive", "multi_lineage_conflict", lineage_families=fams, raw_pos=raw)
    final_pos, removed = list(raw), []
    if len(fams) == 2:
        winner = choose_family(pos_row, margin_row, source_map, family_map, fams)
        for fam in fams:
            if fam != winner:
                for m in family_map[fam]:
                    if m in final_pos:
                        final_pos.remove(m); removed.append(m)
        fams = [winner]
    if "PanCK" in final_pos:
        return make_result("PanCK+ Tumor", final_pos, group, "PanCK+ tumor", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "CD31" in final_pos:
        return make_result("Vascular", final_pos, group, "CD31+ vascular", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "CD20" in final_pos and "CD21" in final_pos:
        return make_result("CD21+ B cell", final_pos, group, "CD20+CD21+ B cell", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "CD20" in final_pos:
        return make_result("B cell", final_pos, group, "CD20+ B cell", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    return make_result("undefined", final_pos, group, "no defining lineage marker", "none", removed, fams, raw_pos=raw)


def phenotype_t_panel(pos_row, margin_row, source_map, group, broad_thr):
    raw = positive_markers(pos_row)
    if len(raw) >= broad_thr:
        return make_result("undefined", raw, group, "broad multipositive", "broad_multipositive_artifact", broad=True, raw_pos=raw)
    final_pos, removed = list(raw), []
    cd4, cd8, foxp3 = bool(pos_row.get("CD4", False)), bool(pos_row.get("CD8", False)), bool(pos_row.get("FOXP3", False))
    if cd4 and cd8 and not foxp3:
        family_map = {"CD4_lineage": ["CD4"], "CD8_lineage": ["CD8"]}
        winner = choose_family(pos_row, margin_row, source_map, family_map, ["CD4_lineage", "CD8_lineage"])
        if winner == "CD4_lineage":
            if "CD8" in final_pos: final_pos.remove("CD8"); removed.append("CD8")
            cd8 = False
        else:
            if "CD4" in final_pos: final_pos.remove("CD4"); removed.append("CD4")
            cd4 = False
    if cd4 and cd8 and foxp3:
        if "CD8" in final_pos: final_pos.remove("CD8"); removed.append("CD8")
        cd8 = False
    if cd4 and foxp3:
        return make_result("Treg", final_pos, group, "CD4+FOXP3+ Treg", "family_conflict_resolved" if removed else "none", removed, ["CD4_lineage"], raw_pos=raw)
    if cd4:
        return make_result("CD4 T cell", final_pos, group, "CD4+FOXP3- CD4 T cell", "family_conflict_resolved" if removed else "none", removed, ["CD4_lineage"], raw_pos=raw)
    if cd8 and (bool(pos_row.get("PD1", False)) or bool(pos_row.get("TIM3", False))):
        return make_result("CD8 exhausted T cell", final_pos, group, "CD8+PD1/TIM3+ exhausted T cell", "family_conflict_resolved" if removed else "none", removed, ["CD8_lineage"], raw_pos=raw)
    if cd8:
        return make_result("CD8 T cell", final_pos, group, "CD8+ PD1/TIM3- CD8 T cell", "family_conflict_resolved" if removed else "none", removed, ["CD8_lineage"], raw_pos=raw)
    return make_result("undefined", final_pos, group, "no CD4/CD8-defined T lineage; CD3-only not used", "none", removed, [], raw_pos=raw)


def phenotype_group3(pos_row, margin_row, source_map, broad_thr, max_families):
    group = "v1_panel03_HER2_myeloid"
    raw = positive_markers(pos_row)
    family_map = {"tumor": ["HER2"], "myeloid": ["CD68", "CD163"], "DC": ["CD11c"]}
    fams = lineages_present(pos_row, family_map)
    if len(raw) >= broad_thr:
        return make_result("undefined", raw, group, "broad multipositive", "broad_multipositive_artifact", lineage_families=fams, broad=True, raw_pos=raw)
    if len(fams) >= max_families:
        return make_result("undefined", raw, group, ">=3 lineage families positive", "multi_lineage_conflict", lineage_families=fams, raw_pos=raw)
    final_pos, removed = list(raw), []
    if "tumor" in fams and (("myeloid" in fams) or ("DC" in fams)):
        winner = choose_family(pos_row, margin_row, source_map, family_map, fams)
        for fam in fams:
            if fam != winner:
                for m in family_map[fam]:
                    if m in final_pos: final_pos.remove(m); removed.append(m)
        fams = [winner]
    if "HER2" in final_pos:
        return make_result("HER2+ Tumor", final_pos, group, "HER2+ tumor", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "CD68" in final_pos and "CD163" in final_pos:
        return make_result("CD163+ Macrophage", final_pos, group, "CD68+CD163+ macrophage", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "CD68" in final_pos:
        return make_result("CD163- Macrophage", final_pos, group, "CD68+CD163- macrophage", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "CD11c" in final_pos:
        return make_result("DC", final_pos, group, "CD11c+ DC", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    return make_result("undefined", final_pos, group, "no defining lineage marker; HLA-DR alone not used", "none", removed, fams, raw_pos=raw)


def phenotype_group4(pos_row, margin_row, source_map, broad_thr, max_families):
    group = "v2_panel01_myeloid_B"
    raw = positive_markers(pos_row)
    family_map = {"B": ["CD20", "CD21"], "myeloid": ["CD68", "CD163"], "DC": ["CD11c"]}
    fams = lineages_present(pos_row, family_map)
    if len(raw) >= broad_thr:
        return make_result("undefined", raw, group, "broad multipositive", "broad_multipositive_artifact", lineage_families=fams, broad=True, raw_pos=raw)
    if len(fams) >= max_families:
        return make_result("undefined", raw, group, ">=3 lineage families positive", "multi_lineage_conflict", lineage_families=fams, raw_pos=raw)
    final_pos, removed = list(raw), []
    if "B" in fams and (("myeloid" in fams) or ("DC" in fams)):
        winner = choose_family(pos_row, margin_row, source_map, family_map, fams)
        for fam in fams:
            if fam != winner:
                for m in family_map[fam]:
                    if m in final_pos: final_pos.remove(m); removed.append(m)
        fams = [winner]
    if "CD68" in final_pos and "CD163" in final_pos:
        return make_result("CD163+ Macrophage", final_pos, group, "CD68+CD163+ macrophage", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "CD68" in final_pos:
        return make_result("CD163- Macrophage", final_pos, group, "CD68+CD163- macrophage", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "CD11c" in final_pos:
        return make_result("DC", final_pos, group, "CD11c+ DC", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "CD20" in final_pos and "CD21" in final_pos:
        return make_result("CD21+ B cell", final_pos, group, "CD20+CD21+ B cell", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "CD20" in final_pos:
        return make_result("B cell", final_pos, group, "CD20+ B cell", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    return make_result("undefined", final_pos, group, "no defining lineage marker; HLA-DR alone not used", "none", removed, fams, raw_pos=raw)


def phenotype_group5(pos_row, margin_row, source_map, broad_thr, max_families):
    group = "v2_panel02_epi_stroma"
    raw = positive_markers(pos_row)
    family_map = {"tumor": ["PanCK", "HER2"], "vascular": ["CD31", "SMA"]}
    fams = lineages_present(pos_row, family_map)
    if len(raw) >= broad_thr:
        return make_result("undefined", raw, group, "broad multipositive", "broad_multipositive_artifact", lineage_families=fams, broad=True, raw_pos=raw)
    if len(fams) >= max_families:
        return make_result("undefined", raw, group, ">=3 lineage families positive", "multi_lineage_conflict", lineage_families=fams, raw_pos=raw)
    final_pos, removed = list(raw), []
    if "tumor" in fams and "vascular" in fams:
        winner = choose_family(pos_row, margin_row, source_map, family_map, ["tumor", "vascular"])
        loser = "vascular" if winner == "tumor" else "tumor"
        for m in family_map[loser]:
            if m in final_pos: final_pos.remove(m); removed.append(m)
        fams = [winner]
    if "HER2" in final_pos:
        return make_result("HER2+ Tumor", final_pos, group, "HER2+ tumor", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "PanCK" in final_pos:
        return make_result("PanCK+ Tumor", final_pos, group, "PanCK+ tumor", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "CD31" in final_pos:
        return make_result("Vascular", final_pos, group, "CD31+ vascular; SMA allowed", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    if "SMA" in final_pos:
        return make_result("Fibroblast / SMA+ stromal", final_pos, group, "CD31-SMA+ stromal", "family_conflict_resolved" if removed else "none", removed, fams, raw_pos=raw)
    return make_result("undefined", final_pos, group, "no defining lineage marker", "none", removed, fams, raw_pos=raw)


def phenotype_one_row(group, pos_row, margin_row, source_map, n_markers, broad_custom, max_families):
    bt = broad_threshold(n_markers, broad_custom)
    if group == "v1_panel01":
        return phenotype_group1(pos_row, margin_row, source_map, bt, max_families)
    if group == "v1_panel02_T":
        return phenotype_t_panel(pos_row, margin_row, source_map, group, bt)
    if group == "v2_panel03_T":
        return phenotype_t_panel(pos_row, margin_row, source_map, group, bt)
    if group == "v1_panel03_HER2_myeloid":
        return phenotype_group3(pos_row, margin_row, source_map, bt, max_families)
    if group == "v2_panel01_myeloid_B":
        return phenotype_group4(pos_row, margin_row, source_map, bt, max_families)
    if group == "v2_panel02_epi_stroma":
        return phenotype_group5(pos_row, margin_row, source_map, bt, max_families)
    raw = positive_markers(pos_row)
    return make_result("undefined", raw, group, "unknown or no-PanCK panel group", "unknown_panel_group", raw_pos=raw)


def phenotype_sample(h5ad_path, output_h5ad, layer, require_pancK, broad_positive_threshold, max_lineage_families):
    sample = h5ad_path.parent.name
    adata = ad.read_h5ad(h5ad_path)

    adata = add_cellid_and_spatial_columns(adata, sample)

    if require_pancK and not sample_has_real_pancK(adata, layer):
        raise RuntimeError("Sample skipped: no real PanCK measurement detected")
    if "panel" not in adata.obs.columns:
        raise KeyError("obs['panel'] not found")

    phenotype_cols = [
        "celltype", "subcelltype", "phenotype_panel_group", "phenotype_rule", "phenotype_conflict",
        "phenotype_conflict_markers", "phenotype_removed_markers", "phenotype_positive_markers",
        "phenotype_n_positive_markers", "phenotype_n_positive_lineage_families", "phenotype_lineage_families",
        "phenotype_broad_multipositive",
    ]
    for c in phenotype_cols:
        adata.obs[c] = False if c == "phenotype_broad_multipositive" else (0 if c.startswith("phenotype_n_") else "")

    panels = sorted(adata.obs["panel"].astype(str).unique())
    for panel in panels:
        markers_detected = panel_measured_markers(adata, panel, layer)
        group = detect_panel_group(markers_detected)
        if require_pancK and group == "v1_panel01_no_PanCK":
            raise RuntimeError(f"Sample skipped: panel {panel} matches no-PanCK panel")

        idx = np.where(adata.obs["panel"].astype(str).to_numpy() == str(panel))[0]
        if group in PANEL_GROUPS:
            markers = [m for m in PANEL_GROUPS[group] if get_marker_index(adata, m) is not None]
        else:
            markers = markers_detected

        pos_df, margin_df, source_map, _ = get_positive_matrix_for_panel(adata, idx, markers, layer)
        rows = []
        for obs_idx, pos_row in pos_df.iterrows():
            rows.append(
                phenotype_one_row(
                    group,
                    pos_row,
                    margin_df.loc[obs_idx],
                    source_map,
                    len(markers),
                    broad_positive_threshold,
                    max_lineage_families,
                )
            )
        res_df = pd.DataFrame(rows, index=pos_df.index)
        for c in phenotype_cols:
            adata.obs.loc[res_df.index, c] = res_df[c].values

    adata.obs["celltype"] = pd.Categorical(adata.obs["celltype"])
    adata.obs["phenotype_panel_group"] = pd.Categorical(adata.obs["phenotype_panel_group"])
    adata.obs["phenotype_conflict"] = pd.Categorical(adata.obs["phenotype_conflict"])
    adata.obs["phenotype_broad_multipositive"] = adata.obs["phenotype_broad_multipositive"].astype(bool)

    adata.uns["phenotyping_params"] = {
        "layer": layer,
        "require_pancK": bool(require_pancK),
        "broad_positive_threshold": int(broad_positive_threshold),
        "broad_positive_threshold_rule": "default: >=5 for 6-marker panels, >=4 for <=5-marker panels",
        "max_lineage_families": int(max_lineage_families),
        "CD31_SMA_conflict": False,
        "CD31_SMA_rule": "CD31+SMA+ is classified as Vascular",
        "HER2_rule": "HER2+ is HER2+ Tumor and has priority over PanCK+ Tumor",
        "CD3_only_rule": "CD3-only is not used to define T cells",
        "HLADR_only_rule": "HLA-DR/HLADR alone is not used to define DC/APC",
        "cellid_rule": "cellid = sample + '__' + Cell_ID/cellid/obs_name",
        "spatial_coordinate_source": adata.uns.get("spatial_coordinate_source", "not_available"),
        "no_csv_output": True,
    }

    adata.write_h5ad(output_h5ad)

    return {
        "sample": sample,
        "status": "ok",
        "input_h5ad": str(h5ad_path),
        "output_h5ad": str(output_h5ad),
        "n_cells": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "panels": ";".join(panels),
        "spatial_coordinate_source": adata.uns.get("spatial_coordinate_source", "not_available"),
        "message": "",
    }


def process_one_sample(sample, input_root, output_root, layer, overwrite, require_pancK, broad_positive_threshold, max_lineage_families):
    input_root, output_root = Path(input_root), Path(output_root)
    h5ad_path = input_root / sample / f"{sample}_merged_panels_filtered.h5ad"
    out_dir = output_root / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    output_h5ad = h5ad_path if overwrite else out_dir / f"{sample}_merged_panels_filtered_phenotyped.h5ad"

    if not h5ad_path.exists():
        return {
            "sample": sample,
            "status": "skip_missing_h5ad",
            "input_h5ad": str(h5ad_path),
            "output_h5ad": str(output_h5ad),
            "message": "filtered h5ad not found",
        }

    try:
        log(f"[PHENOTYPE] {sample}")
        return phenotype_sample(
            h5ad_path=h5ad_path,
            output_h5ad=output_h5ad,
            layer=layer,
            require_pancK=require_pancK,
            broad_positive_threshold=broad_positive_threshold,
            max_lineage_families=max_lineage_families,
        )
    except Exception as e:
        log(f"[ERROR] {sample}: {e}")
        traceback.print_exc()
        return {
            "sample": sample,
            "status": "error",
            "input_h5ad": str(h5ad_path),
            "output_h5ad": str(output_h5ad),
            "message": str(e),
        }


def main():
    parser = argparse.ArgumentParser(description="Phenotype filtered h5ad files without writing CSV files.")
    parser.add_argument("--input_root", required=True, type=str)
    parser.add_argument("--output_root", default=None, type=str, help="Default: same as input_root.")
    parser.add_argument("--samples", nargs="+", default=None)
    parser.add_argument("--layer", default="arcsinh_resid_mean", type=str)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--require_pancK", action="store_true")
    parser.add_argument("--allow_no_pancK", action="store_true")
    parser.add_argument("--broad_positive_threshold", type=int, default=-1)
    parser.add_argument("--max_lineage_families", type=int, default=3)

    # No CSV output. Optional combined table is written as parquet, with pickle fallback.
    parser.add_argument("--write_combined_table", action="store_true")
    parser.add_argument("--include_positive_cols_in_combined_table", action="store_true")
    parser.add_argument("--run_summary_json", default="phenotyping_run_summary.json")
    args = parser.parse_args()

    input_root = Path(args.input_root)
    output_root = Path(args.output_root) if args.output_root else input_root
    output_root.mkdir(parents=True, exist_ok=True)
    samples = find_samples(input_root, args.samples)
    require_pancK = bool(args.require_pancK and not args.allow_no_pancK)

    log(f"[INFO] input_root: {input_root}")
    log(f"[INFO] output_root: {output_root}")
    log(f"[INFO] samples: {len(samples)}")
    log(f"[INFO] layer: {args.layer}")
    log(f"[INFO] overwrite: {args.overwrite}")
    log(f"[INFO] require_pancK: {require_pancK}")
    log(f"[INFO] workers: {args.workers}")
    log("[INFO] CSV output: disabled")

    kwargs = dict(
        input_root=str(input_root),
        output_root=str(output_root),
        layer=args.layer,
        overwrite=args.overwrite,
        require_pancK=require_pancK,
        broad_positive_threshold=args.broad_positive_threshold,
        max_lineage_families=args.max_lineage_families,
    )

    results = []
    if args.workers <= 1:
        for s in samples:
            results.append(process_one_sample(sample=s, **kwargs))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(process_one_sample, sample=s, **kwargs): s for s in samples}
            for fut in as_completed(futs):
                try:
                    results.append(fut.result())
                except Exception as e:
                    results.append({"sample": futs[fut], "status": "worker_error", "message": str(e)})

    run_summary = output_root / args.run_summary_json
    write_run_summary_json(results, run_summary)

    combined_path = None
    if args.write_combined_table:
        combined_path = write_combined_table_from_h5ads(
            results,
            output_root,
            include_positive_cols=args.include_positive_cols_in_combined_table,
        )

    log("[DONE]")
    log(f"[RUN SUMMARY JSON] {run_summary}")
    if combined_path:
        log(f"[COMBINED TABLE] {combined_path}")

    status_counts = pd.Series([r.get("status", "unknown") for r in results]).value_counts(dropna=False)
    log(status_counts.to_string())


if __name__ == "__main__":
    main()
