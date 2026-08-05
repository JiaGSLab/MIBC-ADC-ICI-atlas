#!/usr/bin/env python3
"""WP1.4 CN block adjacency/coupling and fragmentation analysis.

The geometry step follows the WP1.4 locked definition:
per-CN Delaunay -> radius200 filtered same-CN components -> block_size>=5
valid CN blocks -> all-cell Delaunay only to decide whether blocks touch.

Cross-panel/section edges are intentionally retained. The permutation null
shuffles valid CN-block labels within each sample on the same block-contact
graph, so pooled section density and cross-section geometry remain in the null.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import re
import time
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from scipy.spatial import Delaunay, QhullError


PROJECT_DIR = Path("/data/home/chenwenhui/mIHC/fig6_pm")
CELLINFO_DIR = Path("/data/home/chenwenhui/mIHC/CellInfo/BLCA_fig6_B0_refresh")
CN_DIR = CELLINFO_DIR / "neighborhood_radius200_B0"
DELIVERABLE_DIR = PROJECT_DIR / "deliverables"
SCRATCH_DIR = PROJECT_DIR / "scratch_wp1.4"

CN_TABLE = CN_DIR / "neighborhood_composition.parquet"
SAMPLE_VALIDATION = DELIVERABLE_DIR / "WP0.3_sample_validation.csv"
CN_QC = DELIVERABLE_DIR / "WP1.1_cn_qc_by_sample.csv"
NICHE_MAP_K10 = DELIVERABLE_DIR / "WP1.2_functional_niche_map_k10.csv"
NICHE_MAP_K11 = DELIVERABLE_DIR / "WP1.2_functional_niche_map_k11.csv"

OUT_ADJACENCY = DELIVERABLE_DIR / "WP1.4_adjacency_sample_metrics.csv"
OUT_FRAGMENTATION = DELIVERABLE_DIR / "WP1.4_fragmentation_coherence.csv"
OUT_STATS = DELIVERABLE_DIR / "WP1.4_CRvsNR_stats.csv"
OUT_QC = DELIVERABLE_DIR / "WP1.4_qc_blocks.csv"
OUT_REPORT = DELIVERABLE_DIR / "WP1.4_adjacency_fragmentation_report.md"

RADIUS = 200.0
MIN_BLOCK_SIZE = 5
PRIMARY_LOW_CELL_THRESHOLD = 30000
SENSITIVITY_LOW_CELL_THRESHOLD = 20000
DEFAULT_PERMUTATIONS = 200
DEFAULT_BOOTSTRAP = 5000
PRESENCE_MIN_PER_GROUP = 5

V1_LIKE_SAMPLES = {
    "pt16-pre",
    "pt17-pre",
    "pt19-pre",
    "pt20-pre",
    "pt21-pre",
    "pt23-post",
    "pt23-pre",
}
V1_LAYOUT = "v1_like_missing_SMA_CD3_PDL1"
V2_LAYOUT = "v2_like"


@dataclass(frozen=True)
class AnalysisConfig:
    name: str
    k: str
    low_cell_threshold: int | None
    exclude_v1_all_features: bool = False


ANALYSIS_CONFIGS = [
    AnalysisConfig("primary_k10_threshold30k_panelaware", "10", PRIMARY_LOW_CELL_THRESHOLD),
    AnalysisConfig("sensitivity_k10_include_lowcells_panelaware", "10", None),
    AnalysisConfig("sensitivity_k10_threshold20k_panelaware", "10", SENSITIVITY_LOW_CELL_THRESHOLD),
    AnalysisConfig("sensitivity_k10_threshold30k_exclude_v1_all", "10", PRIMARY_LOW_CELL_THRESHOLD, True),
    AnalysisConfig("sensitivity_k11_threshold30k_panelaware", "11", PRIMARY_LOW_CELL_THRESHOLD),
]


def log(message: str) -> None:
    print(message, flush=True)


def cn_sort_key(value: str) -> tuple[int, str]:
    match = re.search(r"(\d+)$", str(value))
    return (int(match.group(1)) if match else 10**9, str(value))


def pair_name(a: str, b: str) -> str:
    a, b = str(a), str(b)
    return f"{a}__{b}" if cn_sort_key(a) <= cn_sort_key(b) else f"{b}__{a}"


def normalize_cn_label(value) -> str:
    if pd.isna(value):
        return ""
    text = str(value)
    if text.startswith("CN"):
        return text
    try:
        return f"CN{int(float(text))}"
    except ValueError:
        return text


class UnionFind:
    def __init__(self, n: int):
        self.parent = np.arange(n, dtype=np.int64)
        self.rank = np.zeros(n, dtype=np.int8)

    def find(self, x: int) -> int:
        x = int(x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = int(self.parent[x])
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1


def connected_component_labels(n_nodes: int, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    uf = UnionFind(int(n_nodes))
    if edges.size:
        for a, b in np.asarray(edges, dtype=np.int64):
            uf.union(int(a), int(b))
    roots = np.array([uf.find(i) for i in range(int(n_nodes))], dtype=np.int64)
    unique_roots, labels = np.unique(roots, return_inverse=True)
    sizes = np.bincount(labels, minlength=len(unique_roots)).astype(np.int64)
    return labels.astype(np.int64), sizes


def unique_edges_from_simplices(simplices: np.ndarray) -> np.ndarray:
    if simplices.size == 0:
        return np.empty((0, 2), dtype=np.int64)
    simplices = np.asarray(simplices, dtype=np.int64)
    edges = np.concatenate(
        [simplices[:, [0, 1]], simplices[:, [0, 2]], simplices[:, [1, 2]]],
        axis=0,
    )
    edges.sort(axis=1)
    return np.unique(edges, axis=0).astype(np.int64, copy=False)


def build_delaunay_edges(coords: np.ndarray, max_edge_length: float = RADIUS, seed: int = 0) -> tuple[np.ndarray, dict]:
    coords = np.asarray(coords, dtype=np.float64)
    n = int(coords.shape[0])
    info = {
        "n_points": n,
        "delaunay_status": "not_run",
        "n_edges_before_filter": 0,
        "n_edges_after_filter": 0,
        "edge_length_median": np.nan,
        "edge_length_q95": np.nan,
        "edge_length_q99": np.nan,
        "max_edge_length": float(max_edge_length),
        "filtered_edge_fraction": np.nan,
    }
    if n < 2:
        info["delaunay_status"] = "too_few_points"
        return np.empty((0, 2), dtype=np.int64), info
    if n == 2:
        length = float(np.linalg.norm(coords[0] - coords[1]))
        keep = length <= float(max_edge_length)
        info.update(
            {
                "delaunay_status": "two_points",
                "n_edges_before_filter": 1,
                "n_edges_after_filter": int(keep),
                "edge_length_median": length,
                "edge_length_q95": length,
                "edge_length_q99": length,
                "filtered_edge_fraction": 0.0 if keep else 1.0,
            }
        )
        if keep:
            return np.array([[0, 1]], dtype=np.int64), info
        return np.empty((0, 2), dtype=np.int64), info

    rng = np.random.default_rng(seed)
    coords_j = coords + rng.normal(0.0, 1e-7, size=coords.shape)
    try:
        tri = Delaunay(coords_j, qhull_options="QJ Qbb Qc Q12")
        edges_all = unique_edges_from_simplices(tri.simplices)
        info["delaunay_status"] = "ok"
    except QhullError as exc:
        info["delaunay_status"] = f"qhull_error:{str(exc)[:150]}"
        return np.empty((0, 2), dtype=np.int64), info
    except Exception as exc:
        info["delaunay_status"] = f"error:{str(exc)[:150]}"
        return np.empty((0, 2), dtype=np.int64), info

    if edges_all.size == 0:
        info["n_edges_before_filter"] = 0
        info["n_edges_after_filter"] = 0
        return np.empty((0, 2), dtype=np.int64), info

    diffs = coords[edges_all[:, 0]] - coords[edges_all[:, 1]]
    lengths = np.sqrt(np.einsum("ij,ij->i", diffs, diffs))
    keep = lengths <= float(max_edge_length)
    edges = edges_all[keep].astype(np.int64, copy=False)
    info.update(
        {
            "n_edges_before_filter": int(len(edges_all)),
            "n_edges_after_filter": int(len(edges)),
            "edge_length_median": float(np.nanmedian(lengths)),
            "edge_length_q95": float(np.nanquantile(lengths, 0.95)),
            "edge_length_q99": float(np.nanquantile(lengths, 0.99)),
            "filtered_edge_fraction": float(1.0 - len(edges) / len(edges_all)) if len(edges_all) else np.nan,
        }
    )
    return edges, info


def build_blocks_per_cn(
    coords: np.ndarray,
    cn_labels: np.ndarray,
    global_cns: list[str],
    *,
    k: str,
    sample: str,
    min_block_size: int = MIN_BLOCK_SIZE,
    max_edge_length: float = RADIUS,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    cn_labels = np.asarray(cn_labels).astype(str)
    n = len(cn_labels)
    component_labels = np.full(n, -1, dtype=np.int64)
    valid_block_index_by_cell = np.full(n, -1, dtype=np.int64)
    block_records = []
    valid_block_cn = []
    next_component = 0
    next_valid_block = 0

    for cn in sorted(global_cns, key=cn_sort_key):
        idx = np.where(cn_labels == cn)[0]
        if len(idx) == 0:
            block_records.append(
                {
                    "sample": sample,
                    "k": k,
                    "CN": cn,
                    "block_id": -1,
                    "valid_block_index": -1,
                    "block_size": 0,
                    "retained_for_adjacency": False,
                    "block_centroid_x": np.nan,
                    "block_centroid_y": np.nan,
                    "same_cn_delaunay_status": "absent",
                    "same_cn_edges_after_filter": 0,
                    "same_cn_edges_before_filter": 0,
                    "same_cn_filtered_edge_fraction": np.nan,
                }
            )
            continue
        local_edges, edge_info = build_delaunay_edges(
            coords[idx], max_edge_length=max_edge_length, seed=seed + cn_sort_key(cn)[0] + 101
        )
        local_components, local_sizes = connected_component_labels(len(idx), local_edges)
        for local_component_id, size in enumerate(local_sizes):
            local_member = np.where(local_components == local_component_id)[0]
            global_member = idx[local_member]
            block_id = next_component
            next_component += 1
            retained = int(size) >= int(min_block_size)
            component_labels[global_member] = block_id
            valid_index = -1
            if retained:
                valid_index = next_valid_block
                next_valid_block += 1
                valid_block_cn.append(cn)
                valid_block_index_by_cell[global_member] = valid_index
            block_records.append(
                {
                    "sample": sample,
                    "k": k,
                    "CN": cn,
                    "block_id": int(block_id),
                    "valid_block_index": int(valid_index),
                    "block_size": int(size),
                    "retained_for_adjacency": bool(retained),
                    "block_centroid_x": float(coords[global_member, 0].mean()),
                    "block_centroid_y": float(coords[global_member, 1].mean()),
                    "same_cn_delaunay_status": edge_info.get("delaunay_status", ""),
                    "same_cn_edges_after_filter": int(edge_info.get("n_edges_after_filter", 0)),
                    "same_cn_edges_before_filter": int(edge_info.get("n_edges_before_filter", 0)),
                    "same_cn_filtered_edge_fraction": edge_info.get("filtered_edge_fraction", np.nan),
                }
            )

    return (
        component_labels,
        valid_block_index_by_cell,
        np.asarray(valid_block_cn, dtype=object),
        pd.DataFrame(block_records),
    )


def compute_block_contact_edges(valid_block_index_by_cell: np.ndarray, cn_labels: np.ndarray, all_edges: np.ndarray) -> np.ndarray:
    if all_edges.size == 0:
        return np.empty((0, 2), dtype=np.int64)
    block_a = valid_block_index_by_cell[all_edges[:, 0]]
    block_b = valid_block_index_by_cell[all_edges[:, 1]]
    mask = (block_a >= 0) & (block_b >= 0) & (block_a != block_b)
    if not np.any(mask):
        return np.empty((0, 2), dtype=np.int64)
    cn_labels = np.asarray(cn_labels).astype(str)
    mask &= cn_labels[all_edges[:, 0]] != cn_labels[all_edges[:, 1]]
    if not np.any(mask):
        return np.empty((0, 2), dtype=np.int64)
    out = np.column_stack([np.minimum(block_a[mask], block_b[mask]), np.maximum(block_a[mask], block_b[mask])])
    return np.unique(out.astype(np.int64, copy=False), axis=0)


def _pair_index_matrix(global_cns: list[str]) -> tuple[dict[str, int], np.ndarray]:
    global_cns = sorted(list(map(str, global_cns)), key=cn_sort_key)
    pair_names = [pair_name(a, b) for a, b in itertools.combinations(global_cns, 2)]
    pair_to_idx = {p: i for i, p in enumerate(pair_names)}
    matrix = np.full((len(global_cns), len(global_cns)), -1, dtype=np.int64)
    for i, a in enumerate(global_cns):
        for j, b in enumerate(global_cns):
            if i == j:
                continue
            matrix[i, j] = pair_to_idx[pair_name(a, b)]
    return pair_to_idx, matrix


def compute_adjacency_metrics_from_block_graph(
    sample: str,
    k: str,
    block_cn: np.ndarray,
    block_contact_edges: np.ndarray,
    global_cns: list[str],
    n_permutations: int,
    random_seed: int,
) -> pd.DataFrame:
    global_cns = sorted(list(map(str, global_cns)), key=cn_sort_key)
    pair_to_idx, pair_idx_matrix = _pair_index_matrix(global_cns)
    cn_to_code = {cn: i for i, cn in enumerate(global_cns)}
    code_to_cn = {i: cn for cn, i in cn_to_code.items()}
    pair_list = [None] * len(pair_to_idx)
    for name, idx in pair_to_idx.items():
        pair_list[idx] = name

    block_cn = np.asarray(block_cn).astype(str)
    block_codes = np.array([cn_to_code.get(cn, -1) for cn in block_cn], dtype=np.int64)
    valid_edges = np.asarray(block_contact_edges, dtype=np.int64)
    if valid_edges.size == 0:
        valid_edges = np.empty((0, 2), dtype=np.int64)

    observed_counts = np.zeros(len(pair_list), dtype=np.int64)
    if len(valid_edges):
        a = block_codes[valid_edges[:, 0]]
        b = block_codes[valid_edges[:, 1]]
        mask = (a >= 0) & (b >= 0) & (a != b)
        idx = pair_idx_matrix[a[mask], b[mask]]
        observed_counts = np.bincount(idx[idx >= 0], minlength=len(pair_list)).astype(np.int64)

    n_blocks_by_code = np.bincount(block_codes[block_codes >= 0], minlength=len(global_cns)).astype(np.int64)
    possible = np.zeros(len(pair_list), dtype=np.int64)
    for pair, idx in pair_to_idx.items():
        cn_a, cn_b = pair.split("__")
        possible[idx] = n_blocks_by_code[cn_to_code[cn_a]] * n_blocks_by_code[cn_to_code[cn_b]]
    observed_rates = np.divide(
        observed_counts.astype(float),
        possible.astype(float),
        out=np.full(len(pair_list), np.nan, dtype=float),
        where=possible > 0,
    )

    null_rates = np.full((int(n_permutations), len(pair_list)), np.nan, dtype=float)
    if n_permutations > 0 and len(valid_edges) and len(block_codes) > 1:
        rng = np.random.default_rng(random_seed)
        for perm_i in range(int(n_permutations)):
            perm_codes = rng.permutation(block_codes)
            a = perm_codes[valid_edges[:, 0]]
            b = perm_codes[valid_edges[:, 1]]
            mask = (a >= 0) & (b >= 0) & (a != b)
            idx = pair_idx_matrix[a[mask], b[mask]]
            counts = np.bincount(idx[idx >= 0], minlength=len(pair_list)).astype(float)
            null_rates[perm_i] = np.divide(
                counts,
                possible.astype(float),
                out=np.full(len(pair_list), np.nan, dtype=float),
                where=possible > 0,
            )

    null_mean = np.full(len(pair_list), np.nan, dtype=float)
    null_sd = np.full(len(pair_list), np.nan, dtype=float)
    if n_permutations > 0:
        for idx in range(len(pair_list)):
            vals = null_rates[:, idx]
            vals = vals[np.isfinite(vals)]
            if len(vals):
                null_mean[idx] = float(np.mean(vals))
            if len(vals) > 1:
                null_sd[idx] = float(np.std(vals, ddof=1))
    roe = np.divide(
        observed_rates,
        null_mean,
        out=np.full(len(pair_list), np.nan, dtype=float),
        where=np.isfinite(null_mean) & (null_mean > 0),
    )
    z_score = np.divide(
        observed_rates - null_mean,
        null_sd,
        out=np.full(len(pair_list), np.nan, dtype=float),
        where=np.isfinite(null_sd) & (null_sd > 0),
    )

    records = []
    for idx, pair in enumerate(pair_list):
        cn_i, cn_j = pair.split("__")
        records.append(
            {
                "sample": sample,
                "k": k,
                "CN_i": cn_i,
                "CN_j": cn_j,
                "CN_pair": pair,
                "observed_contact_block_pairs": int(observed_counts[idx]),
                "possible_block_pairs": int(possible[idx]),
                "pair_contact_rate": float(observed_rates[idx]) if np.isfinite(observed_rates[idx]) else np.nan,
                "null_mean_contact_rate": float(null_mean[idx]) if np.isfinite(null_mean[idx]) else np.nan,
                "null_sd_contact_rate": float(null_sd[idx]) if np.isfinite(null_sd[idx]) else np.nan,
                "RoE": float(roe[idx]) if np.isfinite(roe[idx]) else np.nan,
                "z_score": float(z_score[idx]) if np.isfinite(z_score[idx]) else np.nan,
                "n_permutations": int(n_permutations),
                "n_blocks_CN_i": int(n_blocks_by_code[cn_to_code[cn_i]]),
                "n_blocks_CN_j": int(n_blocks_by_code[cn_to_code[cn_j]]),
                "n_valid_blocks_total": int(len(block_cn)),
                "n_block_contact_edges_total": int(len(valid_edges)),
            }
        )
    return pd.DataFrame(records)


def summarize_fragmentation_for_k(
    sample: str,
    k: str,
    cn_labels: np.ndarray,
    component_labels: np.ndarray,
    all_edges: np.ndarray,
    global_cns: list[str],
    min_block_size: int = MIN_BLOCK_SIZE,
) -> pd.DataFrame:
    cn_labels = np.asarray(cn_labels).astype(str)
    component_labels = np.asarray(component_labels, dtype=np.int64)
    records = []
    for cn in sorted(global_cns, key=cn_sort_key):
        idx = np.where(cn_labels == cn)[0]
        if len(idx) == 0:
            records.append(
                {
                    "sample": sample,
                    "k": k,
                    "CN": cn,
                    "cn_cells": 0,
                    "component_count_all": 0,
                    "valid_block_count": 0,
                    "small_block_count": 0,
                    "small_block_cell_fraction": np.nan,
                    "lcc_fraction": np.nan,
                    "frag_components_per_cell": np.nan,
                    "mean_patch_size": np.nan,
                    "median_patch_size": np.nan,
                    "cn_coherence": np.nan,
                    "same_cn_edges": 0,
                    "cn_incident_edges": 0,
                }
            )
            continue
        comp_sizes = (
            pd.Series(component_labels[idx])
            .value_counts()
            .sort_index()
            .astype(int)
            .to_numpy()
        )
        same_edges = 0
        incident_edges = 0
        if all_edges.size:
            edge_a = all_edges[:, 0]
            edge_b = all_edges[:, 1]
            a_is = cn_labels[edge_a] == cn
            b_is = cn_labels[edge_b] == cn
            incident = a_is | b_is
            incident_edges = int(incident.sum())
            same_edges = int((a_is & b_is).sum())
        valid_block_count = int((comp_sizes >= int(min_block_size)).sum())
        small_block_count = int((comp_sizes < int(min_block_size)).sum())
        small_cells = int(comp_sizes[comp_sizes < int(min_block_size)].sum()) if len(comp_sizes) else 0
        records.append(
            {
                "sample": sample,
                "k": k,
                "CN": cn,
                "cn_cells": int(len(idx)),
                "component_count_all": int(len(comp_sizes)),
                "valid_block_count": valid_block_count,
                "small_block_count": small_block_count,
                "small_block_cell_fraction": small_cells / len(idx) if len(idx) else np.nan,
                "lcc_fraction": float(comp_sizes.max() / len(idx)) if len(comp_sizes) else np.nan,
                "frag_components_per_cell": float(len(comp_sizes) / len(idx)) if len(idx) else np.nan,
                "mean_patch_size": float(np.mean(comp_sizes)) if len(comp_sizes) else np.nan,
                "median_patch_size": float(np.median(comp_sizes)) if len(comp_sizes) else np.nan,
                "cn_coherence": same_edges / incident_edges if incident_edges > 0 else np.nan,
                "same_cn_edges": same_edges,
                "cn_incident_edges": incident_edges,
            }
        )
    return pd.DataFrame(records)


def read_sample_meta() -> pd.DataFrame:
    validation = pd.read_csv(SAMPLE_VALIDATION)
    qc = pd.read_csv(CN_QC)
    qc_cols = [
        "sample",
        "cn_participant_cells",
        "cn_filtered_cells",
        "cn_participant_fraction",
        "low_cell_lt30000",
        "neighbor_count_median",
        "neighbor_count_mean",
    ]
    meta = validation.merge(qc[qc_cols], on="sample", how="left", validate="one_to_one")
    meta["panel_layout"] = np.where(meta["sample"].isin(V1_LIKE_SAMPLES), V1_LAYOUT, V2_LAYOUT)
    meta["is_v1_like"] = meta["sample"].isin(V1_LIKE_SAMPLES)
    meta["passes_lowcell_30000"] = meta["cn_participant_cells"].ge(PRIMARY_LOW_CELL_THRESHOLD)
    meta["passes_lowcell_20000"] = meta["cn_participant_cells"].ge(SENSITIVITY_LOW_CELL_THRESHOLD)
    return meta


def read_niche_maps() -> tuple[pd.DataFrame, pd.DataFrame]:
    k10 = pd.read_csv(NICHE_MAP_K10)
    k11 = pd.read_csv(NICHE_MAP_K11)
    for df in [k10, k11]:
        df["CN"] = df["CN"].map(normalize_cn_label)
        df["k"] = df["k"].astype(str)
    return k10, k11


def process_one_sample(
    sample: str,
    sample_df: pd.DataFrame,
    global_cns_by_k: dict[str, list[str]],
    *,
    n_permutations: int,
    random_seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    start = time.time()
    sub = sample_df.loc[sample_df["valid_neighborhood"].astype(bool)].copy()
    sub = sub.dropna(subset=["spatial_x", "spatial_y"])
    coords = sub[["spatial_x", "spatial_y"]].to_numpy(dtype=np.float64)
    sample_seed = (random_seed + zlib.crc32(sample.encode("utf-8"))) & 0xFFFFFFFF
    all_edges, all_edge_info = build_delaunay_edges(coords, max_edge_length=RADIUS, seed=sample_seed + 999)

    adjacency_parts = []
    frag_parts = []
    qc_parts = []
    for k, cn_col in [("10", "CN_k10"), ("11", "CN_k11")]:
        cn_labels = sub[cn_col].map(normalize_cn_label).to_numpy(dtype=str)
        global_cns = global_cns_by_k[k]
        component_labels, valid_block_by_cell, valid_block_cn, block_df = build_blocks_per_cn(
            coords,
            cn_labels,
            global_cns,
            k=k,
            sample=sample,
            min_block_size=MIN_BLOCK_SIZE,
            max_edge_length=RADIUS,
            seed=sample_seed + int(k) * 1000,
        )
        block_contact_edges = compute_block_contact_edges(valid_block_by_cell, cn_labels, all_edges)
        adj = compute_adjacency_metrics_from_block_graph(
            sample=sample,
            k=k,
            block_cn=valid_block_cn,
            block_contact_edges=block_contact_edges,
            global_cns=global_cns,
            n_permutations=n_permutations,
            random_seed=sample_seed + int(k) * 2000,
        )
        frag = summarize_fragmentation_for_k(
            sample=sample,
            k=k,
            cn_labels=cn_labels,
            component_labels=component_labels,
            all_edges=all_edges,
            global_cns=global_cns,
            min_block_size=MIN_BLOCK_SIZE,
        )
        if not block_df.empty:
            qc = (
                block_df.groupby(["sample", "k", "CN"], as_index=False)
                .agg(
                    cn_cells=("block_size", "sum"),
                    block_count_all=("block_id", lambda s: int((pd.Series(s) >= 0).sum())),
                    valid_block_count=("retained_for_adjacency", "sum"),
                    small_block_count=("retained_for_adjacency", lambda s: int((~pd.Series(s).astype(bool)).sum())),
                    median_block_size=("block_size", "median"),
                    mean_block_size=("block_size", "mean"),
                    max_block_size=("block_size", "max"),
                    same_cn_edges_after_filter=("same_cn_edges_after_filter", "sum"),
                    same_cn_edges_before_filter=("same_cn_edges_before_filter", "sum"),
                )
            )
            qc["small_block_fraction"] = qc["small_block_count"] / qc["block_count_all"].replace(0, np.nan)
            qc["small_block_cell_fraction"] = qc.merge(
                frag[["sample", "k", "CN", "small_block_cell_fraction"]],
                on=["sample", "k", "CN"],
                how="left",
            )["small_block_cell_fraction"].to_numpy()
        else:
            qc = pd.DataFrame()
        for frame in [adj, frag, qc]:
            frame["n_sample_valid_cells"] = int(len(sub))
            frame["all_cell_delaunay_status"] = all_edge_info.get("delaunay_status", "")
            frame["all_cell_edges_after_filter"] = int(all_edge_info.get("n_edges_after_filter", 0))
            frame["all_cell_edges_before_filter"] = int(all_edge_info.get("n_edges_before_filter", 0))
            frame["all_cell_filtered_edge_fraction"] = all_edge_info.get("filtered_edge_fraction", np.nan)
            frame["max_edge_length"] = RADIUS
            frame["min_block_size"] = MIN_BLOCK_SIZE
        adj["n_block_contact_edges_observed_graph"] = int(len(block_contact_edges))
        adjacency_parts.append(adj)
        frag_parts.append(frag)
        qc_parts.append(qc)

    adjacency = pd.concat(adjacency_parts, ignore_index=True)
    fragmentation = pd.concat(frag_parts, ignore_index=True)
    qc_blocks = pd.concat(qc_parts, ignore_index=True)
    for frame in [adjacency, fragmentation, qc_blocks]:
        frame["sample_runtime_seconds"] = time.time() - start
    return adjacency, fragmentation, qc_blocks


def cohens_d_nr_minus_cr(cr_values: Iterable[float], nr_values: Iterable[float]) -> float:
    cr = np.asarray(list(cr_values), dtype=float)
    nr = np.asarray(list(nr_values), dtype=float)
    cr = cr[np.isfinite(cr)]
    nr = nr[np.isfinite(nr)]
    if len(cr) == 0 or len(nr) == 0:
        return np.nan
    diff = float(np.mean(nr) - np.mean(cr))
    if len(cr) + len(nr) < 3:
        return np.nan
    var_cr = np.var(cr, ddof=1) if len(cr) > 1 else 0.0
    var_nr = np.var(nr, ddof=1) if len(nr) > 1 else 0.0
    denom_df = len(cr) + len(nr) - 2
    pooled_var = ((len(cr) - 1) * var_cr + (len(nr) - 1) * var_nr) / denom_df
    if pooled_var <= 0 or not np.isfinite(pooled_var):
        return 0.0 if diff == 0 else np.nan
    return diff / math.sqrt(pooled_var)


def bootstrap_cohens_d_ci(
    cr_values: Iterable[float],
    nr_values: Iterable[float],
    *,
    n_boot: int = DEFAULT_BOOTSTRAP,
    seed: int = 0,
) -> tuple[float, float]:
    cr = np.asarray(list(cr_values), dtype=float)
    nr = np.asarray(list(nr_values), dtype=float)
    cr = cr[np.isfinite(cr)]
    nr = nr[np.isfinite(nr)]
    if len(cr) == 0 or len(nr) == 0:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    boot = np.empty(int(n_boot), dtype=float)
    for i in range(int(n_boot)):
        boot[i] = cohens_d_nr_minus_cr(
            rng.choice(cr, size=len(cr), replace=True),
            rng.choice(nr, size=len(nr), replace=True),
        )
    boot = boot[np.isfinite(boot)]
    if len(boot) == 0:
        return (np.nan, np.nan)
    return (float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975)))


def mannwhitney_p(cr_values: Iterable[float], nr_values: Iterable[float]) -> float:
    cr = np.asarray(list(cr_values), dtype=float)
    nr = np.asarray(list(nr_values), dtype=float)
    cr = cr[np.isfinite(cr)]
    nr = nr[np.isfinite(nr)]
    if len(cr) == 0 or len(nr) == 0:
        return np.nan
    try:
        return float(scipy_stats.mannwhitneyu(nr, cr, alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def bh_fdr(p_values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(p_values), dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return q
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    q[valid] = restored
    return q


def make_feature_table(
    adjacency: pd.DataFrame,
    fragmentation: pd.DataFrame,
    meta: pd.DataFrame,
    niche_maps: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    name_maps = {k: dict(zip(df["CN"], df["functional_niche"])) for k, df in niche_maps.items()}
    family_maps = {k: dict(zip(df["CN"], df["niche_family"])) for k, df in niche_maps.items()}

    feature_parts = []
    for metric in ["pair_contact_rate", "RoE", "z_score"]:
        part = adjacency[
            [
                "sample",
                "k",
                "CN_i",
                "CN_j",
                "CN_pair",
                metric,
                "observed_contact_block_pairs",
                "possible_block_pairs",
                "n_blocks_CN_i",
                "n_blocks_CN_j",
            ]
        ].copy()
        part = part.rename(columns={metric: "value", "CN_pair": "feature_id"})
        part["feature_class"] = "adjacency"
        part["metric"] = metric
        part["feature_name"] = part.apply(
            lambda r: f"{name_maps[str(r['k'])].get(r['CN_i'], r['CN_i'])} | {name_maps[str(r['k'])].get(r['CN_j'], r['CN_j'])}",
            axis=1,
        )
        part["niche_family"] = part.apply(
            lambda r: f"{family_maps[str(r['k'])].get(r['CN_i'], '')}|{family_maps[str(r['k'])].get(r['CN_j'], '')}",
            axis=1,
        )
        part["panel_confounded_direct"] = part["niche_family"].str.contains("fibroblast_stromal", regex=False)
        feature_parts.append(part)

    for metric in ["lcc_fraction", "frag_components_per_cell", "mean_patch_size", "cn_coherence"]:
        part = fragmentation[
            [
                "sample",
                "k",
                "CN",
                metric,
                "cn_cells",
                "component_count_all",
                "valid_block_count",
                "small_block_cell_fraction",
            ]
        ].copy()
        part = part.rename(columns={metric: "value", "CN": "feature_id"})
        part["feature_class"] = "fragmentation"
        part["metric"] = metric
        part["feature_name"] = part.apply(lambda r: name_maps[str(r["k"])].get(r["feature_id"], r["feature_id"]), axis=1)
        part["niche_family"] = part.apply(lambda r: family_maps[str(r["k"])].get(r["feature_id"], ""), axis=1)
        part["panel_confounded_direct"] = part["niche_family"].eq("fibroblast_stromal")
        feature_parts.append(part)

    features = pd.concat(feature_parts, ignore_index=True)
    features["k"] = features["k"].astype(str)
    features = features.merge(meta, on="sample", how="left", validate="many_to_one")
    return features


def _feature_filter_for_config(features: pd.DataFrame, config: AnalysisConfig) -> pd.Series:
    mask = features["Time_short"].eq("Pre") & features["Outcome_short"].isin(["CR", "NR"])
    mask &= features["k"].astype(str).eq(config.k)
    if config.low_cell_threshold is not None:
        mask &= features["cn_participant_cells"].ge(config.low_cell_threshold)
    return mask


def compute_crnr_stats(
    features: pd.DataFrame,
    *,
    bootstrap_n: int = DEFAULT_BOOTSTRAP,
    min_group_presence: int = PRESENCE_MIN_PER_GROUP,
    random_seed: int = 0,
) -> pd.DataFrame:
    records = []
    group_cols = [
        "feature_class",
        "metric",
        "k",
        "feature_id",
        "feature_name",
        "niche_family",
        "panel_confounded_direct",
    ]
    for config in ANALYSIS_CONFIGS:
        analysis = features.loc[_feature_filter_for_config(features, config)].copy()
        for key, sub in analysis.groupby(group_cols, dropna=False, sort=False):
            feature_class, metric, k, feature_id, feature_name, niche_family, panel_confounded = key
            feature_df = sub
            if config.exclude_v1_all_features:
                feature_df = feature_df.loc[feature_df["panel_layout"].ne(V1_LAYOUT)]
                panel_handling = "excluded_v1_all_features"
            elif bool(panel_confounded):
                feature_df = feature_df.loc[feature_df["panel_layout"].ne(V1_LAYOUT)]
                panel_handling = "excluded_v1_for_direct_panel_confounded_feature"
            else:
                panel_handling = "included_all_panel_layouts"

            cr = feature_df.loc[feature_df["Outcome_short"].eq("CR"), "value"].astype(float).to_numpy()
            nr = feature_df.loc[feature_df["Outcome_short"].eq("NR"), "value"].astype(float).to_numpy()
            cr = cr[np.isfinite(cr)]
            nr = nr[np.isfinite(nr)]
            presence_pass = len(cr) >= int(min_group_presence) and len(nr) >= int(min_group_presence)
            seed = (
                random_seed
                + zlib.crc32(f"{config.name}|{feature_class}|{metric}|{k}|{feature_id}".encode("utf-8"))
            ) & 0xFFFFFFFF
            d = cohens_d_nr_minus_cr(cr, nr)
            ci_low, ci_high = bootstrap_cohens_d_ci(cr, nr, n_boot=bootstrap_n, seed=seed)
            p_value = mannwhitney_p(cr, nr)
            if np.isfinite(d):
                direction = "higher_in_NR" if d > 0 else "higher_in_CR" if d < 0 else "no_difference"
            else:
                diff = np.nanmean(nr) - np.nanmean(cr) if len(cr) and len(nr) else np.nan
                direction = "higher_in_NR" if diff > 0 else "higher_in_CR" if diff < 0 else "no_difference"
            records.append(
                {
                    "analysis_set": config.name,
                    "low_cell_threshold": config.low_cell_threshold if config.low_cell_threshold is not None else "none",
                    "exclude_v1_all_features": bool(config.exclude_v1_all_features),
                    "panel_handling": panel_handling,
                    "feature_class": feature_class,
                    "metric": metric,
                    "k": k,
                    "feature_id": feature_id,
                    "feature_name": feature_name,
                    "niche_family": niche_family,
                    "panel_confounded_direct": bool(panel_confounded),
                    "presence_n_CR": int(len(cr)),
                    "presence_n_NR": int(len(nr)),
                    "presence_filter_min_per_group": int(min_group_presence),
                    "presence_filter_pass": bool(presence_pass),
                    "mean_CR": float(np.nanmean(cr)) if len(cr) else np.nan,
                    "mean_NR": float(np.nanmean(nr)) if len(nr) else np.nan,
                    "median_CR": float(np.nanmedian(cr)) if len(cr) else np.nan,
                    "median_NR": float(np.nanmedian(nr)) if len(nr) else np.nan,
                    "diff_NR_minus_CR": float(np.nanmean(nr) - np.nanmean(cr)) if len(cr) and len(nr) else np.nan,
                    "cohens_d_NR_minus_CR": d,
                    "cohens_d_bootstrap95_low": ci_low,
                    "cohens_d_bootstrap95_high": ci_high,
                    "wilcoxon_mannwhitney_p": p_value,
                    "direction": direction,
                }
            )
    stats_df = pd.DataFrame(records)
    stats_df["FDR_BH_within_analysis_set"] = np.nan
    for _, idx in stats_df.groupby("analysis_set").groups.items():
        idx = list(idx)
        fdr_mask = stats_df.loc[idx, "presence_filter_pass"].astype(bool).to_numpy()
        idx_pass = np.array(idx)[fdr_mask]
        if len(idx_pass):
            stats_df.loc[idx_pass, "FDR_BH_within_analysis_set"] = bh_fdr(
                stats_df.loc[idx_pass, "wilcoxon_mannwhitney_p"].to_numpy()
            )
    stats_df["passes_FDR_0.10"] = stats_df["FDR_BH_within_analysis_set"].lt(0.10)
    return stats_df.sort_values(
        ["analysis_set", "FDR_BH_within_analysis_set", "feature_class", "metric", "feature_id"],
        na_position="last",
    ).reset_index(drop=True)


def fmt_float(value: float, digits: int = 3) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def fmt_p(value: float) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    if value < 0.001:
        return f"{value:.2e}"
    return f"{value:.3f}"


def markdown_table(df: pd.DataFrame, columns: list[str], max_rows: int = 12) -> str:
    if df.empty:
        return "None."
    view = df.loc[:, columns].head(max_rows).copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            if "p" in col.lower() or "fdr" in col.lower():
                view[col] = view[col].map(fmt_p)
            else:
                view[col] = view[col].map(lambda x: fmt_float(x, 3))
    view = view.fillna("NA").astype(str)
    view = view.map(lambda x: x.replace("|", "\\|").replace("\n", " "))
    header = "| " + " | ".join(view.columns) + " |"
    sep = "| " + " | ".join(["---"] * len(view.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in view.to_numpy()]
    return "\n".join([header, sep, *rows])


def build_report(
    adjacency: pd.DataFrame,
    fragmentation: pd.DataFrame,
    stats_df: pd.DataFrame,
    qc_blocks: pd.DataFrame,
    n_permutations: int,
) -> str:
    primary = stats_df.loc[stats_df["analysis_set"].eq("primary_k10_threshold30k_panelaware")].copy()
    primary["abs_d"] = primary["cohens_d_NR_minus_CR"].abs()
    sig = primary.loc[primary["FDR_BH_within_analysis_set"].lt(0.10)].sort_values("FDR_BH_within_analysis_set")
    top_coupling = primary.loc[primary["feature_class"].eq("adjacency")].sort_values("abs_d", ascending=False)
    top_fragmentation = primary.loc[primary["feature_class"].eq("fragmentation")].sort_values("abs_d", ascending=False)
    immune_lcc = primary.loc[
        primary["metric"].eq("lcc_fraction")
        & primary["niche_family"].isin(["B_cell_lymphoid", "macrophage_immune_stromal", "CD8_T", "B_cell_immune", "CD4_Treg"])
    ].sort_values("cohens_d_NR_minus_CR")
    summary = (
        stats_df.groupby("analysis_set")
        .agg(
            n_features=("feature_id", "size"),
            n_features_in_FDR=("presence_filter_pass", "sum"),
            n_FDR_lt_0_10=("passes_FDR_0.10", "sum"),
            min_FDR=("FDR_BH_within_analysis_set", "min"),
        )
        .reset_index()
        .sort_values("analysis_set")
    )
    qc_summary = (
        qc_blocks.groupby(["k"], as_index=False)
        .agg(
            median_valid_blocks=("valid_block_count", "median"),
            median_small_block_fraction=("small_block_fraction", "median"),
            median_small_block_cell_fraction=("small_block_cell_fraction", "median"),
        )
        .sort_values("k")
    )
    lines = [
        "# WP1.4 CN Block Adjacency, Fragmentation, and CR-vs-NR Association",
        "",
        "## Method",
        "",
        "WP1.4 used the WP1.1 B0 CN assignments; it did not recluster. For each sample and k, "
        "each CN was triangulated separately with Delaunay, radius200 removed long same-CN edges, "
        "and same-CN connected components were CN blocks. Blocks with <5 cells were excluded "
        "from adjacency statistics. One all-cell Delaunay graph per sample defined whether "
        "retained different-CN blocks touched. Cross-panel/section edges were retained.",
        "",
        f"Coupling metrics are pair_contact_rate, RoE, and z_score. RoE/z_score used {n_permutations} "
        "within-sample permutations of valid CN-block labels on the same observed block-contact graph. "
        "This keeps the pooled-section geometry and density structure in the null.",
        "",
        "## QC",
        "",
        markdown_table(
            qc_summary,
            ["k", "median_valid_blocks", "median_small_block_fraction", "median_small_block_cell_fraction"],
            10,
        ),
        "",
        "## Primary k=10 CR-vs-NR Results",
        "",
    ]
    if sig.empty:
        lines.append("No primary k=10 adjacency or fragmentation/coherence feature passed BH-FDR < 0.10.")
    else:
        lines.append(f"{len(sig)} primary k10 features passed BH-FDR < 0.10:")
        lines.append(
            markdown_table(
                sig,
                [
                    "feature_class",
                    "metric",
                    "feature_id",
                    "feature_name",
                    "cohens_d_NR_minus_CR",
                    "wilcoxon_mannwhitney_p",
                    "FDR_BH_within_analysis_set",
                    "direction",
                ],
                30,
            )
        )
    lines.extend(
        [
            "",
            "Top coupling features by absolute Cohen's d:",
            "",
            markdown_table(
                top_coupling,
                [
                    "metric",
                    "feature_id",
                    "feature_name",
                    "presence_n_CR",
                    "presence_n_NR",
                    "mean_CR",
                    "mean_NR",
                    "cohens_d_NR_minus_CR",
                    "cohens_d_bootstrap95_low",
                    "cohens_d_bootstrap95_high",
                    "wilcoxon_mannwhitney_p",
                    "FDR_BH_within_analysis_set",
                    "direction",
                    "panel_handling",
                ],
                12,
            ),
            "",
            "Top fragmentation/coherence features by absolute Cohen's d:",
            "",
            markdown_table(
                top_fragmentation,
                [
                    "metric",
                    "feature_id",
                    "feature_name",
                    "presence_n_CR",
                    "presence_n_NR",
                    "mean_CR",
                    "mean_NR",
                    "cohens_d_NR_minus_CR",
                    "cohens_d_bootstrap95_low",
                    "cohens_d_bootstrap95_high",
                    "wilcoxon_mannwhitney_p",
                    "FDR_BH_within_analysis_set",
                    "direction",
                    "panel_handling",
                ],
                12,
            ),
            "",
            "Immune-niche LCC direction check:",
            "",
            markdown_table(
                immune_lcc,
                [
                    "feature_id",
                    "feature_name",
                    "mean_CR",
                    "mean_NR",
                    "cohens_d_NR_minus_CR",
                    "FDR_BH_within_analysis_set",
                    "direction",
                ],
                12,
            ),
            "",
            "## Sensitivity Summary",
            "",
            markdown_table(summary, ["analysis_set", "n_features", "n_features_in_FDR", "n_FDR_lt_0_10", "min_FDR"], 20),
            "",
            "## Outputs",
            "",
            f"- `{OUT_ADJACENCY}`",
            f"- `{OUT_FRAGMENTATION}`",
            f"- `{OUT_STATS}`",
            f"- `{OUT_QC}`",
            f"- `{OUT_REPORT}`",
            "",
        ]
    )
    return "\n".join(lines)


def write_metadata_comments(path: Path, input_path: Path, command: str) -> None:
    meta_path = path.with_suffix(path.suffix + ".metadata.json")
    payload = {
        "output": str(path),
        "input": str(input_path),
        "generated_at_unix": time.time(),
        "command": command,
    }
    meta_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="WP1.4 CN block adjacency and fragmentation analysis")
    parser.add_argument("--n-permutations", type=int, default=DEFAULT_PERMUTATIONS)
    parser.add_argument("--bootstrap-n", type=int, default=DEFAULT_BOOTSTRAP)
    parser.add_argument("--random-seed", type=int, default=14)
    parser.add_argument("--samples", nargs="*", default=None)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    DELIVERABLE_DIR.mkdir(parents=True, exist_ok=True)
    intermediate_dir = SCRATCH_DIR / "sample_intermediates"
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    meta = read_sample_meta()
    map10, map11 = read_niche_maps()
    niche_maps = {"10": map10, "11": map11}
    global_cns_by_k = {
        "10": sorted(map10["CN"].map(normalize_cn_label).tolist(), key=cn_sort_key),
        "11": sorted(map11["CN"].map(normalize_cn_label).tolist(), key=cn_sort_key),
    }

    columns = ["sample", "spatial_x", "spatial_y", "valid_neighborhood", "CN_k10", "CN_k11"]
    labels = pd.read_parquet(CN_TABLE, columns=columns)
    labels = labels.loc[labels["valid_neighborhood"].astype(bool)].copy()
    sample_order = meta["sample"].tolist()
    if args.samples:
        requested = set(args.samples)
        sample_order = [s for s in sample_order if s in requested]
        labels = labels.loc[labels["sample"].isin(sample_order)].copy()

    adjacency_parts = []
    fragmentation_parts = []
    qc_parts = []
    command = "python deliverables/WP1.4_adjacency_fragmentation.py"
    command += f" --n-permutations {args.n_permutations} --bootstrap-n {args.bootstrap_n} --random-seed {args.random_seed}"
    if args.resume:
        command += " --resume"

    for i, sample in enumerate(sample_order, start=1):
        sample_key = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample)
        adj_path = intermediate_dir / f"{sample_key}.adjacency.csv"
        frag_path = intermediate_dir / f"{sample_key}.fragmentation.csv"
        qc_path = intermediate_dir / f"{sample_key}.qc_blocks.csv"
        if args.resume and adj_path.exists() and frag_path.exists() and qc_path.exists():
            log(f"[{i}/{len(sample_order)}] reuse {sample}")
            adjacency_parts.append(pd.read_csv(adj_path))
            fragmentation_parts.append(pd.read_csv(frag_path))
            qc_parts.append(pd.read_csv(qc_path))
            continue
        log(f"[{i}/{len(sample_order)}] process {sample}")
        sample_df = labels.loc[labels["sample"].eq(sample), columns].copy()
        adjacency, fragmentation, qc_blocks = process_one_sample(
            sample,
            sample_df,
            global_cns_by_k,
            n_permutations=args.n_permutations,
            random_seed=args.random_seed,
        )
        adjacency.to_csv(adj_path, index=False)
        fragmentation.to_csv(frag_path, index=False)
        qc_blocks.to_csv(qc_path, index=False)
        adjacency_parts.append(adjacency)
        fragmentation_parts.append(fragmentation)
        qc_parts.append(qc_blocks)

    adjacency = pd.concat(adjacency_parts, ignore_index=True) if adjacency_parts else pd.DataFrame()
    fragmentation = pd.concat(fragmentation_parts, ignore_index=True) if fragmentation_parts else pd.DataFrame()
    qc_blocks = pd.concat(qc_parts, ignore_index=True) if qc_parts else pd.DataFrame()

    adjacency = adjacency.merge(meta, on="sample", how="left", validate="many_to_one")
    fragmentation = fragmentation.merge(meta, on="sample", how="left", validate="many_to_one")
    qc_blocks = qc_blocks.merge(meta, on="sample", how="left", validate="many_to_one")

    features = make_feature_table(adjacency, fragmentation, meta, niche_maps)
    stats_df = compute_crnr_stats(
        features,
        bootstrap_n=args.bootstrap_n,
        min_group_presence=PRESENCE_MIN_PER_GROUP,
        random_seed=args.random_seed,
    )
    report = build_report(adjacency, fragmentation, stats_df, qc_blocks, args.n_permutations)

    adjacency.to_csv(OUT_ADJACENCY, index=False)
    fragmentation.to_csv(OUT_FRAGMENTATION, index=False)
    stats_df.to_csv(OUT_STATS, index=False)
    qc_blocks.to_csv(OUT_QC, index=False)
    OUT_REPORT.write_text(report, encoding="utf-8")
    for path in [OUT_ADJACENCY, OUT_FRAGMENTATION, OUT_STATS, OUT_QC, OUT_REPORT]:
        write_metadata_comments(path, CN_TABLE, command)

    log(f"Wrote {OUT_ADJACENCY}")
    log(f"Wrote {OUT_FRAGMENTATION}")
    log(f"Wrote {OUT_STATS}")
    log(f"Wrote {OUT_QC}")
    log(f"Wrote {OUT_REPORT}")


if __name__ == "__main__":
    main()
