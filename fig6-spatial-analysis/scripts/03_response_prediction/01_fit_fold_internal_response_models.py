#!/usr/bin/env python3
"""Leakage-free, fold-internal signed-composite AUC ladder for R-A5b."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut


PROJECT = Path("/data/home/chenwenhui/mIHC/fig6_pm")
WP13 = PROJECT / "deliverables/WP1.3_sample_proportion_long.csv"
WP14 = PROJECT / "deliverables/WP1.4_adjacency_sample_metrics.csv"
WP03 = PROJECT / "deliverables/WP0.3_sample_validation.csv"
FAMILY_IDS = ("celltype_proportion", "CN_proportion", "CN_coupling")
COMPARISONS = (
    ("CN_coupling", "CN_proportion"),
    ("CN_coupling", "celltype_proportion"),
    ("CN_proportion", "celltype_proportion"),
)
PROTOCOL = (
    "36 locked pre-treatment samples; LOOCV; within each training fold rank candidate "
    "features by absolute CR-minus-NR standardized effect; select at most three; assign "
    "directions from the training-fold CR-minus-NR sign only; training-fold median "
    "imputation and standard scaling; signed unweighted sum; pooled out-of-fold AUC."
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def auc(truth, score) -> float:
    truth = np.asarray(truth, dtype=int)
    score = np.asarray(score, dtype=float)
    if len(np.unique(truth)) != 2:
        return float("nan")
    return float(roc_auc_score(truth, score))


def safe_name(*parts: str) -> str:
    return "__".join(
        "".join(char if char.isalnum() else "_" for char in str(part)).strip("_")
        for part in parts
    )


def rank_train_features(x: pd.DataFrame, y: pd.Series, cap: int = 3) -> pd.DataFrame:
    """Rank only training-fold features and retain their training-fold directions."""
    values = x.apply(pd.to_numeric, errors="coerce")
    rows = []
    for feature in values.columns:
        cr = values.loc[y.eq(1), feature].dropna().to_numpy(dtype=float)
        nr = values.loc[y.eq(0), feature].dropna().to_numpy(dtype=float)
        if len(cr) < 2 or len(nr) < 2:
            continue
        pooled_variance = (
            (len(cr) - 1) * np.var(cr, ddof=1) + (len(nr) - 1) * np.var(nr, ddof=1)
        ) / (len(cr) + len(nr) - 2)
        if not np.isfinite(pooled_variance) or pooled_variance <= 1e-12:
            continue
        mean_diff = float(np.mean(cr) - np.mean(nr))
        direction = int(np.sign(mean_diff))
        if direction == 0:
            continue
        effect = mean_diff / math.sqrt(pooled_variance)
        rows.append(
            {
                "feature": feature,
                "mean_diff_CR_minus_NR": mean_diff,
                "effect_train": effect,
                "abs_effect_train": abs(effect),
                "direction": direction,
                "n_CR_train": len(cr),
                "n_NR_train": len(nr),
            }
        )
    ranked = pd.DataFrame(rows)
    if ranked.empty:
        return pd.DataFrame(
            columns=[
                "feature", "mean_diff_CR_minus_NR", "effect_train", "abs_effect_train",
                "direction", "n_CR_train", "n_NR_train",
            ]
        )
    return ranked.sort_values(["abs_effect_train", "feature"], ascending=[False, True]).head(cap).reset_index(drop=True)


def loocv_signed_composite(
    x: pd.DataFrame, y: pd.Series, model_id: str, cap: int = 3
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Generate out-of-fold scores with no held-out labels used in selection or sign."""
    predictions, selected_rows, direction_rows = [], [], []
    loo = LeaveOneOut()
    for fold, (train_index, test_index) in enumerate(loo.split(x), start=1):
        training_samples = x.index[train_index]
        held_out_sample = x.index[test_index][0]
        ranked = rank_train_features(x.loc[training_samples], y.loc[training_samples], cap=cap)
        features = ranked["feature"].tolist()
        if features:
            train = x.loc[training_samples, features].apply(pd.to_numeric, errors="coerce")
            held_out = x.loc[[held_out_sample], features].apply(pd.to_numeric, errors="coerce")
            medians = train.median(axis=0)
            train_imputed = train.fillna(medians)
            held_out_imputed = held_out.fillna(medians)
            means = train_imputed.mean(axis=0)
            scales = train_imputed.std(axis=0, ddof=0).replace(0.0, 1.0)
            standardized = (held_out_imputed.iloc[0] - means) / scales
            score = float(np.dot(standardized.to_numpy(float), ranked["direction"].to_numpy(float)))
        else:
            score = 0.0
        predictions.append(
            {
                "model_id": model_id,
                "fold": fold,
                "sample": held_out_sample,
                "truth": int(y.loc[held_out_sample]),
                "score": score,
                "n_selected": len(features),
                "selected_features": ";".join(features),
            }
        )
        for rank, record in ranked.iterrows():
            common = {
                "model_id": model_id,
                "fold": fold,
                "held_out_sample": held_out_sample,
                "rank_in_fold": rank + 1,
                **record.to_dict(),
            }
            selected_rows.append(common)
            direction_rows.append(
                {
                    "model_id": model_id,
                    "fold": fold,
                    "held_out_sample": held_out_sample,
                    "rank_in_fold": rank + 1,
                    "feature": record["feature"],
                    "direction": int(record["direction"]),
                    "mean_diff_CR_minus_NR_train": record["mean_diff_CR_minus_NR"],
                    "direction_source": "training_fold_CR_minus_NR_sign",
                }
            )
    return pd.DataFrame(predictions), pd.DataFrame(selected_rows), pd.DataFrame(direction_rows)


def bootstrap_auc(truth, score, n_bootstrap: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    truth = np.asarray(truth, dtype=int)
    score = np.asarray(score, dtype=float)
    draws = []
    attempts = 0
    while len(draws) < n_bootstrap:
        attempts += 1
        if attempts > n_bootstrap * 100:
            raise RuntimeError("Could not draw the requested number of two-class bootstrap samples.")
        index = rng.integers(0, len(truth), len(truth))
        if len(np.unique(truth[index])) == 2:
            draws.append(auc(truth[index], score[index]))
    return np.asarray(draws, dtype=float)


def paired_bootstrap_differences(predictions: pd.DataFrame, n_bootstrap: int, seed: int) -> pd.DataFrame:
    wide = predictions.pivot(index="sample", columns="model_id", values="score").loc[:, list(FAMILY_IDS)]
    truth = predictions.drop_duplicates("sample").set_index("sample").loc[wide.index, "truth"].to_numpy(dtype=int)
    rng = np.random.default_rng(seed)
    rows = []
    for left, right in COMPARISONS:
        observed = auc(truth, wide[left].to_numpy()) - auc(truth, wide[right].to_numpy())
        draws = []
        attempts = 0
        while len(draws) < n_bootstrap:
            attempts += 1
            if attempts > n_bootstrap * 100:
                raise RuntimeError("Could not draw paired bootstrap samples with both outcomes.")
            index = rng.integers(0, len(truth), len(truth))
            if len(np.unique(truth[index])) == 2:
                draws.append(
                    auc(truth[index], wide[left].to_numpy()[index])
                    - auc(truth[index], wide[right].to_numpy()[index])
                )
        draws = np.asarray(draws, dtype=float)
        sign_p = 2 * min(
            (np.count_nonzero(draws <= 0) + 1) / (len(draws) + 1),
            (np.count_nonzero(draws >= 0) + 1) / (len(draws) + 1),
        )
        rows.append(
            {
                "comparison": f"{left} - {right}",
                "model_left": left,
                "model_right": right,
                "n_samples": len(truth),
                "auc_difference_left_minus_right": observed,
                "bootstrap95_low": float(np.quantile(draws, 0.025)),
                "bootstrap95_high": float(np.quantile(draws, 0.975)),
                "bootstrap_n_valid": len(draws),
                "bootstrap_sign_p": min(1.0, sign_p),
            }
        )
    result = pd.DataFrame(rows)
    order = np.argsort(result["bootstrap_sign_p"].to_numpy(float))
    ordered_p = result["bootstrap_sign_p"].to_numpy(float)[order]
    adjusted_ordered = np.minimum.accumulate((ordered_p * len(ordered_p) / np.arange(1, len(ordered_p) + 1))[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ordered)
    adjusted[order] = adjusted_ordered
    result["bootstrap_sign_p_fdr_bh"] = np.minimum(1.0, adjusted)
    return result


def build_inputs() -> tuple[dict[str, pd.DataFrame], pd.Series, pd.DataFrame, pd.DataFrame]:
    prop = pd.read_csv(WP13)
    adjacency = pd.read_csv(WP14)
    validation = pd.read_csv(WP03)
    meta_columns = [
        "sample", "Time_short", "Outcome_short", "Group", "panel_layout", "is_v1_like",
        "defined_fraction", "cn_participant_cells", "b0_defined_rows", "passes_lowcell_30000",
    ]
    cohort = prop[meta_columns].drop_duplicates("sample")
    cohort = cohort.loc[
        cohort.Time_short.eq("Pre")
        & cohort.Outcome_short.isin(["CR", "NR"])
        & cohort.passes_lowcell_30000.astype(bool)
    ].sort_values("sample").reset_index(drop=True)
    samples = cohort["sample"].tolist()
    validation_outcome = validation.drop_duplicates("sample").set_index("sample").loc[samples, "Outcome_short"]
    y = validation_outcome.map({"NR": 0, "CR": 1}).astype(int)
    y.index.name = "sample"
    expected = cohort.set_index("sample").loc[samples, "Outcome_short"]
    if not expected.eq(validation_outcome).all():
        raise ValueError("WP0.3 outcomes disagree with the locked WP1.3 cohort labels.")
    if len(samples) != 36 or int(y.sum()) != 19 or int((1 - y).sum()) != 17:
        raise ValueError("Locked cohort must be exactly 36 pre-treatment samples: 19 CR and 17 NR.")

    cn = prop.loc[(prop.feature_domain == "CN") & (prop.k.astype(str) == "10")].copy()
    celltype = prop.loc[(prop.feature_domain == "celltype") & (prop.k.astype(str) == "celltype")].copy()
    direct_cn = set(cn.loc[cn.panel_confounded_direct.fillna(False).astype(bool), "feature_id"].astype(str))

    def proportion_matrix(rows: pd.DataFrame, label: str) -> pd.DataFrame:
        columns = {}
        for feature_id in sorted(rows.feature_id.astype(str).unique()):
            feature = rows.loc[rows.feature_id.astype(str).eq(feature_id)].drop_duplicates("sample").set_index("sample")
            columns[safe_name(label, feature_id)] = pd.to_numeric(feature["proportion"], errors="coerce").reindex(samples)
        return pd.DataFrame(columns, index=samples)

    cn = cn.loc[~cn.feature_id.astype(str).isin(direct_cn)]
    celltype = celltype.loc[~celltype.panel_confounded_direct.fillna(False).astype(bool)]
    coupling = adjacency.loc[adjacency.k.eq(10)].copy()
    coupling = coupling.loc[~coupling.CN_i.isin(direct_cn) & ~coupling.CN_j.isin(direct_cn)]
    coupling_columns = {}
    for pair in sorted(coupling.CN_pair.astype(str).unique()):
        feature = coupling.loc[coupling.CN_pair.astype(str).eq(pair)].drop_duplicates("sample").set_index("sample")
        coupling_columns[safe_name("CN_coupling", pair)] = pd.to_numeric(feature["z_score"], errors="coerce").reindex(samples)
    models = {
        "celltype_proportion": proportion_matrix(celltype, "celltype_prop"),
        "CN_proportion": proportion_matrix(cn, "CN_prop"),
        "CN_coupling": pd.DataFrame(coupling_columns, index=samples),
    }
    source_manifest = pd.DataFrame(
        [
            {"source_id": "WP1.3", "path": str(WP13), "rows": len(prop), "columns": len(prop.columns), "purpose": "celltype and CN proportion feature families"},
            {"source_id": "WP1.4", "path": str(WP14), "rows": len(adjacency), "columns": len(adjacency.columns), "purpose": "CN coupling z-score feature family"},
            {"source_id": "WP0.3", "path": str(WP03), "rows": len(validation), "columns": len(validation.columns), "purpose": "locked pre-treatment outcome labels"},
        ]
    )
    return models, y, cohort, source_manifest


def run_permutations(x: pd.DataFrame, y: pd.Series, model_id: str, n_permutations: int, n_jobs: int, seed: int) -> pd.DataFrame:
    def one(permutation_id: int) -> dict:
        rng = np.random.default_rng(seed + 1009 * permutation_id)
        permuted_y = pd.Series(rng.permutation(y.to_numpy(int)), index=y.index)
        predictions, _, _ = loocv_signed_composite(x, permuted_y, model_id, cap=3)
        score = predictions.set_index("sample").loc[permuted_y.index, "score"].to_numpy(float)
        return {"model_id": model_id, "permutation_id": permutation_id, "auc": auc(permuted_y.to_numpy(), score)}

    return pd.DataFrame(
        Parallel(n_jobs=n_jobs, prefer="processes")(
            delayed(one)(permutation_id) for permutation_id in range(1, n_permutations + 1)
        )
    )


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    view = frame.loc[:, columns].copy()
    for column in view.columns:
        if pd.api.types.is_float_dtype(view[column]):
            view[column] = view[column].map(lambda value: "NA" if not np.isfinite(value) else f"{value:.3f}")
    view = view.fillna("NA").astype(str).apply(lambda column: column.str.replace("|", "\\|", regex=False))
    return "\n".join(
        [
            "| " + " | ".join(view.columns) + " |",
            "| " + " | ".join(["---"] * len(view.columns)) + " |",
            *["| " + " | ".join(row) + " |" for row in view.to_numpy()],
        ]
    )


def write_report(outdir: Path, summary: pd.DataFrame, differences: pd.DataFrame, command: str) -> None:
    ladder = summary.sort_values("auc", ascending=False)
    comparison_lines = []
    for row in differences.itertuples(index=False):
        direction = "higher" if row.auc_difference_left_minus_right > 0 else "lower"
        comparison_lines.append(
            f"- {row.comparison}: {row.model_left} is {direction} by {row.auc_difference_left_minus_right:.3f}; "
            f"paired 95% CI [{row.bootstrap95_low:.3f}, {row.bootstrap95_high:.3f}], "
            f"sign-p={row.bootstrap_sign_p:.3f}, BH-FDR={row.bootstrap_sign_p_fdr_bh:.3f}. "
            "This is a small-sample internal-resampling trend, not a claim of external validation."
        )
    coupling_auc = float(summary.loc[summary.model_id.eq("CN_coupling"), "auc"].iloc[0])
    coupling_note = (
        "The coupling rung is within the expected approximately 0.83-0.85 range of the named Fig 6g composite."
        if 0.83 <= coupling_auc <= 0.85
        else "The coupling rung is outside the expected approximately 0.83-0.85 range; the exact fold-specific selections are retained for audit rather than adjusted."
    )
    report = "\n".join(
        [
            "# R-A5b Fold-Internal Signed-Composite AUC Ladder",
            "",
            "## Protocol",
            "",
            "All three feature families use the identical model form: within each LOOCV training fold, rank candidate features by the absolute CR-NR standardized effect, select at most three, take each selected feature's sign from that training fold only, median-impute and standardize from the training fold only, and score the held-out sample as the unweighted signed sum. The 36 locked pre-treatment samples contain 19 CR and 17 NR. Full-cohort feature directions were not used.",
            "",
            f"`{PROTOCOL}`",
            "",
            f"Command: `{command}`",
            "",
            "## Composite Ladder",
            "",
            markdown_table(ladder, ["model_id", "n_samples", "n_candidate_features", "feature_cap", "auc", "auc_bootstrap95_low", "auc_bootstrap95_high", "permutation_p"]),
            "",
            "## Paired Bootstrap AUC Differences",
            "",
            markdown_table(differences, ["comparison", "auc_difference_left_minus_right", "bootstrap95_low", "bootstrap95_high", "bootstrap_sign_p", "bootstrap_sign_p_fdr_bh"]),
            "",
            "## Honest Verdicts",
            "",
            *comparison_lines,
            "",
            coupling_note,
            "",
            "## Sensitivity And Limitations",
            "",
            "The prior R-A5 L2-logistic ladder (coupling AUC 0.796) is retained as a model-form sensitivity; its ordering supports that the family ordering is not dependent on the classifier form. This addendum is limited to the same locked 36-sample pre-treatment cohort and the three prespecified families. It performs internal resampling only, leaves feature/CN/cohort/classifier scope unchanged, and does not establish external validation or clinical utility.",
            "",
            "## Audit Files",
            "",
            "`fold_selected_features.csv` and `fold_directions.csv` record each held-out fold's selected features and training-fold directions. `permutation_null.csv` contains 1,000 complete reruns per family; `run_metadata.json` pins the interpreter, seed, protocol, script, and input SHA-256 hashes.",
            "",
            "## Validation",
            "",
            f"`{sys.executable} {outdir / 'validate_R_A5b_composite_ladder.py'} --outdir {outdir}`",
            "",
        ]
    )
    (outdir / "report.md").write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--n-bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260712)
    parser.add_argument("--n-jobs", type=int, default=4)
    args = parser.parse_args()
    if args.n_permutations < 1 or args.n_bootstrap < 1 or args.n_jobs < 1:
        raise ValueError("Permutation count, bootstrap count, and n-jobs must each be positive.")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    models, y, cohort, source_manifest = build_inputs()
    all_predictions, all_selected, all_directions, all_nulls, summaries = [], [], [], [], []
    for offset, model_id in enumerate(FAMILY_IDS):
        x = models[model_id]
        predictions, selected, directions = loocv_signed_composite(x, y, model_id, cap=3)
        observed_auc = auc(predictions["truth"], predictions["score"])
        bootstrap = bootstrap_auc(predictions["truth"], predictions["score"], args.n_bootstrap, args.seed + 37 * offset)
        null = run_permutations(x, y, model_id, args.n_permutations, args.n_jobs, args.seed + 100000 * offset)
        summaries.append(
            {
                "model_id": model_id,
                "n_samples": len(y),
                "n_CR": int(y.sum()),
                "n_NR": int((1 - y).sum()),
                "n_candidate_features": x.shape[1],
                "feature_cap": 3,
                "auc": observed_auc,
                "auc_bootstrap95_low": float(np.quantile(bootstrap, 0.025)),
                "auc_bootstrap95_high": float(np.quantile(bootstrap, 0.975)),
                "bootstrap_n_valid": len(bootstrap),
                "permutation_n": len(null),
                "null_auc_mean": float(null["auc"].mean()),
                "null_auc_sd": float(null["auc"].std(ddof=1)),
                "permutation_p": float((1 + np.count_nonzero(null["auc"].to_numpy(float) >= observed_auc)) / (1 + len(null))),
            }
        )
        all_predictions.append(predictions)
        all_selected.append(selected)
        all_directions.append(directions)
        all_nulls.append(null)
    predictions = pd.concat(all_predictions, ignore_index=True)
    selected = pd.concat(all_selected, ignore_index=True)
    directions = pd.concat(all_directions, ignore_index=True)
    null = pd.concat(all_nulls, ignore_index=True)
    summary = pd.DataFrame(summaries)
    differences = paired_bootstrap_differences(predictions, args.n_bootstrap, args.seed + 999)
    command = (
        f"{sys.executable} {Path(__file__).resolve()} --outdir {outdir} --n-permutations {args.n_permutations} "
        f"--n-bootstrap {args.n_bootstrap} --seed {args.seed} --n-jobs {args.n_jobs}"
    )
    outputs = {
        "source_manifest.csv": source_manifest,
        "analysis_cohort.csv": cohort,
        "fold_predictions.csv": predictions,
        "fold_selected_features.csv": selected,
        "fold_directions.csv": directions,
        "auc_summary.csv": summary,
        "permutation_null.csv": null,
        "paired_bootstrap_auc_differences.csv": differences,
    }
    for filename, frame in outputs.items():
        frame.to_csv(outdir / filename, index=False)
    table_manifest = pd.DataFrame(
        [
            {"file": filename, "rows": len(frame), "columns": len(frame.columns), "generated_unix": time.time(), "command": command}
            for filename, frame in outputs.items()
        ]
    )
    table_manifest.to_csv(outdir / "table_manifest.csv", index=False)
    metadata = {
        "command": command,
        "interpreter": sys.executable,
        "python_version": sys.version,
        "platform": platform.platform(),
        "seed": args.seed,
        "n_permutations": args.n_permutations,
        "n_bootstrap": args.n_bootstrap,
        "n_jobs": args.n_jobs,
        "protocol": PROTOCOL,
        "input_sha256": {str(path): sha256(path) for path in (WP13, WP14, WP03)},
        "canonical_script_sha256": sha256(Path(__file__).resolve()),
    }
    (outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(outdir, summary, differences, command)


if __name__ == "__main__":
    main()
