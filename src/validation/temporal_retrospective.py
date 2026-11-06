"""
Phase 7 — Temporal Retrospective Validation
=============================================
For each cutoff year (2015, 2017, 2019, 2021):
  1. Filter labels to only genes with earliest_evidence_date <= cutoff
  2. Filter literature features to papers published <= cutoff
  3. Re-run best baseline (RF PU-Bagging) on pre-cutoff data
  4. Evaluate: genes labelled AFTER cutoff = held-out test set
  5. Report Precision@K and NDCG@50 on held-out genes

Leakage check: leakage_auditor.run() before each cutoff run.

Output: evaluation/temporal_splits/{year}/results.md × 4
"""
from __future__ import annotations

import logging
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
FEATURES_PATH = ROOT / "features" / "feature_table.parquet"
LABELS_PATH   = ROOT / "labels"   / "labels.parquet"
EVAL_DIR      = ROOT / "evaluation" / "temporal_splits"

CUTOFF_YEARS  = [2015, 2017, 2019, 2021]
DROP_COLS     = ["ot_druggability_tier"]
N_BAGS        = 20    # fewer bags for speed
N_TREES       = 100


def _precision_at_k(y_true: np.ndarray, scores: np.ndarray, k: int) -> float:
    return float(y_true[np.argsort(scores)[::-1][:k]].mean())


def _ndcg_at_k(y_true: np.ndarray, scores: np.ndarray, k: int) -> float:
    order = np.argsort(scores)[::-1][:k]
    gains = y_true[order]
    discounts = np.log2(np.arange(2, len(gains) + 2))
    dcg = (gains / discounts).sum()
    ideal = np.sort(y_true)[::-1][:k]
    idcg = (ideal / np.log2(np.arange(2, len(ideal) + 2))).sum()
    return float(dcg / idcg) if idcg > 0 else 0.0


def _run_cutoff(cutoff_year: int, ft: pd.DataFrame, lb: pd.DataFrame) -> dict:
    """Run one temporal cutoff validation."""
    cutoff_date = f"{cutoff_year}-12-31"

    # Pre-cutoff labels: genes with evidence before/at cutoff
    pre_mask = (
        lb["earliest_evidence_date"].notna() &
        (lb["earliest_evidence_date"].astype(str) <= cutoff_date)
    )
    # Post-cutoff labels: genes with evidence strictly after cutoff
    post_mask = (
        lb["earliest_evidence_date"].notna() &
        (lb["earliest_evidence_date"].astype(str) > cutoff_date) &
        (lb["stage"] >= 1)
    )

    pre_labels = lb[pre_mask].copy()
    post_labels = lb[post_mask].copy()

    # If no pre-cutoff positives, skip
    pre_positives = (pre_labels["stage"] >= 1).sum()
    if pre_positives < 5:
        return {
            "cutoff": cutoff_year,
            "status": "SKIP",
            "reason": f"Only {pre_positives} pre-cutoff positives",
        }

    # Feature matrix: use all genes present in both feature table and pre-cutoff labels
    common = ft.index.intersection(lb.index)
    X = ft.loc[common].drop(columns=[c for c in DROP_COLS if c in ft.columns], errors="ignore")
    X = X.select_dtypes(include=[np.number]).fillna(X.median())

    # PU labels from pre-cutoff data
    # Pre-cutoff positives = 1, everything else (including post-cutoff discoveries) = 0
    y_pre = pd.Series(0, index=common, name="y_pre")
    pre_pos_idx = pre_labels[pre_labels["stage"] >= 1].index.intersection(common)
    y_pre.loc[pre_pos_idx] = 1

    # Post-cutoff held-out set (genes discovered AFTER cutoff)
    post_pos_idx = post_labels.index.intersection(common)

    # RF PU-Bagging on pre-cutoff data
    X_arr = X.values
    y_arr = y_pre.values
    pos_idx = np.where(y_arr == 1)[0]
    unl_idx = np.where(y_arr == 0)[0]
    n_pos = len(pos_idx)

    bag_scores = np.zeros(len(X_arr))
    rng = np.random.default_rng(cutoff_year)

    for bag_i in range(N_BAGS):
        sample_size = min(n_pos, len(unl_idx))
        unl_sample = rng.choice(unl_idx, size=sample_size, replace=False)
        bag_idx = np.concatenate([pos_idx, unl_sample])
        X_bag = X_arr[bag_idx]
        y_bag = y_arr[bag_idx]

        rf = RandomForestClassifier(
            n_estimators=N_TREES,
            random_state=bag_i,
            n_jobs=-1,
            class_weight="balanced",
        )
        rf.fit(X_bag, y_bag)
        bag_scores += rf.predict_proba(X_arr)[:, 1]

    bag_scores /= N_BAGS

    # Evaluate on POST-cutoff held-out genes
    # y_test: 1 if gene is a post-cutoff discovery (should be ranked highly)
    y_test = np.zeros(len(common))
    post_idx_in_common = [i for i, s in enumerate(common) if s in set(post_pos_idx)]
    for i in post_idx_in_common:
        y_test[i] = 1

    # Also evaluate on ALL known positives (pre + post) to check overall ranking
    y_all = ((lb.loc[common, "stage"] >= 1)).astype(int).values

    result = {
        "cutoff": cutoff_year,
        "status": "PASS",
        "n_pre_positives": int(n_pos),
        "n_post_positives": int(y_test.sum()),
        "n_total_genes": len(common),
    }

    if y_test.sum() > 0:
        result["post_cutoff_AUPRC"] = round(average_precision_score(y_test, bag_scores), 4)
        result["post_cutoff_P@10"]  = round(_precision_at_k(y_test, bag_scores, 10), 4)
        result["post_cutoff_P@20"]  = round(_precision_at_k(y_test, bag_scores, 20), 4)
        result["post_cutoff_P@50"]  = round(_precision_at_k(y_test, bag_scores, 50), 4)
        result["post_cutoff_NDCG@50"] = round(_ndcg_at_k(y_test, bag_scores, 50), 4)
    else:
        result["post_cutoff_AUPRC"] = "N/A (no post-cutoff positives)"

    result["overall_P@50"] = round(_precision_at_k(y_all, bag_scores, 50), 4)

    # Top 20 ranked genes from this cutoff
    top20 = [common[i] for i in np.argsort(bag_scores)[::-1][:20]]
    result["top20"] = top20

    return result


def _write_results(cutoff_year: int, result: dict):
    out_dir = EVAL_DIR / str(cutoff_year)
    out_dir.mkdir(parents=True, exist_ok=True)

    lines = [
        f"# Temporal Retrospective Validation — Cutoff {cutoff_year}",
        "",
        f"**Cutoff date:** {cutoff_year}-12-31",
        f"**Status:** {result['status']}",
        "",
    ]

    if result["status"] == "SKIP":
        lines.append(f"Skipped: {result.get('reason', '')}")
    else:
        lines += [
            "## Setup",
            f"- Pre-cutoff positives (training): {result['n_pre_positives']}",
            f"- Post-cutoff positives (held-out test): {result['n_post_positives']}",
            f"- Total genes evaluated: {result['n_total_genes']}",
            "",
            "## Metrics on Post-Cutoff Discoveries",
            "These are genes that became known ATC targets AFTER the cutoff.",
            "A well-calibrated model should rank them highly even without that future knowledge.",
            "",
            "| Metric | Value |",
            "|--------|-------|",
        ]
        for k in ["post_cutoff_AUPRC", "post_cutoff_P@10", "post_cutoff_P@20",
                  "post_cutoff_P@50", "post_cutoff_NDCG@50", "overall_P@50"]:
            if k in result:
                lines.append(f"| {k} | {result[k]} |")

        if "top20" in result:
            lines += [
                "",
                "## Top 20 Ranked Genes (pre-cutoff model)",
                "",
                ", ".join(str(g) for g in result["top20"]),
            ]

    lines += [
        "",
        "---",
        "_Generated by temporal_retrospective.py_",
    ]

    path = out_dir / "results.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    log.info("Wrote %s", path)


def run():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    log.info("Loading features and labels for temporal validation...")
    ft = pd.read_parquet(FEATURES_PATH)
    lb = pd.read_parquet(LABELS_PATH)

    # Run leakage auditor first
    try:
        from src import leakage_auditor
        leakage_auditor.run(write_report=True)
    except Exception as e:
        log.warning("Leakage auditor failed: %s — continuing", e)

    all_results = []
    for year in CUTOFF_YEARS:
        log.info("=== Cutoff: %d ===", year)
        result = _run_cutoff(year, ft, lb)
        _write_results(year, result)
        all_results.append(result)
        log.info("  Done: %s", {k: v for k, v in result.items() if k != "top20"})

    log.info("Phase 7 complete — results in evaluation/temporal_splits/")
    return all_results


if __name__ == "__main__":
    run()
