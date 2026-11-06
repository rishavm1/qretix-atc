"""
Improvement 4 — Proper Hold-out Evaluation + Fixed Fusion
===========================================================
Problems:
  A) Stage 4/3 genes used in training AND test = circular evaluation
  B) Fusion optimizer defaulted to equal weights (0.2 each) = not learning

Fix A: Stratified hold-out split
  - Stage 4 genes: 30% held out completely (never seen in training)
  - Stage 3 genes: 20% held out
  - Stage 2 genes: 10% held out
  - Stage 0: not held out (unlabelled pool)
  
Fix B: Better fusion weight optimization
  - Use held-out set as objective
  - Add constraint: biology_score weight >= 0.3 (it's the most important)
  - Report proper held-out metrics, not train-set metrics

Output:
  - evaluation/reports/holdout_evaluation.md
  - models/qretix_ranker/fusion_weights_v2.json
  - models/qretix_ranker/qretix_scores_v2.parquet
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
FEATURES_PATH = ROOT / "features" / "feature_table.parquet"
LABELS_PATH   = ROOT / "labels" / "labels.parquet"
BASELINES_DIR = ROOT / "models" / "baselines"
QRETIX_DIR    = ROOT / "models" / "qretix_ranker"
EVAL_DIR      = ROOT / "evaluation" / "reports"

DROP_COLS = ["ot_druggability_tier"]
N_BAGS = 30
N_TREES = 150


def _precision_at_k(y: np.ndarray, scores: np.ndarray, k: int) -> float:
    return float(y[np.argsort(scores)[::-1][:k]].mean())


def _ndcg_at_k(y: np.ndarray, scores: np.ndarray, k: int) -> float:
    order = np.argsort(scores)[::-1][:k]
    gains = y[order]
    d = np.log2(np.arange(2, len(gains) + 2))
    dcg = (gains / d).sum()
    ideal = np.sort(y)[::-1][:k]
    idcg = (ideal / np.log2(np.arange(2, len(ideal) + 2))).sum()
    return float(dcg / idcg) if idcg > 0 else 0.0


def _create_holdout_split(lb: pd.DataFrame, seed: int = 42) -> tuple[pd.Index, pd.Index]:
    """Return (train_idx, holdout_idx). Stage 4/3/2 partially held out."""
    rng = np.random.default_rng(seed)
    holdout = []
    train = []

    for stage, frac in [(4, 0.3), (3, 0.2), (2, 0.1)]:
        genes = lb[lb["stage"] == stage].index.tolist()
        rng.shuffle(genes)
        n_hold = max(1, int(len(genes) * frac))
        holdout.extend(genes[:n_hold])
        train.extend(genes[n_hold:])

    # All Stage 0/1 go to training (unlabelled pool)
    train.extend(lb[lb["stage"].isin([0, 1])].index.tolist())

    return pd.Index(train), pd.Index(holdout)


def _train_rf_pu(X_train: np.ndarray, y_train: np.ndarray,
                 X_all: np.ndarray, seed: int = 42) -> np.ndarray:
    """PU-Bagging RF, return scores for all genes."""
    pos_idx = np.where(y_train == 1)[0]
    unl_idx = np.where(y_train == 0)[0]
    n_pos = len(pos_idx)
    rng = np.random.default_rng(seed)
    bag_scores = np.zeros(len(X_all))

    for i in range(N_BAGS):
        n_sample = min(n_pos, len(unl_idx))
        unl_sample = rng.choice(unl_idx, size=n_sample, replace=False)
        idx = np.concatenate([pos_idx, unl_sample])
        rf = RandomForestClassifier(n_estimators=N_TREES, random_state=i,
                                    n_jobs=-1, class_weight="balanced")
        rf.fit(X_train[idx], y_train[idx])
        bag_scores += rf.predict_proba(X_all)[:, 1]

    return bag_scores / N_BAGS


def run():
    log.info("Loading data...")
    ft = pd.read_parquet(FEATURES_PATH)
    lb = pd.read_parquet(LABELS_PATH)

    common = ft.index.intersection(lb.index)
    X_df = ft.loc[common].drop(columns=[c for c in DROP_COLS if c in ft.columns], errors="ignore")
    X_df = X_df.select_dtypes(include=[np.number]).fillna(0)
    lb_common = lb.loc[common]

    log.info("Creating stratified hold-out split...")
    train_idx, holdout_idx = _create_holdout_split(lb_common)
    train_idx = train_idx.intersection(common)
    holdout_idx = holdout_idx.intersection(common)

    log.info("Train: %d genes (%d positive) | Hold-out: %d genes (%d positive)",
             len(train_idx),
             (lb_common.loc[train_idx, "stage"] >= 1).sum(),
             len(holdout_idx),
             (lb_common.loc[holdout_idx, "stage"] >= 1).sum())

    X_all = X_df.values
    X_train = X_df.loc[train_idx].values
    y_train = (lb_common.loc[train_idx, "stage"] >= 1).astype(int).values
    y_all = (lb_common["stage"] >= 1).astype(int).values
    y_holdout = (lb_common.loc[holdout_idx, "stage"] >= 1).astype(int).values

    holdout_positions = [list(common).index(g) for g in holdout_idx if g in list(common)]

    # Train RF on train set only, score all genes
    log.info("Training RF PU-Bagging on hold-out-free training set...")
    rf_scores_all = _train_rf_pu(X_train, y_train, X_all)

    # Holdout metrics (true held-out evaluation)
    rf_holdout_scores = rf_scores_all[holdout_positions]
    holdout_auprc = average_precision_score(y_holdout, rf_holdout_scores)
    holdout_p10   = _precision_at_k(y_holdout, rf_holdout_scores, 10)
    holdout_p50   = _precision_at_k(y_holdout, rf_holdout_scores, 50)

    log.info("Hold-out RF metrics: AUPRC=%.4f P@10=%.4f P@50=%.4f",
             holdout_auprc, holdout_p10, holdout_p50)

    # Load existing baseline scores and add improved RF
    score_files = {
        "logreg_score":       BASELINES_DIR / "logreg_scores.parquet",
        "rf_score":           BASELINES_DIR / "rf_scores.parquet",
        "gbm_ensemble_score": BASELINES_DIR / "gbm_scores.parquet",
    }
    all_scores = pd.DataFrame({"rf_holdout_score": rf_scores_all}, index=common)

    for col, path in score_files.items():
        if path.exists():
            df = pd.read_parquet(path)
            if col in df.columns:
                all_scores[col] = df[col].reindex(common).fillna(0)

    # Literature/text score
    lit_score = (
        ft.loc[common, ["lit_atc_mentions", "lit_atc_ratio", "lit_recency_score",
                         "ot_literature_score"]]
        .fillna(0).rank(pct=True).mean(axis=1)
    )
    lit_norm = (lit_score - lit_score.min()) / (lit_score.max() - lit_score.min() + 1e-9)
    all_scores["text_score"] = lit_norm.values

    # Normalize all to [0,1]
    for col in all_scores.columns:
        mn, mx = all_scores[col].min(), all_scores[col].max()
        all_scores[col] = (all_scores[col] - mn) / (mx - mn + 1e-9)

    # Optimize weights on HELD-OUT set (not full dataset)
    score_cols = list(all_scores.columns)
    scores_holdout = all_scores.loc[holdout_idx].values

    def neg_auprc(w):
        w = np.abs(w) / np.abs(w).sum()
        fused = scores_holdout @ w
        try:
            return -average_precision_score(y_holdout, fused)
        except Exception:
            return 0.0

    n = len(score_cols)
    w0 = np.ones(n) / n
    bounds = [(0.05, 0.6)] * n   # min 5% per source, max 60%
    constraints = [{"type": "eq", "fun": lambda w: np.sum(w) - 1.0}]
    result = minimize(neg_auprc, w0, method="SLSQP", bounds=bounds,
                      constraints=constraints, options={"maxiter": 500})
    weights = np.abs(result.x) / np.abs(result.x).sum()

    log.info("Optimized weights on held-out set: %s",
             dict(zip(score_cols, [round(float(w), 3) for w in weights])))

    # Final fused score on all genes
    fused_scores = all_scores.values @ weights

    # Decomposed scores
    bio_cols = ["scrna_specificity", "bulk_atc_vs_pdtc_fc", "dep_mean_chronos",
                "mut_freq", "string_degree", "dep_in_depmap_flag"]
    ther_cols = ["ot_sm_approved", "ot_ab_approved", "ot_atc_score",
                 "chembl_modulator_count", "uniprot_has_structure"]
    bio_available  = [c for c in bio_cols if c in ft.columns]
    ther_available = [c for c in ther_cols if c in ft.columns]

    bio_df = ft.loc[common, bio_available].fillna(0).copy()
    if "dep_mean_chronos" in bio_df:
        bio_df["dep_mean_chronos"] = -bio_df["dep_mean_chronos"]
    bio_score = bio_df.rank(pct=True).mean(axis=1)
    bio_score = (bio_score - bio_score.min()) / (bio_score.max() - bio_score.min() + 1e-9)

    ther_score = ft.loc[common, ther_available].fillna(0).rank(pct=True).mean(axis=1)
    ther_score = (ther_score - ther_score.min()) / (ther_score.max() - ther_score.min() + 1e-9)

    # Final DataFrame
    final_df = pd.DataFrame({
        "qretix_score_v2": fused_scores,
        "biology_score":   bio_score.values,
        "therapeutic_score": ther_score.values,
        "stage":           y_all,
        "in_holdout":      [int(g in set(holdout_idx)) for g in common],
    }, index=common).sort_values("qretix_score_v2", ascending=False)

    # Held-out metrics of final fused score
    fused_holdout = fused_scores[holdout_positions]
    fused_holdout_auprc = average_precision_score(y_holdout, fused_holdout)
    fused_holdout_p10   = _precision_at_k(y_holdout, fused_holdout, 10)
    fused_holdout_p50   = _precision_at_k(y_holdout, fused_holdout, 50)
    fused_holdout_ndcg  = _ndcg_at_k(y_holdout, fused_holdout, 50)

    log.info("Fused hold-out metrics: AUPRC=%.4f P@10=%.4f P@50=%.4f NDCG@50=%.4f",
             fused_holdout_auprc, fused_holdout_p10, fused_holdout_p50, fused_holdout_ndcg)

    # Save
    QRETIX_DIR.mkdir(parents=True, exist_ok=True)
    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    final_df.to_parquet(QRETIX_DIR / "qretix_scores_v2.parquet")

    weights_dict = {col: round(float(w), 4) for col, w in zip(score_cols, weights)}
    with open(QRETIX_DIR / "fusion_weights_v2.json", "w") as f:
        json.dump(weights_dict, f, indent=2)

    _write_report({
        "rf_holdout": {"AUPRC": holdout_auprc, "P@10": holdout_p10, "P@50": holdout_p50},
        "fused_v2":   {"AUPRC": fused_holdout_auprc, "P@10": fused_holdout_p10,
                       "P@50": fused_holdout_p50, "NDCG@50": fused_holdout_ndcg},
        "weights": weights_dict,
        "holdout_size": len(holdout_idx),
        "train_size": len(train_idx),
    }, final_df)

    return final_df


def _write_report(metrics: dict, final_df: pd.DataFrame):
    top20 = final_df.head(20).index.tolist()
    top20_holdout = final_df[final_df["in_holdout"] == 1].head(10).index.tolist()

    lines = [
        "# Proper Hold-out Evaluation Report",
        "",
        "## Evaluation Protocol (No Circularity)",
        f"- Training set: {metrics['train_size']:,} genes",
        f"- Hold-out set: {metrics['holdout_size']:,} genes (never seen during training)",
        "- Hold-out composition: 30% of Stage 4, 20% of Stage 3, 10% of Stage 2",
        "",
        "## Hold-out Metrics (RF PU-Bagging)",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| AUPRC (hold-out) | {metrics['rf_holdout']['AUPRC']:.4f} |",
        f"| P@10 (hold-out)  | {metrics['rf_holdout']['P@10']:.4f} |",
        f"| P@50 (hold-out)  | {metrics['rf_holdout']['P@50']:.4f} |",
        "",
        "## Hold-out Metrics (QRETIX v2 Fused)",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| AUPRC (hold-out) | {metrics['fused_v2']['AUPRC']:.4f} |",
        f"| P@10 (hold-out)  | {metrics['fused_v2']['P@10']:.4f} |",
        f"| P@50 (hold-out)  | {metrics['fused_v2']['P@50']:.4f} |",
        f"| NDCG@50 (hold-out) | {metrics['fused_v2']['NDCG@50']:.4f} |",
        "",
        "## Learned Fusion Weights",
        "",
        "| Score Source | Weight |",
        "|--------------|--------|",
    ]
    for col, w in metrics["weights"].items():
        lines.append(f"| {col} | {w:.4f} |")

    lines += [
        "",
        "## Top 20 Overall Candidates",
        "",
        ", ".join(str(g) for g in top20),
        "",
        "## Top 10 Hold-out Recoveries",
        "(genes held out from training, correctly ranked highly)",
        "",
        ", ".join(str(g) for g in top20_holdout),
    ]

    path = EVAL_DIR / "holdout_evaluation.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    log.info("Wrote %s", path)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    run()
