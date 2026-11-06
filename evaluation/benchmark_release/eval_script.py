#!/usr/bin/env python3
"""
QRETIX-ATC Benchmark — Standalone Evaluator
============================================
Usage:
    python eval_script.py --predictions your_predictions.csv \
                          --ground_truth ground_truth.parquet \
                          --holdout_genes holdout_genes.txt

Prediction file format (CSV):
    hgnc_symbol,score
    BRD4,0.98
    CDK1,0.95
    ...

One row per gene. Score must be numeric (higher = more likely ATC target).
All ~19k–41k protein-coding genes should be included; missing genes receive
the minimum score.

Outputs metrics to stdout and writes results.json.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score


def precision_at_k(y, scores, k):
    return float(y[np.argsort(scores)[::-1][:k]].mean())


def ndcg_at_k(y, scores, k):
    order = np.argsort(scores)[::-1][:k]
    g = y[order].astype(float)
    d = np.log2(np.arange(2, len(g) + 2))
    dcg = (g / d).sum()
    ig = np.sort(y).astype(float)[::-1][:k]
    idcg = (ig / np.log2(np.arange(2, len(ig) + 2))).sum()
    return float(dcg / idcg) if idcg > 0 else 0.0


def evaluate(pred_path: str, gt_path: str, holdout_path: str | None) -> dict:
    preds = pd.read_csv(pred_path).set_index("hgnc_symbol")["score"]
    gt    = pd.read_parquet(gt_path)

    # If holdout gene list provided, evaluate only on those
    if holdout_path:
        holdout_genes = Path(holdout_path).read_text().splitlines()
        holdout_genes = [g.strip() for g in holdout_genes if g.strip()]
        gt = gt.loc[gt.index.intersection(holdout_genes)]

    common = gt.index.intersection(preds.index)
    missing = len(gt) - len(common)
    if missing > 0:
        print(f"WARNING: {missing} genes in ground truth have no prediction — "
              f"assigned minimum score.", file=sys.stderr)

    scores = preds.reindex(gt.index).fillna(preds.min() - 1).values
    y      = (gt["stage"] >= 2).astype(int).values

    if y.sum() == 0:
        print("ERROR: No positives in evaluation set.", file=sys.stderr)
        sys.exit(1)

    metrics = {
        "n_genes_evaluated": int(len(gt)),
        "n_positives":       int(y.sum()),
        "n_missing_preds":   int(missing),
        "AUPRC":   round(float(average_precision_score(y, scores)), 4),
        "P@10":    round(precision_at_k(y, scores, 10), 4),
        "P@20":    round(precision_at_k(y, scores, 20), 4),
        "P@50":    round(precision_at_k(y, scores, 50), 4),
        "NDCG@50": round(ndcg_at_k(y, scores, 50), 4),
    }

    print("\n=== QRETIX-ATC Benchmark Results ===")
    for k, v in metrics.items():
        print(f"  {k:<25} {v}")

    Path("results.json").write_text(json.dumps(metrics, indent=2))
    print("\nWrote results.json")
    return metrics


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions",   required=True)
    parser.add_argument("--ground_truth",  required=True)
    parser.add_argument("--holdout_genes", default=None)
    args = parser.parse_args()
    evaluate(args.predictions, args.ground_truth, args.holdout_genes)
