"""
Phase 9 — Benchmarking Suite
==============================
Compares QRETIX against 4 baseline methods on the same gene universe.

Methods:
  1. DE-only          — rank by bulk ATC vs PDTC fold-change alone
  2. DepMap-only       — rank by mean Chronos score alone (inverted)
  3. Open Targets     — rank by OT ATC association score alone
  4. Network centrality — rank by STRING PPI degree alone
  5. QRETIX           — multimodal fusion score

Evaluation: Precision@K (10, 20, 50), AUPRC, NDCG@50
Ground truth: Stage >= 2 genes (strong evidence labels)

Output: evaluation/reports/benchmark_comparison.md
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
FEATURES_PATH = ROOT / "features" / "feature_table.parquet"
LABELS_PATH   = ROOT / "labels" / "labels.parquet"
QRETIX_PATH   = ROOT / "models" / "qretix_ranker" / "qretix_scores.parquet"
EVAL_DIR      = ROOT / "evaluation" / "reports"


def _precision_at_k(y: np.ndarray, scores: np.ndarray, k: int) -> float:
    return float(y[np.argsort(scores)[::-1][:k]].mean())


def _ndcg_at_k(y: np.ndarray, scores: np.ndarray, k: int) -> float:
    order = np.argsort(scores)[::-1][:k]
    gains = y[order]
    discounts = np.log2(np.arange(2, len(gains) + 2))
    dcg = (gains / discounts).sum()
    ideal = np.sort(y)[::-1][:k]
    idcg = (ideal / np.log2(np.arange(2, len(ideal) + 2))).sum()
    return float(dcg / idcg) if idcg > 0 else 0.0


def _metrics(name: str, y: np.ndarray, scores: np.ndarray) -> dict:
    # Handle NaN scores
    nan_mask = np.isnan(scores)
    scores = scores.copy()
    scores[nan_mask] = scores[~nan_mask].min() - 1 if (~nan_mask).any() else 0.0
    return {
        "method": name,
        "AUPRC":    round(average_precision_score(y, scores), 4),
        "P@10":     round(_precision_at_k(y, scores, 10), 4),
        "P@20":     round(_precision_at_k(y, scores, 20), 4),
        "P@50":     round(_precision_at_k(y, scores, 50), 4),
        "NDCG@50":  round(_ndcg_at_k(y, scores, 50), 4),
    }


def run() -> list[dict]:
    log.info("Loading data for benchmarking...")
    ft = pd.read_parquet(FEATURES_PATH)
    lb = pd.read_parquet(LABELS_PATH)

    common = ft.index.intersection(lb.index)
    ft = ft.loc[common]
    lb = lb.loc[common]

    # Ground truth: Stage >= 2 (strong evidence only)
    y = (lb["stage"] >= 2).astype(int).values
    log.info("Benchmark ground truth: %d Stage>=2 positives / %d total genes",
             y.sum(), len(y))

    all_metrics = []

    # 1. DE-only: rank by |ATC vs PDTC fold-change|
    de_scores = ft["bulk_atc_vs_pdtc_fc"].fillna(0).abs().values
    all_metrics.append(_metrics("DE-only", y, de_scores))

    # 2. DepMap-only: rank by inverted Chronos (lower = more essential = higher rank)
    dep_scores = (-ft["dep_mean_chronos"]).fillna(0).values
    all_metrics.append(_metrics("DepMap-only", y, dep_scores))

    # 3. Open Targets
    ot_scores = ft["ot_atc_score"].fillna(0).values
    all_metrics.append(_metrics("Open Targets", y, ot_scores))

    # 4. Network centrality
    net_scores = ft["string_degree"].fillna(0).values
    all_metrics.append(_metrics("Network centrality", y, net_scores))

    # 5. QRETIX
    if QRETIX_PATH.exists():
        qretix_df = pd.read_parquet(QRETIX_PATH)
        qretix_scores = (
            qretix_df["qretix_score"]
            .reindex(common)
            .fillna(0)
            .values
        )
        all_metrics.append(_metrics("QRETIX", y, qretix_scores))
    else:
        log.warning("QRETIX scores not found — using placeholder")
        all_metrics.append({"method": "QRETIX", "AUPRC": "N/A",
                             "P@10": "N/A", "P@20": "N/A",
                             "P@50": "N/A", "NDCG@50": "N/A"})

    # Print summary
    log.info("\nBenchmark results:")
    for m in all_metrics:
        log.info("  %s: AUPRC=%.4f P@50=%.4f",
                 m["method"],
                 float(m["AUPRC"]) if m["AUPRC"] != "N/A" else 0,
                 float(m["P@50"]) if m["P@50"] != "N/A" else 0)

    _write_report(all_metrics, y, ft, common)
    return all_metrics


def _write_report(all_metrics: list[dict], y: np.ndarray,
                  ft: pd.DataFrame, common: pd.Index):
    EVAL_DIR.mkdir(parents=True, exist_ok=True)

    # Top 20 genes per method
    method_top20 = {}
    for m in all_metrics:
        name = m["method"]
        score_map = {
            "DE-only": ft["bulk_atc_vs_pdtc_fc"].fillna(0).abs(),
            "DepMap-only": (-ft["dep_mean_chronos"]).fillna(0),
            "Open Targets": ft["ot_atc_score"].fillna(0),
            "Network centrality": ft["string_degree"].fillna(0),
        }
        if name in score_map:
            top = score_map[name].nlargest(20).index.tolist()
            method_top20[name] = top
        elif name == "QRETIX" and QRETIX_PATH.exists():
            q = pd.read_parquet(QRETIX_PATH)
            method_top20[name] = q.head(20).index.tolist()

    lines = [
        "# Benchmark Comparison Report — Phase 9",
        "",
        "## Evaluation Setup",
        f"- Ground truth: Stage ≥ 2 genes ({y.sum()} positives / {len(y)} total)",
        "- Methods compared: DE-only, DepMap-only, Open Targets, Network centrality, QRETIX",
        "- Metrics: AUPRC, Precision@K, NDCG@50",
        "",
        "## Results",
        "",
        "| Method | AUPRC | P@10 | P@20 | P@50 | NDCG@50 |",
        "|--------|-------|------|------|------|---------|",
    ]
    for m in all_metrics:
        lines.append(
            f"| {m['method']} | {m['AUPRC']} | {m['P@10']} | "
            f"{m['P@20']} | {m['P@50']} | {m['NDCG@50']} |"
        )

    lines += ["", "## Top 20 Genes per Method", ""]
    for name, genes in method_top20.items():
        lines.append(f"### {name}")
        lines.append(", ".join(str(g) for g in genes))
        lines.append("")

    lines += [
        "---",
        "_Generated by benchmark_suite.py_",
    ]

    path = EVAL_DIR / "benchmark_comparison.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    log.info("Wrote %s", path)


def export_benchmark_package() -> None:
    """Generate the Q3DC submission package after benchmarking completes."""
    from src.benchmarking.benchmark_export import run as export_run
    export_run()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    run()
    export_benchmark_package()
