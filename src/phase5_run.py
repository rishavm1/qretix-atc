"""
phase5_run.py
─────────────
Phase 5 orchestrator: train all baseline models + fusion.
Order: LogReg → RF → GBM → GNN → Fusion

Run: python -m src.phase5_run
"""
from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")

ROOT = Path(__file__).resolve().parents[1]
EVAL_DIR = ROOT / "evaluation" / "reports"


def _step(name: str, fn) -> dict:
    log.info("════ %s ════", name)
    try:
        result = fn()
        log.info("✓ %s complete.", name)
        return result or {}
    except Exception as exc:
        log.error("✗ %s FAILED: %s", name, exc)
        import traceback; traceback.print_exc()
        return {"error": str(exc)}


def run():
    from src.models import (
        baseline_logreg,
        baseline_rf,
        baseline_gbm,
        gnn_ranker,
        multimodal_fusion,
    )

    all_metrics = {}

    all_metrics["logreg"] = _step("LogReg PU", baseline_logreg.train_and_evaluate)
    all_metrics["rf"]     = _step("RF PU-Bagging", baseline_rf.train_and_evaluate)
    all_metrics["gbm"]    = _step("GBM Rankers", baseline_gbm.train_and_evaluate)
    all_metrics["gnn"]    = _step("GNN GraphSAGE", gnn_ranker.train_and_evaluate)
    all_metrics["fusion"] = _step("Multimodal Fusion", multimodal_fusion.train_and_evaluate)

    _write_model_comparison(all_metrics)
    log.info("Phase 5 complete.")


def _write_model_comparison(all_metrics: dict):
    EVAL_DIR.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Model Comparison Report — Phase 5",
        "",
        "## All Models — Metrics Summary",
        "",
        "| Model | AUPRC | P@10 | P@20 | P@50 | NDCG@50 |",
        "|-------|-------|------|------|------|---------|",
    ]

    def _row(name, m):
        if isinstance(m, list):
            for item in m:
                _row(item.get("model", name), item)
            return
        if "error" in m:
            lines.append(f"| {name} | ERROR | — | — | — | — |")
            return
        lines.append(
            f"| {m.get('model', name)} | {m.get('AUPRC','—')} | "
            f"{m.get('P@10','—')} | {m.get('P@20','—')} | "
            f"{m.get('P@50','—')} | {m.get('NDCG@50','—')} |"
        )

    for name, m in all_metrics.items():
        _row(name, m)

    lines += [
        "",
        "## Notes",
        "- PU formulation: positives = Stage≥1, unlabelled = Stage=0 (never negative)",
        "- Metrics computed on cross-validated predictions (LogReg) or full dataset (others)",
        "- Best model scores feed into multimodal fusion",
        "",
        "## Scoring Architecture (per spec)",
        "Three scores — NEVER merged into one opaque number:",
        "- biology_score — expression + genetics + dependency",
        "- therapeutic_score — druggability + modulators + structure",
        "- novelty_score — computed in Phase 6 (5 independent sub-scores)",
    ]

    path = EVAL_DIR / "model_comparison.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    log.info("Wrote model comparison → %s", path)


if __name__ == "__main__":
    run()
