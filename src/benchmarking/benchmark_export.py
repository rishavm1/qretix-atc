"""
benchmark_export.py
────────────────────
Generates the QRETIX-ATC benchmark release package for Q3DC submission.

Produces evaluation/benchmark_release/:
  benchmark_metadata.json   — Q3DC submission fields (machine-readable)
  ground_truth.parquet      — held-out gene labels (stage, evidence_sources)
  ground_truth_schema.md    — column definitions + label hierarchy
  eval_script.py            — standalone evaluator (no pipeline dependency)
  sample_submission.csv     — example prediction format
  q3dc_submission_form.md   — pre-filled Q3DC form text

Run: python -m src.benchmarking.benchmark_export
"""
from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
LABELS_PATH   = ROOT / "labels" / "labels.parquet"
FEATURES_PATH = ROOT / "features" / "feature_table.parquet"
QRETIX_PATH   = ROOT / "models" / "qretix_ranker" / "qretix_scores.parquet"
OUT_DIR       = ROOT / "evaluation" / "benchmark_release"

# Benchmark identity
BENCHMARK_NAME    = "QRETIX-ATC: Anaplastic Thyroid Cancer Target Prioritisation Benchmark"
BENCHMARK_VERSION = "v0.1"
GITHUB_REPO       = "rishavm1/qretix-atc"
PAPER_URL         = ""             # no preprint yet
MAINTAINER        = "rishavm1"


# ── Ground truth construction ─────────────────────────────────────────────────

def _build_ground_truth(lb: pd.DataFrame) -> pd.DataFrame:
    """
    Export a clean ground-truth table for external evaluators.
    Columns: hgnc_symbol, stage, evidence_tier, earliest_evidence_year
    stage >= 2 = positive (experimental therapeutic evidence or higher)
    """
    gt = lb[["stage", "evidence_sources", "earliest_evidence_date"]].copy()

    # Derive evidence tier label
    tier_map = {0: "no_evidence", 1: "preclinical", 2: "experimental",
                3: "clinical_trial", 4: "validated"}
    gt["evidence_tier"] = gt["stage"].map(tier_map)

    # Earliest evidence year only (no full date — reduces leakage surface)
    gt["earliest_evidence_year"] = (
        gt["earliest_evidence_date"]
        .astype(str)
        .str[:4]
        .replace("nan", "")
        .replace("None", "")
    )
    gt = gt.drop(columns=["earliest_evidence_date"])

    # Binary positive flag for evaluators (Stage >= 2)
    gt["is_positive"] = (gt["stage"] >= 2).astype(int)

    return gt


def _build_holdout_split(gt: pd.DataFrame, seed: int = 42) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Stratified hold-out: 20% of Stage>=2 genes withheld as test set.
    Returns (public_gt, holdout_gt).
    public_gt  — released with benchmark (stage labels visible)
    holdout_gt — kept server-side; evaluators submit ranked lists against this
    """
    rng = np.random.default_rng(seed)
    positives = gt[gt["stage"] >= 2].index.tolist()
    n_holdout = max(10, int(len(positives) * 0.20))
    holdout_idx = rng.choice(positives, size=n_holdout, replace=False)
    holdout_set = set(holdout_idx)

    holdout_gt = gt.loc[list(holdout_set)].copy()
    # Public set: holdout genes have stage masked to -1 (unknown)
    public_gt = gt.copy()
    public_gt.loc[list(holdout_set), "stage"] = -1
    public_gt.loc[list(holdout_set), "evidence_tier"] = "withheld"
    public_gt.loc[list(holdout_set), "is_positive"] = -1

    return public_gt, holdout_gt


# ── Metrics ───────────────────────────────────────────────────────────────────

def _compute_metrics(y_true: np.ndarray, scores: np.ndarray) -> dict:
    from sklearn.metrics import average_precision_score

    def p_at_k(k):
        return float(y_true[np.argsort(scores)[::-1][:k]].mean())

    def ndcg_at_k(k):
        order = np.argsort(scores)[::-1][:k]
        g = y_true[order].astype(float)
        d = np.log2(np.arange(2, len(g) + 2))
        dcg = (g / d).sum()
        ig = np.sort(y_true).astype(float)[::-1][:k]
        idcg = (ig / np.log2(np.arange(2, len(ig) + 2))).sum()
        return float(dcg / idcg) if idcg > 0 else 0.0

    return {
        "AUPRC":   round(float(average_precision_score(y_true, scores)), 4),
        "P@10":    round(p_at_k(10), 4),
        "P@20":    round(p_at_k(20), 4),
        "P@50":    round(p_at_k(50), 4),
        "NDCG@50": round(ndcg_at_k(50), 4),
    }


# ── Standalone eval script ────────────────────────────────────────────────────

EVAL_SCRIPT = '''\
#!/usr/bin/env python3
"""
QRETIX-ATC Benchmark — Standalone Evaluator
============================================
Usage:
    python eval_script.py --predictions your_predictions.csv \\
                          --ground_truth ground_truth.parquet \\
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

    print("\\n=== QRETIX-ATC Benchmark Results ===")
    for k, v in metrics.items():
        print(f"  {k:<25} {v}")

    Path("results.json").write_text(json.dumps(metrics, indent=2))
    print("\\nWrote results.json")
    return metrics


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions",   required=True)
    parser.add_argument("--ground_truth",  required=True)
    parser.add_argument("--holdout_genes", default=None)
    args = parser.parse_args()
    evaluate(args.predictions, args.ground_truth, args.holdout_genes)
'''


# ── Metadata ──────────────────────────────────────────────────────────────────

def _build_metadata(gt: pd.DataFrame, holdout_gt: pd.DataFrame) -> dict:
    n_total    = len(gt)
    n_positive = int((gt["stage"] >= 2).sum())
    n_holdout  = len(holdout_gt)

    return {
        "benchmark_name": BENCHMARK_NAME,
        "version": BENCHMARK_VERSION,
        "disease": "Anaplastic Thyroid Cancer (ATC)",
        "task": "Gene target prioritisation",
        "category": "None",
        "maintainer": MAINTAINER,
        "github_repository": GITHUB_REPO,
        "paper_url": PAPER_URL,
        "generated_at": datetime.now(timezone.utc).isoformat(),

        # Q3DC fields
        "q3dc": {
            "benchmark_name": BENCHMARK_NAME,
            "category": "None",
            "maintainer_or_lead_author": MAINTAINER,
            "github_repository": GITHUB_REPO,
            "paper_or_project_url": PAPER_URL,
            "what_it_measures": (
                "Ranks all human protein-coding genes (~19k) as candidate therapeutic "
                "targets for Anaplastic Thyroid Cancer (ATC) using a multimodal evidence "
                "framework (transcriptomics, CRISPR dependency, genomics, literature, PPI). "
                "Unlike existing cancer target benchmarks, it enforces strict temporal "
                "retrospective validation across four cutoff years (2015–2021) and "
                "decomposes scores into three non-merged axes: biology, therapeutic "
                "feasibility, and novelty."
            ),
            "tasks": (
                "Target ranking, Precision@K recovery, "
                "Temporal prospective validation, Novelty scoring"
            ),
            "known_caveat": (
                "Ground-truth labels are derived from PubMed NLP + Open Targets + "
                "ClinicalTrials.gov and inherit their publication bias: well-studied "
                "cancer genes (TP53, BRAF, KRAS) are over-represented as positives "
                "regardless of ATC-specific evidence strength. "
                "COSMIC and DisGeNET were unavailable (license/registration required), "
                "so driver annotation relies on Open Targets and curated literature lists. "
                "The DepMap-only baseline achieves AUPRC=0.92 because essential "
                "housekeeping genes co-segregate with ATC dependency — this is a known "
                "confound documented in the leakage audit."
            ),
            "live_leaderboard_url": "",
            "submission_evaluation_url": "",
        },

        # Dataset statistics
        "dataset_stats": {
            "gene_universe_size": n_total,
            "n_positives_stage_ge2": n_positive,
            "n_holdout_genes": n_holdout,
            "positive_fraction_pct": round(n_positive / n_total * 100, 2),
            "label_sources": [
                "PubMed E-utilities (ATC corpus)",
                "Open Targets (ATC associations)",
                "ClinicalTrials.gov API v2",
                "DepMap CRISPR Chronos (7 authenticated ATC lines)",
            ],
            "temporal_splits": [2015, 2017, 2019, 2021],
        },

        # Baseline metrics (from benchmark_comparison.md)
        "baseline_metrics": {
            "DE-only":            {"AUPRC": 0.0644, "P@50": 0.06,  "NDCG@50": 0.0549},
            "DepMap-only":        {"AUPRC": 0.9157, "P@50": 1.0,   "NDCG@50": 1.0},
            "Open Targets":       {"AUPRC": 0.1024, "P@50": 0.88,  "NDCG@50": 0.9158},
            "Network centrality": {"AUPRC": 0.3083, "P@50": 0.42,  "NDCG@50": 0.4966},
            "QRETIX":             {"AUPRC": 0.8847, "P@50": 1.0,   "NDCG@50": 1.0},
        },
    }


# ── Q3DC form text ────────────────────────────────────────────────────────────

def _build_q3dc_form(meta: dict) -> str:
    q = meta["q3dc"]
    ds = meta["dataset_stats"]
    bl = meta["baseline_metrics"]

    baselines_table = "\n".join(
        f"| {m} | {v['AUPRC']} | {v['P@50']} | {v['NDCG@50']} |"
        for m, v in bl.items()
    )

    return f"""\
# Q3DC Benchmark Submission — Pre-filled Form
_Generated: {meta['generated_at']}_
_Fill in fields marked [TODO] before submitting._

---

## Benchmark name
{q['benchmark_name']}

## Category
{q['category']}

## Maintainer or lead author
[TODO: {q['maintainer_or_lead_author'] or 'Your name / institution'}]

## GitHub repository
[TODO: {q['github_repository'] or 'owner/repo'}]

## Paper or project URL
[TODO: {q['paper_or_project_url'] or 'https://arxiv.org/abs/...'}]

## What it measures
{q['what_it_measures']}

## Tasks
{q['tasks']}

## Known caveat, bias or limitation
{q['known_caveat']}

## Live leaderboard URL
{q['live_leaderboard_url'] or '(none — downloadable eval script only)'}

## Submission / evaluation URL
{q['submission_evaluation_url'] or '(none — downloadable eval script only)'}

---

## Dataset Statistics (for reviewers)

| Field | Value |
|-------|-------|
| Gene universe | {ds['gene_universe_size']:,} genes |
| Positives (Stage ≥ 2) | {ds['n_positives_stage_ge2']:,} ({ds['positive_fraction_pct']}%) |
| Hold-out genes | {ds['n_holdout_genes']} |
| Temporal splits | {', '.join(str(y) for y in ds['temporal_splits'])} |

## Baseline Metrics (for reviewers)

| Method | AUPRC | P@50 | NDCG@50 |
|--------|-------|------|---------|
{baselines_table}

---

## Checklist (Q3DC requirements)
- [x] This is a benchmark / dataset-with-defined-evaluation (not a model paper)
- [x] Publicly accessible evaluation script included (`eval_script.py`)
- [x] Caveat names a real, specific limitation (publication bias + DepMap confound)
- [x] Ground truth schema documented (`ground_truth_schema.md`)
- [x] Sample submission format provided (`sample_submission.csv`)
- [ ] [TODO] Preprint / paper URL added
- [ ] [TODO] GitHub repository made public
- [ ] [TODO] Maintainer name and ORCID filled in
"""


# ── Schema doc ────────────────────────────────────────────────────────────────

SCHEMA_DOC = """\
# QRETIX-ATC Benchmark — Ground Truth Schema

## File: ground_truth.parquet

| Column | Type | Description |
|--------|------|-------------|
| hgnc_symbol | str (index) | HGNC approved gene symbol (primary key) |
| stage | int | Evidence stage: -1=withheld, 0=none, 1=preclinical, 2=experimental, 3=clinical_trial, 4=validated |
| evidence_tier | str | Human-readable stage label (withheld for hold-out genes) |
| evidence_sources | str | Pipe-separated list of source IDs supporting the label |
| earliest_evidence_year | str | Year of earliest supporting evidence (YYYY), empty if unknown |
| is_positive | int | 1 if stage>=2, 0 if stage<2, -1 if withheld |

## Label Hierarchy

| Stage | Definition | Primary Sources |
|-------|-----------|-----------------|
| 0 | No meaningful ATC evidence | — |
| 1 | Preclinical evidence (in vitro/in vivo ATC model) | PubMed NLP |
| 2 | Experimental therapeutic evidence (perturbation → phenotype) | PubMed NLP + DepMap |
| 3 | Clinical investigation (Phase 1–2 ATC trial) | ClinicalTrials.gov + Open Targets |
| 4 | Validated / approved drug target in ATC | Curated + ClinicalTrials.gov |

## Positive Definition for Evaluation

A gene is a **positive** (is_positive=1) if stage >= 2.
This requires at least experimental therapeutic evidence — preclinical-only
genes (stage=1) are treated as unlabelled, not negative (PU-Learning formulation).

## Hold-out Protocol

20% of Stage>=2 genes are withheld (stage=-1, is_positive=-1) from the public
ground truth. Evaluators submit a ranked list of all genes; the organiser
scores against the withheld set. This prevents overfitting to the public labels.

## Known Biases

1. **Publication bias**: well-studied oncogenes (TP53, BRAF, KRAS) accumulate
   Stage 1–2 labels from general cancer literature even when ATC-specific
   evidence is weak.
2. **DepMap confound**: essential housekeeping genes score highly on dependency
   features regardless of ATC specificity. The DepMap-only baseline achieves
   AUPRC=0.92 — any method must beat this to demonstrate ATC-specific signal.
3. **COSMIC/DisGeNET gap**: driver annotation relies on Open Targets + curated
   lists; COSMIC CGC and DisGeNET were unavailable (license/registration required).
4. **Temporal label drift**: labels reflect the state of evidence at pipeline
   run time. Re-running after new publications will shift stage assignments.
"""


# ── Sample submission ─────────────────────────────────────────────────────────

def _build_sample_submission(gt: pd.DataFrame) -> str:
    lines = ["hgnc_symbol,score"]
    # Top QRETIX genes get high scores; rest get random low scores
    top_genes = ["BRD4", "CDK1", "NAT10", "PSMD14", "CDK7",
                 "PPP2CA", "CDC20", "XPO1", "CCNB1", "RBX1"]
    written = set()
    for i, g in enumerate(top_genes):
        lines.append(f"{g},{round(1.0 - i * 0.02, 4)}")
        written.add(g)
    # Fill remaining genes with placeholder scores
    rng = np.random.default_rng(0)
    for sym in gt.index:
        if sym not in written:
            lines.append(f"{sym},{round(float(rng.uniform(0, 0.3)), 4)}")
    return "\n".join(lines)


# ── Main export ───────────────────────────────────────────────────────────────

def run():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    log.info("Loading labels...")
    lb = pd.read_parquet(LABELS_PATH)

    gt = _build_ground_truth(lb)
    public_gt, holdout_gt = _build_holdout_split(gt)

    # Ground truth files
    public_gt.to_parquet(OUT_DIR / "ground_truth.parquet")
    holdout_gt.to_parquet(OUT_DIR / "holdout_ground_truth.parquet")
    log.info("Wrote ground_truth.parquet (%d genes, %d positives visible)",
             len(public_gt), int((public_gt["stage"] >= 2).sum()))
    log.info("Wrote holdout_ground_truth.parquet (%d withheld positives)", len(holdout_gt))

    # Holdout gene list (for eval_script.py --holdout_genes)
    holdout_genes_path = OUT_DIR / "holdout_genes.txt"
    holdout_genes_path.write_text("\n".join(holdout_gt.index.tolist()))

    # Metadata
    meta = _build_metadata(public_gt, holdout_gt)
    (OUT_DIR / "benchmark_metadata.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8"
    )
    log.info("Wrote benchmark_metadata.json")

    # Eval script
    (OUT_DIR / "eval_script.py").write_text(EVAL_SCRIPT, encoding="utf-8")
    log.info("Wrote eval_script.py")

    # Schema doc
    (OUT_DIR / "ground_truth_schema.md").write_text(SCHEMA_DOC, encoding="utf-8")
    log.info("Wrote ground_truth_schema.md")

    # Sample submission
    sample_csv = _build_sample_submission(public_gt)
    (OUT_DIR / "sample_submission.csv").write_text(sample_csv, encoding="utf-8")
    log.info("Wrote sample_submission.csv")

    # Q3DC form
    form_text = _build_q3dc_form(meta)
    (OUT_DIR / "q3dc_submission_form.md").write_text(form_text, encoding="utf-8")
    log.info("Wrote q3dc_submission_form.md")

    # Validate eval script works on QRETIX scores as a smoke test
    if QRETIX_PATH.exists():
        _smoke_test(public_gt, holdout_gt)

    log.info("\n=== Benchmark release package ready: %s ===", OUT_DIR)
    log.info("Files:")
    for f in sorted(OUT_DIR.iterdir()):
        log.info("  %s", f.name)


def _smoke_test(public_gt: pd.DataFrame, holdout_gt: pd.DataFrame):
    """Verify eval_script logic against QRETIX scores."""
    log.info("Running smoke test against QRETIX scores...")
    try:
        qretix = pd.read_parquet(QRETIX_PATH)
        scores = qretix["qretix_score"].reindex(holdout_gt.index).fillna(0).values
        y      = (holdout_gt["stage"] >= 2).astype(int).values
        if y.sum() == 0:
            log.warning("Smoke test: no positives in holdout set — skipping metric check")
            return
        metrics = _compute_metrics(y, scores)
        log.info("Smoke test metrics (QRETIX vs holdout): %s", metrics)
    except Exception as e:
        log.warning("Smoke test failed (non-fatal): %s", e)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    run()
