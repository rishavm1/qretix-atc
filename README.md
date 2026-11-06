# QRETIX-ATC 

**Anaplastic Thyroid Cancer (ATC) Therapeutic Target Prioritisation Benchmark**

An auditable, evidence-provenance-aware pipeline that ranks human protein-coding genes as candidate ATC therapeutic targets, with strict separation of biology, therapeutic feasibility, and novelty scores.

---

## What This Is

QRETIX-ATC is both a **pipeline** and a **benchmark**. It:

- Ranks ~19k–41k human protein-coding genes as ATC drug targets
- Uses multimodal evidence: scRNA-seq, bulk RNA, CRISPR dependency (DepMap), genomics, PPI networks, literature, clinical trials
- Enforces temporal retrospective validation across 4 cutoff years (2015, 2017, 2019, 2021)
- Decomposes every score into 3 non-merged axes: `biology_score`, `therapeutic_score`, `novelty_score`
- Ships a standalone benchmark evaluator for external methods to compare against

---

## Benchmark (Q3DC Submission)

The benchmark package is in `evaluation/benchmark_release/`:

| File | Description |
|------|-------------|
| `ground_truth.parquet` | 41,217 genes with stage labels (Stage 0–4) |
| `eval_script.py` | Standalone evaluator — no pipeline dependency |
| `ground_truth_schema.md` | Column definitions, label hierarchy, known biases |
| `sample_submission.csv` | Example prediction format |
| `benchmark_metadata.json` | Machine-readable benchmark metadata |

### Evaluate your method

```bash
pip install pandas scikit-learn pyarrow

python evaluation/benchmark_release/eval_script.py \
    --predictions your_predictions.csv \
    --ground_truth evaluation/benchmark_release/ground_truth.parquet
```

**Prediction format** (`your_predictions.csv`):
```
hgnc_symbol,score
BRD4,0.98
CDK1,0.95
...
```

One row per gene, higher score = more likely ATC target. Include all ~19k–41k protein-coding genes.

### Baseline metrics

| Method | AUPRC | P@50 | NDCG@50 |
|--------|-------|------|---------|
| DE-only | 0.0644 | 0.06 | 0.0549 |
| DepMap-only | 0.9157 | 1.0 | 1.0 |
| Open Targets | 0.1024 | 0.88 | 0.9158 |
| Network centrality | 0.3083 | 0.42 | 0.4966 |
| **QRETIX** | **0.8847** | **1.0** | **1.0** |

> **Known caveat:** Ground-truth labels inherit publication bias from PubMed/Open Targets — well-studied oncogenes (TP53, BRAF, KRAS) are over-represented. The DepMap-only baseline achieves AUPRC=0.92 due to housekeeping gene co-segregation with ATC dependency, a known confound documented in the leakage audit.

---

## Label Hierarchy

| Stage | Definition |
|-------|-----------|
| 0 | No meaningful ATC evidence |
| 1 | Preclinical evidence (in vitro/in vivo ATC model) |
| 2 | Experimental therapeutic evidence (perturbation → phenotype) |
| 3 | Clinical investigation (Phase 1–2 ATC trial) |
| 4 | Validated / approved drug target in ATC |

Positive = Stage ≥ 2. Stage 1 genes are treated as **unlabelled**, not negative (PU-Learning formulation).

---

## Pipeline Phases

| Phase | Name | Gate Artifact |
|-------|------|---------------|
| 0 | Project Scaffolding | Directory tree + core modules |
| 1 | Data Acquisition | `data/provenance/*.json` + `docs/dataset_inventory.md` |
| 2 | QC & Preprocessing | `docs/qc_report.md` |
| 3 | Feature Engineering | `features/feature_table.parquet` |
| 4 | Label Construction | `labels/labels.parquet` |
| 5 | Modular Scoring | `evaluation/reports/model_comparison.md` |
| 6 | Novelty Engine | `features/novelty_scores.parquet` |
| 7 | Temporal Retrospective Validation | `evaluation/temporal_splits/{year}/results.md` × 4 |
| 8 | Evidence Cards | `evidence_cards/{GENE}.md` × 51 |
| 9 | Benchmarking + Export | `evaluation/benchmark_release/` |

---

## Data Sources

- **Transcriptomics:** GEO GSE193581 (scRNA-seq), GSE76039 (bulk microarray)
- **Genomics:** cBioPortal, DepMap (7 authenticated ATC lines — contaminated lines excluded per Schweppe et al. 2008)
- **Functional dependency:** DepMap CRISPR Chronos scores
- **Knowledge graph:** STRING, UniProt, Reactome, Open Targets, ChEMBL
- **Literature:** PubMed E-utilities, ClinicalTrials.gov API v2
- **Blocked/gated:** COSMIC (license required), DrugBank (license required), DisGeNET (registration required)

---

## Installation

```bash
pip install -r requirements.txt
```

Run the full pipeline:
```bash
python -m src.pipeline run --all
```

Run a single phase:
```bash
python -m src.pipeline run --phase 9
```

Check status:
```bash
python -m src.pipeline status
```

---

## Standing Rules

1. Never invent or simulate a dataset — report gaps instead
2. Never let a post-cutoff data point enter a pre-cutoff retrospective run
3. Never merge biology, novelty, and therapeutic feasibility into one opaque number
4. Every number must be traceable to a named file and source ID
5. Every deviation is logged in `docs/decisions_log.md`

---

## Top Ranked Candidates (QRETIX v0.1)

BRD4, CDK1, NAT10, PSMD14, CDK7, PPP2CA, CDC20, XPO1, CCNB1, RBX1, MYC, CCNA2, POLE, PCNA, KIF11, AURKB, MTOR, BUB1B, AURKA, KIF23

Evidence cards for each candidate are in `evidence_cards/`.

---

## Maintainer

GitHub: [@rishavm1](https://github.com/rishavm1)
