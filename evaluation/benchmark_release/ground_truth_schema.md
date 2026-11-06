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
