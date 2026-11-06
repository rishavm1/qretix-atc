# QRETIX-ATC v0.1 — Living Specification

**Project:** QRETIX-ATC-v0.1  
**Disease:** Anaplastic Thyroid Cancer (ATC)  
**Version:** 0.1  
**Status:** Active development  

---

## Objective

Build an auditable, evidence-provenance-aware pipeline that ranks human
protein-coding genes as candidate ATC therapeutic targets, explicitly separating
KNOWN biology from NOVEL hypotheses.

**Non-goal:** Do not produce a single opaque "target score." Every score must
decompose into traceable evidence components.

---

## Standing Rules (enforced throughout all phases)

1. Never invent or simulate a dataset — report gaps instead.
2. Never let a post-cutoff data point enter a pre-cutoff retrospective run.
3. Never merge biological confidence, novelty, and therapeutic feasibility into one opaque number.
4. Never treat "low publication count" alone as evidence of novelty.
5. Every number in every output must be traceable to a named file and source ID.
6. Every deviation, judgment call, or blocked data source must be logged in `docs/decisions_log.md`.

---

## Phase Summary

| Phase | Name | Gate Artifact |
|-------|------|---------------|
| 0 | Project Scaffolding | Directory tree + core modules |
| 1 | Data Inventory & Acquisition | `data/provenance/*.json` + `docs/dataset_inventory.md` |
| 2 | QC & Preprocessing | `docs/qc_report.md` + processed files |
| 3 | Feature Engineering | `features/feature_table.parquet` + `docs/feature_spec.md` |
| 4 | Label Construction | `labels/labels.parquet` + `labels/label_spec.md` |
| 5 | Modular Scoring Architecture | `evaluation/reports/model_comparison.md` |
| 6 | Novelty Engine | `features/novelty_scores.parquet` (5 sub-scores) |
| 7 | Temporal Retrospective Validation | `evaluation/temporal_splits/{year}/results.md` × 4 |
| 8 | Evidence Cards | `evidence_cards/{GENE}.md` per candidate |
| 9 | Benchmarking | `evaluation/reports/benchmark_comparison.md` |

---

## Label Hierarchy

| Stage | Definition |
|-------|------------|
| 0 | No meaningful ATC evidence |
| 1 | Preclinical evidence (in vitro/in vivo, ATC or ATC-relevant model) |
| 2 | Experimental therapeutic evidence (perturbation → phenotype in ATC) |
| 3 | Clinical investigation (trial exists, any phase) |
| 4 | Clinical-stage / validated target |

---

## Scoring Architecture

Three separate scores — never merged into one number:

- **Biology Score** — disease relevance + tumor-cell specificity + genetic evidence + functional dependency
- **Therapeutic Feasibility Score** — druggability + existing modulators + structural evidence
- **Novelty Score** — five independent sub-scores (see Phase 6)

---

## Novelty Sub-scores (Phase 6)

All five must be populated independently per candidate gene:

1. Disease novelty
2. Target novelty
3. Mechanistic novelty
4. Combination novelty
5. Research saturation

---

## Temporal Retrospective Validation Cutoffs (Phase 7)

2015 · 2017 · 2019 · 2021

Each cutoff gets a fully isolated pre-cutoff dataset rebuild. Post-cutoff PMIDs
must be literally absent from input files, verified by `leakage_auditor.py`.

---

## Data Sources

### Transcriptomics
- GSE193581 (scRNA-seq, 10 ATC / 7 PTC / 6 normal)
- GSE76039 (bulk microarray, PDTC/ATC)
- GEO search: anaplastic thyroid cancer

### Genomics
- cBioPortal: thca_gatci_cellrep_2024 (190 ATC), thca_mskcc_2016
- GDC TCGA-THCA (PTC reference only — zero ATC cases)
- COSMIC: **manual approval required** — flagged as gap

### Epigenomics / Proteomics
- ENCODE (normal thyroid reference only — no ATC biosamples expected)
- PRIDE (sparse ATC proteomics — log sparsity explicitly)

### Functional Dependency
- DepMap: CRISPR Chronos scores, authenticated ATC cell lines only
- Excluded contaminated lines: ARO, DRO, KAT-4, NPA, FRO (Schweppe et al. 2008)

### Knowledge Graph
- STRING, UniProt, Reactome, Open Targets, ChEMBL, Gene Ontology
- BioGRID: requires free API key — decision pending (see decisions_log.md)
- DisGeNET: registration required — flagged as gap
- DrugBank: licensed — flagged as gap

### Literature & Clinical Trials
- PubMed E-utilities (ATC-specific + broad thyroid corpora)
- ClinicalTrials.gov API v2

---

## Key Architectural Decisions

1. **Central Gene Universe** — HGNC protein-coding set → `data/processed/gene_universe.parquet`. All joins go through this table.
2. **PU-Learning formulation** — unlabeled genes are NOT treated as negatives.
3. **GATCI 2024 cohort** — largest unified ATC genomic dataset (190 tumors), added as high-priority source.
4. **TCGA-THCA explicit downgrade** — used as differentiated thyroid baseline only, never as ATC evidence.
5. **Druggability gating** — post-ranking tier assignment (A/B/C), does not affect ranking score.
6. **Leakage Auditor** — automated checks run before every model training step.

---

*This document is updated as decisions are made. See `docs/decisions_log.md` for all deviations.*
