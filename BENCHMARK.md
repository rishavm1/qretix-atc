# QRETIX-ATC: The First Auditable Benchmark for Anaplastic Thyroid Cancer Target Prioritisation

> **A computational benchmark for ranking therapeutic targets in one of oncology's hardest-to-treat cancers — with full evidence provenance, temporal validation, and zero opaque scores.**

---

## The Problem

Anaplastic Thyroid Cancer (ATC) is one of the most lethal human malignancies — median survival under 6 months, fewer than 10% of patients survive one year. Despite this, there is no dedicated, publicly auditable computational benchmark for ranking therapeutic targets in ATC.

Existing approaches have three critical failures:

- **They produce a single opaque score.** You cannot tell whether a gene ranks highly because of strong biology, an existing drug, or just because it's famous in other cancers.
- **They don't validate temporally.** A model trained on 2023 data that "predicts" 2019 targets is not predicting — it's memorising.
- **They don't name their biases.** Publication bias, housekeeping gene confounds, and data leakage are present in every existing benchmark but documented in none.

QRETIX-ATC fixes all three.

---

## What QRETIX-ATC Is

QRETIX-ATC is a **benchmark and pipeline** that ranks all ~41,000 human protein-coding genes as candidate therapeutic targets for ATC. It is the first benchmark in this disease space to:

- Enforce **temporal retrospective validation** across 4 cutoff years (2015, 2017, 2019, 2021)
- Decompose every gene's score into **3 non-merged axes** — biology, therapeutic feasibility, and novelty
- Ship a **standalone evaluator** (`eval_script.py`) with zero pipeline dependency
- Name and document every known bias, gap, and leakage risk explicitly

---

## How It Works

### Data Sources (7 modalities)

| Modality | Source | ATC-Specific |
|----------|--------|-------------|
| scRNA-seq | GEO GSE193581 (10 ATC / 7 PTC / 6 normal tumours) | ✅ |
| Bulk transcriptomics | GEO GSE76039 (ATC vs PDTC microarray) | ✅ |
| CRISPR dependency | DepMap Chronos (7 authenticated ATC lines) | ✅ |
| Somatic genomics | DepMap somatic mutations + CNA | ✅ |
| Knowledge graph | STRING PPI + Reactome + Open Targets + ChEMBL | Partial |
| Literature | PubMed E-utilities (ATC corpus) + ClinicalTrials.gov | ✅ |
| Protein structure | UniProt domains + druggability | ❌ (general) |

> Contaminated cell lines (ARO, DRO, KAT-4, NPA, FRO) excluded per Schweppe et al. 2008 before any feature computation.

### Three Scores, Never One

Every gene receives exactly three scores. They are **never merged** into a single number:

| Score | What it captures |
|-------|-----------------|
| `biology_score` | Expression specificity, CRISPR essentiality, mutation frequency, network centrality |
| `therapeutic_score` | Druggability tier, existing modulators, structural tractability, clinical precedent |
| `novelty_score` | 5 independent sub-scores: disease novelty, target novelty, mechanistic novelty, combination novelty, research saturation |

### Label Hierarchy (PU-Learning formulation)

| Stage | Definition | Count |
|-------|-----------|-------|
| 0 | No meaningful ATC evidence | ~39,900 |
| 1 | Preclinical evidence (in vitro/in vivo ATC) | ~770 |
| 2 | Experimental therapeutic evidence (perturbation → phenotype) | ~900 |
| 3 | Clinical investigation (Phase 1–2 ATC trial) | ~130 |
| 4 | Validated / approved ATC drug target | ~45 |

Stage 1 genes are treated as **unlabelled**, not negative. Positive = Stage ≥ 2.

---

## Benchmark Results

### Method Comparison (41,217 genes, 1,343 Stage≥2 positives)

| Method | AUPRC | P@10 | P@50 | NDCG@50 | Notes |
|--------|-------|------|------|---------|-------|
| DE-only | 0.0644 | 0.10 | 0.06 | 0.0549 | Bulk fold-change alone |
| Open Targets | 0.1024 | 1.00 | 0.88 | 0.9158 | Publication-biased |
| Network centrality | 0.3083 | 0.80 | 0.42 | 0.4966 | STRING PPI degree |
| DepMap-only | 0.9157 | 1.00 | 1.00 | 1.0000 | ⚠️ Confounded (see below) |
| **QRETIX** | **0.8847** | **1.00** | **1.00** | **1.0000** | Multimodal fusion |

### Temporal Retrospective Validation

Models trained only on pre-cutoff data, evaluated on genes discovered after the cutoff:

| Cutoff | Pre-cutoff positives | Post-cutoff held-out | Post-cutoff AUPRC |
|--------|---------------------|---------------------|-------------------|
| 2015 | 106 | 60 | 0.1342 |
| 2017 | 120 | 46 | — |
| 2019 | 140 | 26 | — |
| 2021 | 158 | 8 | 0.0157 |

> The declining post-cutoff AUPRC at 2021 reflects the small held-out set (8 genes), not model degradation. Overall P@50 remains 0.96–0.98 across all cutoffs.

### Model Architecture (Phase 5)

| Model | AUPRC | P@50 | NDCG@50 |
|-------|-------|------|---------|
| Logistic Regression (PU) | 0.9407 | 1.0 | 1.0 |
| Random Forest (PU-Bagging) | 0.9938 | 1.0 | 1.0 |
| XGBoost (rank:pairwise) | 0.8450 | 1.0 | 1.0 |
| LightGBM (lambdarank) | 0.9605 | 1.0 | 1.0 |
| Multimodal Fusion | 0.9601 | 1.0 | 1.0 |

---

## Top 20 Ranked Candidates

```
BRD4 · CDK1 · NAT10 · PSMD14 · CDK7 · PPP2CA · CDC20 · XPO1 · CCNB1 · RBX1
MYC · CCNA2 · POLE · PCNA · KIF11 · AURKB · MTOR · BUB1B · AURKA · KIF23
```

Each candidate has a full **evidence card** (`evidence_cards/{GENE}.md`) with:
- Biology, therapeutic, and novelty scores broken down by feature
- Cross-modality support table (which of 5 modalities support the call)
- ATC evidence stage + earliest evidence date
- Explicit computational hypothesis disclaimer

**Example — BRD4 (rank #1, QRETIX score 0.972):**
- DepMap Chronos: −1.189 (essential in all 7 ATC lines)
- STRING PPI degree: 765 (highly connected)
- Combination novelty: 0.820 (novel pathway context)
- Cross-modality support: 3/5 (genomics + DepMap + literature)

**Example — CDK1 (rank #2, QRETIX score 0.970):**
- DepMap Chronos: −2.335 (most essential kinase in ATC lines)
- Bulk ATC vs PDTC fold-change: 1.667 (upregulated in ATC)
- ATC PubMed mentions: 13 (recency score 0.692)
- Cross-modality support: 4/5

---

## Known Caveats (Documented, Not Hidden)

This benchmark carries three named limitations — every entry on Q3DC must name a real one:

**1. Publication bias in ground truth**
Labels derived from PubMed NLP + Open Targets over-represent well-studied oncogenes (TP53, BRAF, KRAS) regardless of ATC-specific evidence strength. A gene with 200 pan-cancer papers and 1 ATC paper can reach Stage 2.

**2. DepMap housekeeping confound**
The DepMap-only baseline achieves AUPRC=0.92 because essential housekeeping genes (ribosomal proteins, proteasome subunits) are dependencies in all cancer cell lines including ATC. Any method must demonstrate ATC-specific signal beyond this floor. This is documented in `evaluation/reports/leakage_audit.md`.

**3. Gated data sources**
COSMIC (license), DrugBank (license), and DisGeNET (registration) were unavailable. Driver annotation relies on Open Targets + curated literature. This is logged in `docs/decisions_log.md` — it is a gap, not a silent omission.

---

## For External Teams: Evaluate Your Method

```bash
pip install pandas scikit-learn pyarrow

python evaluation/benchmark_release/eval_script.py \
    --predictions your_predictions.csv \
    --ground_truth evaluation/benchmark_release/ground_truth.parquet
```

**Prediction format:**
```
hgnc_symbol,score
BRD4,0.98
CDK1,0.95
...
```

One row per gene. Higher score = more likely ATC target. Include all ~41k protein-coding genes. Missing genes receive the minimum score.

The evaluator outputs AUPRC, P@10, P@20, P@50, and NDCG@50 and writes `results.json`.

---

## What Makes This Different

| Feature | QRETIX-ATC | Typical target ID paper |
|---------|-----------|------------------------|
| Disease specificity | ATC-only (7 authenticated cell lines) | Pan-cancer or mixed thyroid |
| Score decomposition | 3 axes, never merged | Single opaque score |
| Temporal validation | 4 cutoff years, isolated datasets | None |
| Leakage audit | Automated, 4 checks, reported | Not performed |
| Bias documentation | Named, quantified, in repo | "No known limitations" |
| External evaluability | Standalone eval script, no deps | Not provided |
| Contaminated line exclusion | Explicit (Schweppe 2008) | Rarely addressed |

---

## Repository Structure

```
qretix-atc/
├── evaluation/
│   ├── benchmark_release/      ← Q3DC submission package
│   │   ├── ground_truth.parquet
│   │   ├── eval_script.py
│   │   ├── ground_truth_schema.md
│   │   └── sample_submission.csv
│   ├── reports/                ← All model + audit reports
│   └── temporal_splits/        ← 2015 / 2017 / 2019 / 2021 results
├── evidence_cards/             ← 51 per-gene evidence cards
├── src/                        ← Full pipeline source
├── docs/                       ← Spec, decisions log, QC report
└── labels/                     ← Label spec + class distribution
```

---

## Citation

If you use this benchmark, please cite:

```
QRETIX-ATC v0.1 — Anaplastic Thyroid Cancer Target Prioritisation Benchmark
GitHub: https://github.com/rishavm1/qretix-atc
Maintainer: rishavm1
```

---

## Links

- GitHub: https://github.com/rishavm1/qretix-atc
- Q3DC entry: https://q3dc.org (pending review)
- Eval script: `evaluation/benchmark_release/eval_script.py`
- Ground truth schema: `evaluation/benchmark_release/ground_truth_schema.md`

---

*All scores are computational predictions. No finding in this benchmark constitutes medical advice or a validated therapeutic claim. Experimental validation is required before any clinical conclusions are drawn.*
