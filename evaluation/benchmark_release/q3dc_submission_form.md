# Q3DC Benchmark Submission — Pre-filled Form
_Generated: 2026-11-06T12:34:54.790348+00:00_
_Fill in fields marked [TODO] before submitting._

---

## Benchmark name
QRETIX-ATC: Anaplastic Thyroid Cancer Target Prioritisation Benchmark

## Category
None

## Maintainer or lead author
[TODO: rishavm1]

## GitHub repository
[TODO: rishavm1/qretix-atc]

## Paper or project URL
[TODO: https://arxiv.org/abs/...]

## What it measures
Ranks all human protein-coding genes (~19k) as candidate therapeutic targets for Anaplastic Thyroid Cancer (ATC) using a multimodal evidence framework (transcriptomics, CRISPR dependency, genomics, literature, PPI). Unlike existing cancer target benchmarks, it enforces strict temporal retrospective validation across four cutoff years (2015–2021) and decomposes scores into three non-merged axes: biology, therapeutic feasibility, and novelty.

## Tasks
Target ranking, Precision@K recovery, Temporal prospective validation, Novelty scoring

## Known caveat, bias or limitation
Ground-truth labels are derived from PubMed NLP + Open Targets + ClinicalTrials.gov and inherit their publication bias: well-studied cancer genes (TP53, BRAF, KRAS) are over-represented as positives regardless of ATC-specific evidence strength. COSMIC and DisGeNET were unavailable (license/registration required), so driver annotation relies on Open Targets and curated literature lists. The DepMap-only baseline achieves AUPRC=0.92 because essential housekeeping genes co-segregate with ATC dependency — this is a known confound documented in the leakage audit.

## Live leaderboard URL
(none — downloadable eval script only)

## Submission / evaluation URL
(none — downloadable eval script only)

---

## Dataset Statistics (for reviewers)

| Field | Value |
|-------|-------|
| Gene universe | 41,217 genes |
| Positives (Stage ≥ 2) | 1,075 (2.61%) |
| Hold-out genes | 268 |
| Temporal splits | 2015, 2017, 2019, 2021 |

## Baseline Metrics (for reviewers)

| Method | AUPRC | P@50 | NDCG@50 |
|--------|-------|------|---------|
| DE-only | 0.0644 | 0.06 | 0.0549 |
| DepMap-only | 0.9157 | 1.0 | 1.0 |
| Open Targets | 0.1024 | 0.88 | 0.9158 |
| Network centrality | 0.3083 | 0.42 | 0.4966 |
| QRETIX | 0.8847 | 1.0 | 1.0 |

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
