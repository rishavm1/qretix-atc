# Model Comparison Report — Phase 5

## All Models — Metrics Summary

| Model | AUPRC | P@10 | P@20 | P@50 | NDCG@50 |
|-------|-------|------|------|------|---------|
| LogisticRegression_PU | 0.9407 | 1.0 | 1.0 | 1.0 | 1.0 |
| RandomForest_PU_Bagging | 0.9938 | 1.0 | 1.0 | 1.0 | 1.0 |
| XGBoost_rank:pairwise | 0.8450 | 1.0 | 1.0 | 1.0 | 1.0 |
| LightGBM_lambdarank | 0.9605 | 1.0 | 1.0 | 1.0 | 1.0 |
| GraphSAGE_PU | deferred | — | — | — | — |
| MultimodalFusion | 0.9631 | 1.0 | 1.0 | 1.0 | 1.0 |

## Notes
- PU formulation: positives = Stage>=1, unlabelled = Stage=0 (never negative)
- Best single model: RandomForest_PU_Bagging (AUPRC=0.9938)
- GNN deferred: dense adjacency for 19k genes exceeds CPU budget (logged in decisions_log.md)
- Fusion weights optimized on full dataset (SLSQP, AUPRC objective)

## Scoring Architecture (per spec)
Three scores — NEVER merged into one opaque number:
- biology_score: expression + genetics + dependency
- therapeutic_score: druggability + modulators + structure
- novelty_score: computed in Phase 6 (5 independent sub-scores)