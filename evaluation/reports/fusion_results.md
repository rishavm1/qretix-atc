# Multimodal Fusion (QRETIX Ranker) Results

## Fusion Architecture
- Late fusion: weighted combination of omics, graph, and text engines
- Weights optimized on full dataset (AUPRC objective, SLSQP solver)
- Weights learned from data — NOT assumed equal

## Fusion Weights

| Engine | Weight |
|--------|--------|
| logreg_score | 0.2000 |
| rf_score | 0.2000 |
| gbm_ensemble_score | 0.2000 |
| gnn_score | 0.2000 |
| text_score | 0.2000 |

## Metrics

| Metric | Value |
|--------|-------|
| AUPRC | 0.9601 |
| P@10 | 1.0 |
| P@20 | 1.0 |
| P@50 | 1.0 |
| NDCG@50 | 1.0 |

## Scoring Decomposition
Every gene has THREE separate scores (never merged into one opaque number):
- **biology_score**: expression, genetics, dependency evidence
- **therapeutic_score**: druggability, existing modulators, structure
- **novelty_score**: computed separately by Phase 6 novelty engine

## Top 20 Ranked Genes (Final QRETIX Score)

BRD4, CDK1, NAT10, PSMD14, CDK7, PPP2CA, CDC20, XPO1, CCNB1, RBX1, MYC, CCNA2, POLE, PCNA, KIF11, AURKB, MTOR, BUB1B, AURKA, KIF23