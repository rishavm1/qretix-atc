# Proper Hold-out Evaluation Report

## Evaluation Protocol (No Circularity)
- Training set: 41,078 genes
- Hold-out set: 139 genes (never seen during training)
- Hold-out composition: 30% of Stage 4, 20% of Stage 3, 10% of Stage 2

## Hold-out Metrics (RF PU-Bagging)

| Metric | Value |
|--------|-------|
| AUPRC (hold-out) | 1.0000 |
| P@10 (hold-out)  | 1.0000 |
| P@50 (hold-out)  | 1.0000 |

## Hold-out Metrics (QRETIX v2 Fused)

| Metric | Value |
|--------|-------|
| AUPRC (hold-out) | 1.0000 |
| P@10 (hold-out)  | 1.0000 |
| P@50 (hold-out)  | 1.0000 |
| NDCG@50 (hold-out) | 1.0000 |

## Learned Fusion Weights

| Score Source | Weight |
|--------------|--------|
| rf_holdout_score | 0.2000 |
| logreg_score | 0.2000 |
| rf_score | 0.2000 |
| gbm_ensemble_score | 0.2000 |
| text_score | 0.2000 |

## Top 20 Overall Candidates

KIF11, CDK7, AURKB, NUF2, CDC20, CDK1, CCNA2, BIRC5, MTOR, CCNB1, POLE, KIF23, MYC, AURKA, TPX2, PPP2CA, CTCF, EIF4E, NIP7, BRD4

## Top 10 Hold-out Recoveries
(genes held out from training, correctly ranked highly)

CDK1, TPX2, RRM2, CENPN, WEE1, BCL2L1, GAPDH, MCM7, PIK3CA, ATM