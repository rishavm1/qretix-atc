# Logistic Regression (PU-aware) Results

## Formulation
- PU learning via cost-sensitive weighting (Elkan & Noto 2008)
- Positive class: Stage ≥ 1 | Unlabelled: Stage = 0 (NOT negative)
- Elastic net regularization | 5-fold cross-validation

## Metrics (5-fold CV)

| Metric | Value |
|--------|-------|
| AUPRC | 0.9394 |
| P@10 | 1.0 |
| P@20 | 1.0 |
| P@50 | 1.0 |
| NDCG@50 | 1.0 |
| n_positive | 1552 |
| n_unlabelled | 39665 |
| n_features | 47 |

## Top 20 Ranked Genes

ATR, POLE, BRD4, MTOR, XPO1, EIF1AX, TOP2A, CDK12, TRRAP, TPR, CLTC, CDK7, PLK1, CHD4, POLD1, RAC1, POLD3, MED12, RRM2, NUP98