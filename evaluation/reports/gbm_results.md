# GBM Ranker Results (XGBoost / LightGBM / CatBoost)

## Formulation
- Ranking-based PU learning
- XGBoost: rank:pairwise objective
- LightGBM: LambdaMART (lambdarank objective)
- CatBoost: YetiRank objective
- Positives (stage≥1) = relevant; Unlabelled (stage=0) = non-relevant
- Never treats unlabelled as explicit negatives in loss computation

## Metrics

| Model | AUPRC | P@10 | P@20 | P@50 | NDCG@50 |
|-------|-------|------|------|------|---------|
| XGBoost_rank:pairwise | 0.8674 | 1.0 | 1.0 | 1.0 | 1.0 |
| LightGBM_lambdarank | 0.9752 | 1.0 | 1.0 | 1.0 | 1.0 |

## Top 20 Ranked Genes (Ensemble Score)

AARS1, WDR3, VARS1, DHX16, CPSF2, NOP2, POLR1A, POLR3B, IARS1, POLR1B, DHX8, CPSF3, UTP20, DHX37, UTP4, LARS1, UTP15, WDR77, XAB2, HMGCS1