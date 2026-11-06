# Random Forest (PU-Bagging) Results

## Formulation
- PU bagging (Mordelet & Vert 2014)
- 50 bags × 200 trees
- Each bag: all positives + random unlabelled sample (size=n_pos)
- Final score = mean predicted probability across bags

## Metrics

| Metric | Value |
|--------|-------|
| AUPRC | 0.9938 |
| P@10 | 1.0 |
| P@20 | 1.0 |
| P@50 | 1.0 |
| NDCG@50 | 1.0 |
| n_bags | 50 |
| n_trees_per_bag | 200 |
| n_positive | 1552 |
| n_unlabelled | 39665 |
| n_features | 48 |

## Top 20 Ranked Genes

AQR, EIF3G, PSMD11, PSMD12, MED17, PSMA3, PSMA2, PSMA4, PSMC1, MPHOSPH10, RPF1, SNRNP70, CCT2, RPA2, COPS4, LAS1L, PRPF19, SMC3, SF1, RIOK2