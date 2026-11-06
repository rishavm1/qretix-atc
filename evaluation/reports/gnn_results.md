# GraphSAGE GNN Ranker Results

## Architecture
- 2-layer GraphSAGE with mean aggregation
- Hidden dim: 128 → Output dim: 64
- Edges: STRING PPI (combined_score ≥ 400)
- Loss: PU bipartite ranking (margin=0.5)
- Epochs: 100 | LR=0.001

## Metrics

| Metric | Value |
|--------|-------|
| AUPRC | 0.551 |
| P@10 | 0.9 |
| P@20 | 0.75 |
| P@50 | 0.56 |
| NDCG@50 | 0.6457 |
| n_epochs | 100 |
| hidden_dim | 128 |
| out_dim | 64 |
| n_genes | 18743 |
| loss_final | 0.099432 |

## Top 20 Ranked Genes

BRAF, TERT, RET, SUPT5H, PIK3CA, RPS27A, TP53, UBB, UBA52, SUPT16H, DHX16, GTF2F1, TAF5, UBC, DHX8, MED25, H4C6, NCBP1, TAF6, ERCC2