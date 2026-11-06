# Feature Specification — Phase 3

Total genes: 41,217  |  Total features: 45

## Feature Groups

### genetic_evidence
- Genes with data: 41,217
- Columns:
  - `mut_freq` — 0.0% NaN
  - `amp_freq` — 0.0% NaN
  - `del_freq` — 0.0% NaN
  - `any_genetic_event` — 0.0% NaN

### tumor_cell_specificity
- Genes with data: 41,217
- Columns:
  - `scrna_malignant_mean` — 0.0% NaN
  - `scrna_nonmal_mean` — 0.0% NaN
  - `scrna_specificity` — 0.0% NaN
  - `bulk_atc_mean` — 56.9% NaN
  - `bulk_pdtc_mean` — 56.9% NaN
  - `bulk_atc_vs_pdtc_fc` — 56.9% NaN

### functional_dependency
- Genes with data: 17,657
- Columns:
  - `dep_mean_chronos` — 57.2% NaN
  - `dep_min_chronos` — 57.2% NaN
  - `dep_n_lines` — 57.2% NaN
  - `dep_essential_frac` — 57.2% NaN

### literature_features
- Genes with data: 41,217
- Columns:
  - `lit_atc_mentions` — 0.0% NaN
  - `lit_thyroid_mentions` — 0.0% NaN
  - `lit_atc_ratio` — 0.0% NaN
  - `lit_recency_score` — 0.0% NaN

### disease_relevance
- Genes with data: 41,237
- Columns:
  - `ot_atc_score` — 0.0% NaN
  - `ot_literature_score` — 0.0% NaN
  - `ot_somatic_mut_score` — 0.0% NaN
  - `ot_sm_approved` — 0.0% NaN
  - `ot_ab_approved` — 0.0% NaN
  - `ot_druggability_tier` — categorical

### network_relevance
- Genes with data: 41,217
- Columns:
  - `string_degree` — 0.0% NaN
  - `string_degree_hi` — 0.0% NaN
  - `reactome_pathway_count` — 0.0% NaN
  - `ppi_driver_proximity` — 99.9% NaN

### protein_structural
- Genes with data: 41,217
- Columns:
  - `uniprot_length` — 51.6% NaN
  - `uniprot_domain_count` — 0.0% NaN
  - `uniprot_has_structure` — 0.0% NaN
  - `uniprot_is_membrane` — 0.0% NaN
  - `uniprot_is_secreted` — 0.0% NaN
  - `uniprot_is_reviewed` — 0.0% NaN
  - `chembl_modulator_count` — 0.0% NaN
  - `chembl_approved_count` — 0.0% NaN
  - `chembl_clinical_count` — 0.0% NaN

### evidence_independence
- Genes with data: 41,217
- Columns:
  - `ev_n_modalities` — 0.0% NaN
  - `ev_genomics_flag` — 0.0% NaN
  - `ev_bulk_rna_flag` — 0.0% NaN
  - `ev_scrna_flag` — 0.0% NaN
  - `ev_depmap_flag` — 0.0% NaN
  - `ev_literature_flag` — 0.0% NaN
  - `ev_cross_modality` — 0.0% NaN
  - `ev_lit_n_atc_mentions` — 0.0% NaN

## Gap Summary

| Module | Status |
|--------|--------|
| genetic_evidence | COMPLETE — DepMap somatic mutations + CNA |
| tumor_cell_specificity | COMPLETE — GSE193581 scRNA + GSE76039 bulk |
| functional_dependency | COMPLETE — DepMap CRISPR Chronos (7 ATC lines) |
| literature_features | COMPLETE — PubMed ATC + thyroid corpora |
| disease_relevance | COMPLETE — Open Targets ATC associations + tractability |
| network_relevance | COMPLETE — STRING PPI + Reactome + driver proximity |
| protein_structural | COMPLETE — UniProt + ChEMBL drug mechanisms |
| evidence_independence | COMPLETE — Cross-modality evidence flags |

## Leakage Risk Notes

| Feature | Risk | Mitigation |
|---------|------|------------|
| lit_atc_mentions | TEMPORAL | Per-PMID dates stored; filtered in Phase 7 temporal splits |
| ot_atc_score | TEMPORAL | OT snapshot date logged in provenance |
| ev_literature_flag | TEMPORAL | Derived from literature corpus; filtered in Phase 7 |
| dep_mean_chronos | LOW | CRISPR screens don't leak label information |
| genetic_known_driver_flag | HIGH | Kept as isolated column, never merged into score |