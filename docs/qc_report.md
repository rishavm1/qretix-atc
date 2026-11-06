# QRETIX-ATC v0.1 — QC Report

_Generated: 2026-10-01T23:45:03.412902+00:00_

## Summary

## GSE193581 (scRNA-seq)

- Cells after QC: **66016**
- Genes: 38224
- Cell types: 13
- Malignant cells: 18056
- Non-malignant cells: 47960

### Per-sample QC thresholds

| Sample | Cells in | Cells out | % removed | min_genes | max_genes | max_mito% | doublets |
|--------|----------|-----------|-----------|-----------|-----------|-----------|----------|
| PTC01 | 5543 | 5367 | 3.2% | 200 | 5711 | 20.0 | 55 |
| PTC02 | 2855 | 2330 | 18.4% | 910 | 2982 | 20.0 | 24 |
| PTC03 | 996 | 921 | 7.5% | 200 | 4815 | 27.1 | 10 |
| NORM03 | 808 | 610 | 24.5% | 1011 | 2575 | 20.0 | 7 |
| PTC04 | 2796 | 2750 | 1.6% | 200 | 7150 | 29.8 | 28 |
| PTC05 | 2742 | 2573 | 6.2% | 200 | 8000 | 21.7 | 26 |
| PTC06 | 3109 | 2902 | 6.7% | 200 | 7965 | 20.2 | 30 |
| PTC07 | 6928 | 6858 | 1.0% | 200 | 7072 | 20.0 | 70 |
| NORM07 | 1812 | 1503 | 17.1% | 741 | 3529 | 20.2 | 16 |
| ATC08 | 8994 | 8233 | 8.5% | 200 | 3581 | 21.7 | 84 |
| ATC09 | 2091 | 2020 | 3.4% | 200 | 7619 | 20.0 | 21 |
| ATC10 | 6525 | 6057 | 7.2% | 701 | 3797 | 20.0 | 62 |
| ATC11 | 3226 | 2867 | 11.1% | 201 | 3541 | 20.0 | 29 |
| ATC12 | 1789 | 1668 | 6.8% | 200 | 8000 | 20.0 | 17 |
| ATC13 | 5017 | 4441 | 11.5% | 200 | 8000 | 21.6 | 45 |
| ATC14 | 2174 | 2119 | 2.5% | 200 | 6788 | 28.7 | 22 |
| ATC15 | 2049 | 1907 | 6.9% | 514 | 3634 | 21.9 | 20 |
| ATC17 | 2071 | 1934 | 6.6% | 391 | 3327 | 20.0 | 20 |
| ATC18 | 4417 | 3894 | 11.8% | 644 | 3468 | 20.0 | 40 |
| NORM18 | 2476 | 2122 | 14.3% | 587 | 2595 | 20.0 | 22 |
| NORM19 | 1170 | 948 | 19.0% | 379 | 1979 | 20.0 | 10 |
| NORM20 | 1989 | 1782 | 10.4% | 575 | 3383 | 27.8 | 18 |
| NORM21 | 254 | 210 | 17.3% | 1779 | 3919 | 20.0 | 3 |

## GSE76039 (Bulk Microarray)

- Samples in: **37**
- Samples out: **37** (no sample-level exclusions for bulk array)
- % removed: 0%
- Probes raw: 54675
- Genes after probe collapse: **22880**
- Probe collapse method: max_mean (keep probe with highest mean expression per gene)
- Log transform: log2(x+1) if max>100 else none
- Tumor type distribution: {'ATC': 20, 'PDTC': 17}
- Platform: Affymetrix_HGU133Plus2_GPL570
- Note: Pre-called expression values from GEO series matrix — not raw CEL files.

## Genomics (DepMap ATC Cell Lines)

### Mutations
- Rows: 1293
- Genes: 1155
- Source: DepMap precalled (not raw VCFs) — cBioPortal ATC cohorts unavailable (404)

### CNA
- Rows: 269325
- Genes: 38590
- Source: DepMap precalled (not raw VCFs) — cBioPortal ATC cohorts unavailable (404)


## Batch Tagging

All processed files audited for platform/batch columns by batch_tagger.py.

## Documented Data Gaps

| Gap | Reason | Impact |
|-----|--------|--------|
| cBioPortal GATCI 2024 | Study ID thca_gatci_cellrep_2024 returned 404 | Using DepMap somatic mutations as proxy |
| cBioPortal MSKCC 2016 | Study ID thca_mskcc_2016 returned 404 | Same proxy |
| scrublet doublet detection | annoy build requires MSVC (absent) | Using count_per_gene_99pct heuristic |
| ENCODE ATC epigenomics | No ATC biosamples in ENCODE | Normal thyroid reference only |

---
_All thresholds logged per-dataset in data/provenance/*.json_