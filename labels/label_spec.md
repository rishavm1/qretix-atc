# Label Specification — Phase 4

## Label Hierarchy

| Stage | Definition | Source |
|-------|-----------|--------|
| 0 | No meaningful ATC evidence | Default |
| 1 | Preclinical evidence (in vitro/in vivo ATC) | PubMed NLP |
| 2 | Experimental therapeutic evidence (perturbation → phenotype) | PubMed NLP + DepMap |
| 3 | Clinical investigation (Phase 1-2 ATC trial) | ClinicalTrials.gov + OT |
| 4 | Validated / approved drug target in ATC | Curated + ClinicalTrials.gov |

## Assignment Logic

1. **Stage 4/3**: Curated known targets + clinical trial NLP extraction
2. **Stage 3/2/1**: Open Targets ATC association score thresholds (≥0.6/0.3/0.1)
3. **Stage 2/1**: PubMed NLP — stage2 if ≥2 perturbation-language papers; stage1 if ≥3 ATC mentions
4. **Stage 2**: DepMap selective ATC dependency (Chronos < -0.5)
5. **Hard negatives**: Known cancer genes not linked to ATC (from literature review)
6. **True negatives**: Stage 0 genes with ≥2 negative-result ATC papers

## Temporal Metadata

Every labelled gene has `earliest_evidence_date` (YYYY-MM-DD or YYYY-01-01).
This is used by Phase 7 temporal retrospective validation to filter labels by cutoff date.

## Leakage Firewall

- Known driver flag (BRAF, TP53 etc.) is kept in `label_notes` only
- Feature table has NO column encoding 'is known driver' — isolated here
- DepMap dependency used as evidence proxy ONLY (not as a feature value in labelling loop)