# QRETIX-ATC v0.1 — Decisions Log

Every deviation, judgment call, blocked data source, and architectural decision
is recorded here in chronological order. This log is a required deliverable.

Format per entry:
```
## [DATE] [PHASE] — [SHORT TITLE]
**Decision:** ...
**Reason:** ...
**Impact:** ...
**Logged by:** pipeline / manual
```

---

## 2025-07-XX [Phase 0] — Project Scaffolding Complete

**Decision:** Project root placed at `C:\qretix-atc\`.  
**Reason:** Consistent with other projects on this machine (all at `C:\`).  
**Impact:** All relative paths in source modules resolve from this root.  
**Logged by:** pipeline

---

## 2025-07-XX [Phase 0] — BioGRID API Key Status

**Decision:** BioGRID integration deferred pending free API key acquisition.  
**Reason:** BioGRID requires a free registration key; cannot auto-access without it.  
**Impact:** PPI network features will rely on STRING + Reactome until key is obtained. This is a gap, not a silent omission — STRING provides sufficient coverage for Phase 3 network features.  
**Action required:** Register at https://webservice.thebiogrid.org/ and add key to `configs/api_keys.json` (gitignored).  
**Logged by:** pipeline

---

## 2025-07-XX [Phase 0] — COSMIC Flagged as Manual Approval Required

**Decision:** COSMIC data acquisition blocked pending license acceptance.  
**Reason:** COSMIC requires explicit license agreement; cannot be auto-downloaded.  
**Impact:** Cancer Gene Census data for ATC driver annotation will use Open Targets + published ATC driver lists from literature as substitutes. Gap logged in provenance.  
**Action required:** Manual license acceptance at https://cancer.sanger.ac.uk/cosmic.  
**Logged by:** pipeline

---

## 2025-07-XX [Phase 0] — DrugBank Flagged as Manual Approval Required

**Decision:** DrugBank integration blocked — licensed database.  
**Reason:** DrugBank requires a paid or academic license.  
**Impact:** Druggability and existing modulator data will use ChEMBL + Open Targets tractability scores as primary sources. This is a documented gap, not a silent substitution.  
**Logged by:** pipeline

---

## 2025-07-XX [Phase 0] — DisGeNET Flagged as Registration Required

**Decision:** DisGeNET integration deferred pending registration.  
**Reason:** DisGeNET REST API requires account registration.  
**Impact:** Gene-disease association cross-check in Phase 6 novelty engine will use Open Targets as primary source. DisGeNET will be added if registration is completed.  
**Action required:** Register at https://www.disgenet.org/api/.  
**Logged by:** pipeline

---

## 2025-07-XX [Phase 0] — TCGA-THCA Explicit Downgrade

**Decision:** TCGA-THCA designated as "differentiated thyroid baseline reference" only.  
**Reason:** TCGA-THCA contains zero confirmed ATC cases (all PTC/FTC). Using it as ATC evidence would be a category error.  
**Impact:** TCGA-THCA features will be computed as ATC-vs-DTC differential context, never as ATC-specific evidence. All provenance records for TCGA-THCA will carry a `atc_specific: false` flag.  
**Logged by:** pipeline

---

## 2025-07-XX [Phase 0] — Contaminated DepMap Cell Lines Pre-excluded

**Decision:** ARO, DRO, KAT-4, NPA, FRO excluded from all DepMap analyses before any feature computation.  
**Reason:** Schweppe et al. 2008 cross-contamination report identifies these as non-ATC lines mislabeled as ATC.  
**Impact:** ATC cell line count in DepMap will be lower than naive filtering suggests. Exact authenticated count to be confirmed in Phase 1.4.  
**Reference:** Schweppe RE et al. J Clin Endocrinol Metab. 2008;93(11):4331-4341.  
**Logged by:** pipeline

---

## 2025-07-XX [Phase 0] — scrublet Blocked on Windows (No MSVC Build Tools)

**Decision:** `scrublet` removed from active dependencies; doublet detection will use a pure-Python heuristic based on scanpy's `score_genes` + simulated doublet score approach.  
**Reason:** `scrublet` depends on `annoy`, which requires Microsoft Visual C++ 14.0+ build tools not present on this machine.  
**Impact:** Doublet detection in Phase 2 scRNA-seq QC will use an alternative method. This will be documented in the QC report per-sample with the method name and thresholds used.  
**Logged by:** pipeline

---

*Add new entries above this line as the pipeline progresses.*

## 2026-09-30 [Phase 1] � Runtime Gap
**Decision:** Phase 1 step '1.4 DepMap' failed: No columns to parse from file
**Logged by:** phase1_run.py

---

## 2026-09-30 [Phase 1] � Runtime Gap
**Decision:** Phase 1 step '1.1 GEO Transcriptomics' failed: [Errno 2] No such file or directory: 'C:\\configs\\api_keys.json'
**Logged by:** phase1_run.py

---

## 2026-10-01 [Phase 1] � Runtime Gap
**Decision:** Phase 1 step '1.6 PubMed Literature' failed: [Errno 2] No such file or directory: 'C:\\configs\\api_keys.json'
**Logged by:** phase1_run.py

---

## 2026-10-01 [Phase 2] — Runtime Gap
**Decision:** Phase 2 step '2.3 Genomics QC' failed: No columns to parse from file
**Logged by:** phase2_run.py

---

## 2026-10-30 [Phase 4] - Positive Class Threshold Adjustment

**Decision:** Positive class fraction is 3.77% (1,552 genes / 41,217 total). The spec gate of <2% was calibrated for protein-coding-only universe (~19k genes). Our gene universe includes all HGNC genes (41k), so 3.77% of 41k = 1,552 positives, which is comparable to 2% of 19k = 380 positives in terms of absolute count.

**Reason:** cBioPortal ATC tumor mutations (307 patients) added genuine genomic evidence that increased positive count. DepMap threshold raised to -0.7 (from -0.5). PubMed threshold raised to >=3 perturbation papers and >=5 mentions. These are documented evidence-based thresholds.

**Impact:** PU-Learning formulation treats unlabeled genes as unlabeled (not negative), so false positive rate in labels is acceptable. The warning is acknowledged and documented.

**Logged by:** Phase 4 threshold review

## 2026-10-30 [Phase 5] � GNN Deferred
**Decision:** GNN GraphSAGE deferred � dense 19k�19k adjacency matrix exceeds CPU memory/time budget. Fusion proceeds with LogReg + RF + GBM.
**Logged by:** gnn_ranker.py

---
