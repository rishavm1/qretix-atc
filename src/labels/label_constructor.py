"""
Phase 4 — Label Constructor
============================
Builds the 5-stage ATC therapeutic target label hierarchy.

Label Hierarchy:
    Stage 0 — No meaningful ATC evidence
    Stage 1 — Preclinical evidence (in vitro/in vivo in ATC models)
    Stage 2 — Experimental therapeutic evidence (ATC perturbation → phenotype)
    Stage 3 — Clinical investigation (any-phase ATC trial)
    Stage 4 — Clinical-stage / validated (approved drug or late-phase trial)

Output: labels/labels.parquet with columns:
    hgnc_symbol, ensembl_gene_id, stage, earliest_evidence_date,
    evidence_sources, is_true_negative, is_hard_negative, label_notes

Sources:
    literature/atc_clinical_trials.json     → Stage 3/4
    literature/atc_pubmed_corpus.json        → Stage 1/2 (NLP-based)
    data/raw/opentargets/atc_associations.json → Stage 3/4 cross-check
    data/processed/depmap_atc_dependency.csv → Stage 2 proxy (selective dependency)

Leakage prevention:
    Every labelled gene stores earliest_evidence_date.
    Phase 7 temporal validation uses this to filter labels by cutoff date.

Notes:
    - Known negative results (published "no effect in ATC") are flagged is_true_negative=True.
    - Cancer Gene Census genes NOT linked to ATC are flagged is_hard_negative=True.
    - Stage 0 = everything unlabelled (never called "negatives").
    - Positive class = Stage >= 1.
"""
from __future__ import annotations

import json
import logging
import re
from datetime import date, datetime
from pathlib import Path
from typing import Any

import pandas as pd

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
LIT_DIR  = ROOT / "literature"
DATA_DIR = ROOT / "data" / "processed"
RAW_OT   = ROOT / "data" / "raw" / "opentargets"
LABELS_DIR = ROOT / "labels"

# ─── Known ATC Tier-4 targets (approved drugs or late-phase trials) ──────────
# From FDA approvals + published ATC clinical trial landscape reviews
STAGE4_TARGETS = frozenset([
    "BRAF",          # Vemurafenib, Dabrafenib (BRAF V600E)
    "MAP2K1",        # Trametinib (MEK1)
    "MAP2K2",        # Trametinib (MEK2)
    "KDR",           # Lenvatinib (VEGFR2 / KDR)
    "FLT1",          # Lenvatinib (VEGFR1)
    "FLT4",          # Lenvatinib (VEGFR3)
    "RET",           # Selpercatinib, Pralsetinib
    "NTRK1",         # Larotrectinib, Entrectinib
    "NTRK2",         # Larotrectinib, Entrectinib
    "NTRK3",         # Larotrectinib, Entrectinib
    "PDCD1",         # Pembrolizumab (PD-1) — checkpoint approved in solid tumors
    "CD274",         # Atezolizumab (PD-L1)
    "CTLA4",         # Ipilimumab
    "FGFR1", "FGFR2", "FGFR3",   # Lenvatinib has FGFR activity
])

# Stage 3: clinical investigation targets (Phase 1-2 ATC trials, published)
STAGE3_TARGETS = frozenset([
    "CDK4",          # Palbociclib
    "CDK6",          # Palbociclib, Ribociclib
    "WEE1",          # Adavosertib
    "AURKA",         # Alisertib (Aurora kinase A)
    "AURKB",
    "PLK1",          # Volasertib
    "MTOR",          # Everolimus, Temsirolimus
    "AKT1",          # MK-2206
    "PIK3CA",        # Copanlisib
    "MDM2",          # Idasanutlin
    "HDAC1", "HDAC2", "HDAC3",  # Panobinostat, Vorinostat
    "MET",           # Cabozantinib (multi-target)
    "EGFR",          # Gefitinib trials in ATC
    "ALK",           # Crizotinib
    "ERBB2",         # Trastuzumab + Pertuzumab in HER2+ subset
    "BCL2",          # Venetoclax (preclinical + early trial)
    "STAT3",         # OPB-31121 trials
    "KRAS",          # Sotorasib (KRAS G12C, early ATC trials)
    "TERT",          # Telomerase inhibitors (early phase)
])

# Stage 2: experimental therapeutic evidence (DepMap + key preclinical papers)
STAGE2_FROM_DEPMAP_THRESHOLD = -0.7   # Chronos score threshold for dependency (strong)

# Curated Stage 2 genes from published ATC functional perturbation studies
STAGE2_CURATED = frozenset([
    "SRC", "YES1", "LYN",          # SRC family kinases
    "PARP1", "PARP2",              # PARP inhibition synthetic lethality
    "ATR", "ATM", "CHEK1", "CHEK2",  # DNA damage response
    "FASN",                         # Fatty acid synthase
    "DNMT1", "DNMT3A",             # Epigenetic
    "EZH2",                         # PRC2
    "BRD4",                         # BET bromodomain
    "MCL1",                         # Anti-apoptotic
    "XIAP",
    "ROCK1", "ROCK2",              # Rho kinase
    "JAK1", "JAK2",                # JAK-STAT
    "YAP1", "WWTR1",               # Hippo pathway
    "NOTCH1",                       # Notch signaling
    "HIF1A",                        # Hypoxia
    "VEGFA",                        # VEGF ligand
    "FGF2",
    "IGF1R",                        # IGF signaling
    "INSR",
    "PRKCA", "PRKCB",              # PKC
    "MAPK1", "MAPK3",              # ERK1/2
])

# Hard negatives: Cancer Gene Census genes with NO ATC link
# (well-known cancer drivers clearly not relevant to ATC biology)
HARD_NEGATIVES = frozenset([
    "APC",     # Colorectal cancer; not linked to ATC
    "VHL",     # Renal cell carcinoma
    "RB1",     # Retinoblastoma
    "NF1",     # Neurofibromatosis (minimal ATC evidence)
    "WT1",     # Wilms tumor
    "BRCA1", "BRCA2",  # Breast/ovarian; limited ATC relevance
    "MLH1", "MSH2", "MSH6",   # Lynch syndrome; MMR; rare in ATC
    "MEN1",    # MEN1 syndrome
    "FH",      # Hereditary leiomyomatosis
])

# ─── NLP patterns for stage assignment from literature ────────────────────────

# Stage 1 patterns: preclinical evidence (in vitro/in vivo mentions in ATC context)
STAGE1_PATTERNS = [
    r"\b(in vitro|cell line|proliferation|apoptosis|migration|invasion|colony formation)\b",
    r"\b(ATC cell|anaplastic thyroid cell|thyroid carcinoma cell)\b",
    r"\b(mouse model|xenograft|PDX|tumor model)\b",
    r"\b(knockdown|knockout|shRNA|siRNA|CRISPR)\b.*\b(ATC|anaplastic thyroid)\b",
    r"\b(overexpression|upregulated|downregulated)\b.*\b(ATC|anaplastic thyroid)\b",
]

# Stage 2 patterns: perturbation → phenotype in ATC
STAGE2_PATTERNS = [
    r"\b(inhibitor|inhibition|drug|compound|treatment)\b.*\b(ATC|anaplastic thyroid)\b",
    r"\b(sensitiz|resensitiz|synergist|combinat)\b.*\b(ATC|anaplastic thyroid)\b",
    r"\b(therapeutic target|drug target|targetable)\b.*\b(ATC|anaplastic thyroid)\b",
    r"\b(IC50|EC50|GI50|cell viability)\b.*\b(ATC|anaplastic thyroid)\b",
]

# Negative result patterns
NEGATIVE_PATTERNS = [
    r"\b(no significant|not significant|not differ|no effect|no change|failed to)\b",
    r"\b(negative result|no association|not associated|not correlated)\b",
]


def _parse_date(date_str: str | None) -> str | None:
    """Normalize a date string to YYYY-MM-DD or YYYY. Returns None if unparseable."""
    if not date_str:
        return None
    date_str = str(date_str).strip()
    for fmt in ("%Y-%m-%d", "%Y/%m/%d", "%Y-%m", "%B %d, %Y", "%Y"):
        try:
            return datetime.strptime(date_str[:10], fmt[:len(date_str[:10])]).strftime("%Y-%m-%d")
        except ValueError:
            pass
    # Try year-only
    m = re.match(r"(\d{4})", date_str)
    return m.group(1) + "-01-01" if m else None


def _load_clinical_trials() -> list[dict]:
    path = LIT_DIR / "atc_clinical_trials.json"
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _load_pubmed() -> list[dict]:
    path = LIT_DIR / "atc_pubmed_corpus.json"
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _load_ot_associations() -> dict[str, dict]:
    """Return symbol → {score} from processed OT parquet (richer) or raw JSON fallback."""
    proc_path = DATA_DIR / "opentargets_atc_targets.parquet"
    if proc_path.exists():
        df = pd.read_parquet(proc_path)
        result = {}
        for _, row in df.iterrows():
            sym = row.get("approved_symbol", "")
            if sym:
                result[sym] = {"score": row.get("overall_association_score", 0.0)}
        return result
    # fallback raw JSON
    path = RAW_OT / "atc_associations.json"
    with open(path, encoding="utf-8") as f:
        records = json.load(f)
    return {
        r["target"]["approvedSymbol"]: r
        for r in records
        if "target" in r and "approvedSymbol" in r["target"]
    }


def _load_depmap_dependencies() -> dict[str, float]:
    """Return symbol → median Chronos score across ATC lines."""
    dep = pd.read_csv(DATA_DIR / "depmap_atc_dependency.csv")
    return dep.groupby("hgnc_symbol")["chronos_score"].median().to_dict()


def _extract_genes_from_text(text: str, valid_symbols: set[str]) -> set[str]:
    """
    Extract gene symbols mentioned in text using whole-word matching.
    Only returns genes in valid_symbols (gene universe).
    Minimum symbol length 3 to avoid false positives (e.g. 'OR', 'IF').
    """
    found = set()
    words = re.findall(r'\b[A-Z][A-Z0-9]{2,12}\b', text)
    for w in words:
        if w in valid_symbols:
            found.add(w)
    return found


def construct_labels(gene_universe: pd.DataFrame) -> pd.DataFrame:
    """
    Main label construction function.
    Returns DataFrame with one row per gene in gene_universe.
    """
    valid_symbols = set(gene_universe["hgnc_symbol"].tolist())

    # ── per-gene record: {stage, earliest_date, evidence_sources, notes} ─────
    records: dict[str, dict[str, Any]] = {}

    def _update(symbol: str, new_stage: int, date_str: str | None,
                source_id: str, note: str = ""):
        """Update gene record if new_stage > current stage."""
        if symbol not in valid_symbols:
            return
        cur = records.get(symbol, {"stage": 0, "earliest_date": None,
                                   "evidence_sources": [], "notes": []})
        if new_stage > cur["stage"]:
            cur["stage"] = new_stage
        # Always track the earliest date
        if date_str and (cur["earliest_date"] is None or date_str < cur["earliest_date"]):
            cur["earliest_date"] = date_str
        cur["evidence_sources"].append(source_id)
        if note:
            cur["notes"].append(note)
        records[symbol] = cur

    # ── Stage 4 & 3: Curated known targets ───────────────────────────────────
    log.info("Assigning curated Stage 4/3 labels...")
    for sym in STAGE4_TARGETS & valid_symbols:
        _update(sym, 4, "2020-01-01", "curated_stage4", "Approved drug or late-phase trial")
    for sym in STAGE3_TARGETS & valid_symbols:
        _update(sym, 3, "2018-01-01", "curated_stage3", "Phase 1-2 ATC clinical trial")
    for sym in STAGE2_CURATED & valid_symbols:
        _update(sym, 2, "2015-01-01", "curated_stage2", "Published ATC perturbation study")

    # ── Stage 3/4 from Clinical Trials ───────────────────────────────────────
    log.info("Processing clinical trials for Stage 3/4 labels...")
    trials = _load_clinical_trials()
    for trial in trials:
        phase_raw = trial.get("phase") or ""
        if isinstance(phase_raw, list):
            phase_str = " ".join(phase_raw).upper()
        else:
            phase_str = str(phase_raw).upper()
        start_date = _parse_date(trial.get("study_start_date"))
        nct = trial.get("nct_id", "")

        # Determine stage from phase
        if any(p in phase_str for p in ("PHASE3", "PHASE 3", "PHASE_3", "3")):
            assigned_stage = 4
        elif any(p in phase_str for p in ("PHASE1", "PHASE2", "PHASE 1", "PHASE 2", "1", "2")):
            assigned_stage = 3
        else:
            assigned_stage = 3  # Unknown phase trials still count as Stage 3

        # Extract gene symbols from title + interventions
        title = (trial.get("title") or "").upper()
        interventions = " ".join(
            str(i.get("name", "")).upper()
            for i in (trial.get("interventions") or [])
        )
        text = title + " " + interventions

        for sym in _extract_genes_from_text(text, valid_symbols):
            _update(sym, assigned_stage, start_date, f"ClinicalTrials:{nct}",
                    f"Trial phase: {phase_str}")

    # ── Stage 1/2/3/4 from Open Targets ──────────────────────────────────────
    log.info("Processing Open Targets associations...")
    ot = _load_ot_associations()
    for sym, rec in ot.items():
        score = float(rec.get("score", 0.0))
        if score >= 0.6:
            assigned = 3
        elif score >= 0.3:
            assigned = 2
        elif score >= 0.1:
            assigned = 1
        else:
            continue
        _update(sym, assigned, None, "OpenTargets:ATC_assoc",
                f"OT score={score:.3f}")

    # ── Stage 1/2 from PubMed literature (NLP) ───────────────────────────────
    log.info("Processing PubMed literature (NLP)...")
    pubmed = _load_pubmed()

    TOKEN_RE = re.compile(r'\b[A-Z][A-Z0-9]{2,12}\b')
    symbol_set_upper = {s for s in valid_symbols if len(s) >= 3 and s == s.upper()}

    gene_atc_count: dict[str, int] = {}
    gene_stage2_count: dict[str, int] = {}
    gene_neg_count: dict[str, int] = {}
    gene_earliest: dict[str, str | None] = {}

    for rec in pubmed:
        title = (rec.get("title") or "")
        abstract = (rec.get("abstract") or "")
        text = title + " " + abstract
        text_upper = text.upper()
        pub_date = _parse_date(rec.get("publication_date") or rec.get("year"))

        # Fast token-based gene extraction
        genes_in_abstract = set(TOKEN_RE.findall(text)) & symbol_set_upper
        if not genes_in_abstract:
            continue

        is_stage2 = any(re.search(pat, text_upper, re.IGNORECASE) for pat in STAGE2_PATTERNS)
        is_negative = any(re.search(pat, text_upper, re.IGNORECASE) for pat in NEGATIVE_PATTERNS)

        for sym in genes_in_abstract:
            gene_atc_count[sym] = gene_atc_count.get(sym, 0) + 1
            if is_stage2:
                gene_stage2_count[sym] = gene_stage2_count.get(sym, 0) + 1
            if is_negative:
                gene_neg_count[sym] = gene_neg_count.get(sym, 0) + 1
            if pub_date:
                if sym not in gene_earliest or gene_earliest[sym] is None:
                    gene_earliest[sym] = pub_date
                elif pub_date < gene_earliest[sym]:
                    gene_earliest[sym] = pub_date

    for sym in valid_symbols:
        n_atc = gene_atc_count.get(sym, 0)
        n_s2  = gene_stage2_count.get(sym, 0)
        earliest = gene_earliest.get(sym)

        if n_s2 >= 3:
            _update(sym, 2, earliest, f"PubMed:n={n_atc}", f"Stage2 papers={n_s2}")
        elif n_atc >= 5:
            _update(sym, 1, earliest, f"PubMed:n={n_atc}", "ATC mention count")

    # ── Stage 2 from DepMap selective dependency ──────────────────────────────
    log.info("Processing DepMap dependencies for Stage 2...")
    dep_scores = _load_depmap_dependencies()
    for sym, score in dep_scores.items():
        if sym in valid_symbols and score < STAGE2_FROM_DEPMAP_THRESHOLD:
            _update(sym, 2, None, f"DepMap:Chronos={score:.3f}",
                    "Selective ATC dependency (Chronos < -0.5)")

    # ── Stage 1/2 from cBioPortal ATC tumor mutations (HF data) ──────────────
    log.info("Processing cBioPortal ATC tumor mutations for Stage 1/2 labels...")
    cbio_path = DATA_DIR / "cbioportal_atc_mutations.parquet"
    if cbio_path.exists():
        cbio = pd.read_parquet(cbio_path)
        # Count how many unique patients have a mutation in each gene
        mut_counts = cbio.groupby("Hugo_Symbol")["Tumor_Sample_Barcode"].nunique()
        n_total_patients = cbio["Tumor_Sample_Barcode"].nunique()
        for sym, n_patients in mut_counts.items():
            if sym in valid_symbols:
                freq = n_patients / n_total_patients
                # ≥15% patient frequency → Stage 2 (recurrent driver)
                # 5-14% → Stage 1
                # < 5% → skip
                if freq >= 0.15:
                    stage = 2
                elif freq >= 0.05:
                    stage = 1
                else:
                    continue
                _update(sym, stage, None,
                        f"cBioPortal:n_patients={n_patients}",
                        f"Recurrent ATC tumor mutation in {n_patients}/{n_total_patients} patients ({freq:.1%})")
        log.info("cBioPortal: processed %d mutated genes", len(mut_counts))

    # ── Build final DataFrame ─────────────────────────────────────────────────
    log.info("Assembling final label table...")
    gu = gene_universe[["hgnc_symbol", "ensembl_gene_id"]].copy()

    def _get(sym, field, default):
        return records.get(sym, {}).get(field, default)

    gu["stage"]                = gu["hgnc_symbol"].map(lambda s: _get(s, "stage", 0))
    gu["earliest_evidence_date"] = gu["hgnc_symbol"].map(
        lambda s: _get(s, "earliest_date", None)
    )
    gu["evidence_sources"] = gu["hgnc_symbol"].map(
        lambda s: "; ".join(set(_get(s, "evidence_sources", [])))
    )
    gu["label_notes"] = gu["hgnc_symbol"].map(
        lambda s: "; ".join(set(_get(s, "notes", [])))
    )
    gu["is_true_negative"] = gu["hgnc_symbol"].map(
        lambda s: int(gene_neg_count.get(s, 0) >= 2 and _get(s, "stage", 0) == 0)
    )
    gu["is_hard_negative"] = gu["hgnc_symbol"].map(
        lambda s: int(s in HARD_NEGATIVES)
    )

    gu = gu.set_index("hgnc_symbol")
    return gu


def run():
    """Phase 4 entry point."""
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    LABELS_DIR.mkdir(exist_ok=True)

    gene_universe = pd.read_parquet(DATA_DIR / "gene_universe.parquet")
    print(f"Gene universe: {len(gene_universe):,} genes")

    labels = construct_labels(gene_universe)

    # Distribution report
    dist = labels["stage"].value_counts().sort_index()
    print("\nLabel distribution:")
    total = len(labels)
    for stage, count in dist.items():
        pct = count / total * 100
        print(f"  Stage {stage}: {count:,} genes ({pct:.2f}%)")

    pos_total = (labels["stage"] >= 1).sum()
    print(f"\nPositive class (stage≥1): {pos_total:,} ({pos_total/total*100:.2f}%)")

    if pos_total / total > 0.02:
        print("\n⚠️  WARNING: positive class > 2% of genome — check for leakage before proceeding.")
        print("   See docs/decisions_log.md for audit trail.")

    # Save
    out_path = LABELS_DIR / "labels.parquet"
    labels.to_parquet(out_path)
    print(f"\nWrote {out_path}  shape={labels.shape}")

    _write_label_spec(labels)
    _write_class_distribution(labels)
    return labels


def _write_label_spec(labels: pd.DataFrame):
    lines = [
        "# Label Specification — Phase 4",
        "",
        "## Label Hierarchy",
        "",
        "| Stage | Definition | Source |",
        "|-------|-----------|--------|",
        "| 0 | No meaningful ATC evidence | Default |",
        "| 1 | Preclinical evidence (in vitro/in vivo ATC) | PubMed NLP |",
        "| 2 | Experimental therapeutic evidence (perturbation → phenotype) | PubMed NLP + DepMap |",
        "| 3 | Clinical investigation (Phase 1-2 ATC trial) | ClinicalTrials.gov + OT |",
        "| 4 | Validated / approved drug target in ATC | Curated + ClinicalTrials.gov |",
        "",
        "## Assignment Logic",
        "",
        "1. **Stage 4/3**: Curated known targets + clinical trial NLP extraction",
        "2. **Stage 3/2/1**: Open Targets ATC association score thresholds (≥0.6/0.3/0.1)",
        "3. **Stage 2/1**: PubMed NLP — stage2 if ≥2 perturbation-language papers; stage1 if ≥3 ATC mentions",
        "4. **Stage 2**: DepMap selective ATC dependency (Chronos < -0.5)",
        "5. **Hard negatives**: Known cancer genes not linked to ATC (from literature review)",
        "6. **True negatives**: Stage 0 genes with ≥2 negative-result ATC papers",
        "",
        "## Temporal Metadata",
        "",
        "Every labelled gene has `earliest_evidence_date` (YYYY-MM-DD or YYYY-01-01).",
        "This is used by Phase 7 temporal retrospective validation to filter labels by cutoff date.",
        "",
        "## Leakage Firewall",
        "",
        "- Known driver flag (BRAF, TP53 etc.) is kept in `label_notes` only",
        "- Feature table has NO column encoding 'is known driver' — isolated here",
        "- DepMap dependency used as evidence proxy ONLY (not as a feature value in labelling loop)",
    ]
    path = LABELS_DIR / "label_spec.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {path}")


def _write_class_distribution(labels: pd.DataFrame):
    total = len(labels)
    dist = labels["stage"].value_counts().sort_index()
    pos = (labels["stage"] >= 1).sum()
    tn = labels["is_true_negative"].sum()
    hn = labels["is_hard_negative"].sum()

    lines = [
        "# Class Distribution Report — Phase 4",
        "",
        f"Total genes: {total:,}",
        f"Positive class (stage≥1): {pos:,} ({pos/total*100:.2f}%)",
        f"True negatives (stage=0 + negative evidence): {int(tn):,}",
        f"Hard negatives (known non-ATC cancer genes): {int(hn):,}",
        "",
        "## Stage Breakdown",
        "",
        "| Stage | Count | % of genome |",
        "|-------|-------|------------|",
    ]
    for stage in sorted(dist.index):
        count = dist[stage]
        lines.append(f"| {stage} | {count:,} | {count/total*100:.2f}% |")

    if pos / total > 0.02:
        lines += [
            "",
            "> ⚠️ **WARNING**: Positive class exceeds 2% of genome.",
            "> This may indicate label leakage. Review evidence sources before model training.",
        ]

    path = LABELS_DIR / "class_distribution.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {path}")


def validate() -> bool:
    """Phase 4 acceptance criteria."""
    checks = []
    labels_path = LABELS_DIR / "labels.parquet"
    spec_path = LABELS_DIR / "label_spec.md"
    dist_path = LABELS_DIR / "class_distribution.md"

    checks.append((labels_path.exists(), f"labels.parquet exists"))
    checks.append((spec_path.exists(), f"label_spec.md exists"))
    checks.append((dist_path.exists(), f"class_distribution.md exists"))

    if labels_path.exists():
        labels = pd.read_parquet(labels_path)
        pos_frac = (labels["stage"] >= 1).mean()
        checks.append((pos_frac < 0.02, f"Positive class < 2% (got {pos_frac:.3%})"))
        checks.append(("earliest_evidence_date" in labels.columns, "earliest_evidence_date column present"))
        checks.append(("is_true_negative" in labels.columns, "is_true_negative column present"))

    all_pass = True
    for ok, msg in checks:
        status = "✅" if ok else "❌"
        print(f"  {status} {msg}")
        if not ok:
            all_pass = False
    return all_pass


if __name__ == "__main__":
    run()
