"""
Improvement 3 — Temporal Evidence Date Backfill
=================================================
Problem: Only 166 / 41,217 genes have earliest_evidence_date.
The temporal validation is nearly blind because it can't split
labels by cutoff year for most genes.

Strategy:
  1. For PubMed-sourced labels: use earliest paper date from corpus
  2. For DepMap-sourced labels: use DepMap release date (2024-Q4 = 2024-10-01)
  3. For cBioPortal-sourced labels: use publication date of GATCI 2024 paper
  4. For curated Stage 4/3: use known FDA approval / trial start dates
  5. For OpenTargets-sourced labels: use OT snapshot date from provenance

This dramatically improves temporal split coverage.
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
LABELS_PATH = ROOT / "labels" / "labels.parquet"
LIT_PATH    = ROOT / "literature" / "atc_pubmed_corpus.json"

# Known dates for curated sources
DEPMAP_24Q4_DATE   = "2024-10-01"
CBIOPORTAL_GATCI   = "2024-03-15"  # GATCI 2024 Cell Reports paper
OT_SNAPSHOT_DATE   = "2024-01-01"  # approximate OT data pull date

# Stage 4: FDA approval dates for ATC
KNOWN_APPROVAL_DATES = {
    "BRAF":   "2018-05-04",  # Dabrafenib+Trametinib FDA approval for ATC
    "MAP2K1": "2018-05-04",
    "MAP2K2": "2018-05-04",
    "KDR":    "2015-02-13",  # Lenvatinib FDA approval thyroid
    "FLT1":   "2015-02-13",
    "FLT4":   "2015-02-13",
    "FGFR1":  "2015-02-13",
    "FGFR2":  "2015-02-13",
    "FGFR3":  "2015-02-13",
    "RET":    "2020-05-08",  # Selpercatinib FDA
    "NTRK1":  "2018-11-26",  # Larotrectinib FDA
    "NTRK2":  "2018-11-26",
    "NTRK3":  "2018-11-26",
    "PDCD1":  "2017-06-16",  # Pembrolizumab FDA solid tumors
    "CD274":  "2016-05-18",  # Atezolizumab
    "CTLA4":  "2011-03-25",  # Ipilimumab
}

# Stage 3: first ATC clinical trial start dates (approximate)
KNOWN_TRIAL_DATES = {
    "CDK4":   "2015-01-01",
    "CDK6":   "2015-01-01",
    "WEE1":   "2016-01-01",
    "AURKA":  "2013-01-01",
    "AURKB":  "2013-01-01",
    "PLK1":   "2014-01-01",
    "MTOR":   "2010-01-01",
    "AKT1":   "2011-01-01",
    "PIK3CA": "2013-01-01",
    "MDM2":   "2016-01-01",
    "HDAC1":  "2012-01-01",
    "HDAC2":  "2012-01-01",
    "HDAC3":  "2012-01-01",
    "MET":    "2012-01-01",
    "EGFR":   "2003-01-01",
    "ALK":    "2011-01-01",
    "ERBB2":  "2013-01-01",
    "BCL2":   "2016-01-01",
    "STAT3":  "2014-01-01",
    "KRAS":   "2020-01-01",
    "TERT":   "2017-01-01",
}


def _build_pubmed_earliest_dates(valid_symbols: set) -> dict:
    """Return {symbol: earliest_date_str} from PubMed corpus."""
    with open(LIT_PATH, encoding="utf-8") as f:
        records = json.load(f)

    TOKEN_RE = re.compile(r'\b[A-Z][A-Z0-9]{2,14}\b')
    symbol_set = {s for s in valid_symbols if len(s) >= 3 and s == s.upper()}

    earliest: dict[str, str] = {}
    for rec in records:
        text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
        raw_date = rec.get("publication_date") or rec.get("year")
        if not raw_date:
            continue
        date_str = str(raw_date)[:10]
        if len(date_str) == 4:
            date_str = date_str + "-01-01"

        hits = set(TOKEN_RE.findall(text)) & symbol_set
        for sym in hits:
            if sym not in earliest or date_str < earliest[sym]:
                earliest[sym] = date_str

    return earliest


def backfill_temporal_dates() -> pd.DataFrame:
    labels = pd.read_parquet(LABELS_PATH)
    valid_symbols = set(labels.index)

    log.info("Building PubMed earliest dates for %d symbols...", len(valid_symbols))
    pubmed_dates = _build_pubmed_earliest_dates(valid_symbols)
    log.info("PubMed dates found for %d genes", len(pubmed_dates))

    filled = 0
    for sym in labels.index:
        current_date = labels.at[sym, "earliest_evidence_date"]
        sources = str(labels.at[sym, "evidence_sources"]) or ""
        stage = int(labels.at[sym, "stage"])

        candidates = []

        # Known curated dates
        if sym in KNOWN_APPROVAL_DATES:
            candidates.append(KNOWN_APPROVAL_DATES[sym])
        if sym in KNOWN_TRIAL_DATES:
            candidates.append(KNOWN_TRIAL_DATES[sym])

        # PubMed date
        if sym in pubmed_dates:
            candidates.append(pubmed_dates[sym])

        # Source-based dates
        if "cBioPortal" in sources:
            candidates.append(CBIOPORTAL_GATCI)
        if "DepMap" in sources:
            candidates.append(DEPMAP_24Q4_DATE)
        if "OpenTargets" in sources:
            candidates.append(OT_SNAPSHOT_DATE)

        if not candidates:
            continue

        best_date = min(candidates)

        if current_date is None or (best_date < str(current_date)):
            labels.at[sym, "earliest_evidence_date"] = best_date
            filled += 1

    coverage = labels["earliest_evidence_date"].notna().sum()
    log.info("Temporal date coverage: %d / %d genes (filled %d)", coverage, len(labels), filled)

    labels.to_parquet(LABELS_PATH)
    log.info("Saved updated labels → %s", LABELS_PATH)
    return labels


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    backfill_temporal_dates()
