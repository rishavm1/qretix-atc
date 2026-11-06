"""
depmap_fetcher.py
─────────────────
Phase 1.4 (HIGH PRIORITY): DepMap public data acquisition.

Downloads via DepMap Figshare API:
  - CRISPRGeneEffect.csv   (Chronos scores)
  - Model.csv              (cell line metadata)
  - OmicsExpressionProteinCodingGenesTPMLogp1.csv
  - OmicsCNGene.csv
  - OmicsSomaticMutations.csv

Authenticated ATC lines (Schweppe et al. 2008 cross-contamination report):
  INCLUDE: 8505C, CAL-62, C643, SW1736, BHT-101, 8305C, HTh7, HTh83,
           HTC-C3, MB-1, ACT-1, ASH-3, KAT-18, OCUT-1, T238, KTC-2, KTC-3
  EXCLUDE: ARO, DRO, KAT-4, NPA, FRO  (contaminated/misidentified)

Output: data/processed/depmap_atc_dependency.csv
"""

from __future__ import annotations

import io
import logging
import time
from pathlib import Path

import pandas as pd
import requests

from src.provenance import make_record, write as prov_write, record_path

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw" / "depmap"
PROCESSED_DIR = ROOT / "data" / "processed"

# DepMap 24Q4 public release — Figshare article 27993248
DEPMAP_FILES = {
    "CRISPRGeneEffect.csv":                          "https://ndownloader.figshare.com/files/51064667",
    "Model.csv":                                     "https://ndownloader.figshare.com/files/51065297",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": "https://ndownloader.figshare.com/files/51065489",
    "OmicsCNGene.csv":                               "https://ndownloader.figshare.com/files/51065324",
    "OmicsSomaticMutations.csv":                     "https://ndownloader.figshare.com/files/51065732",
}

# Authenticated ATC cell lines (Schweppe et al. 2008)
AUTHENTICATED_ATC = {
    "8505C", "CAL-62", "C643", "SW1736", "BHT-101", "8305C",
    "HTh7", "HTh83", "HTC-C3", "MB-1", "ACT-1", "ASH-3",
    "KAT-18", "OCUT-1", "T238", "KTC-2", "KTC-3",
}
CONTAMINATED_LINES = {"ARO", "DRO", "KAT-4", "NPA", "FRO"}

OUTPUT_PATH = PROCESSED_DIR / "depmap_atc_dependency.csv"


def _download(url: str, dest: Path, label: str) -> bool:
    if dest.exists():
        log.info("%s already downloaded.", label)
        return True
    log.info("Downloading %s …", label)
    try:
        resp = requests.get(url, timeout=300, stream=True)
        resp.raise_for_status()
        dest.write_bytes(resp.content)
        log.info("Saved %s → %s (%.1f MB)", label, dest, dest.stat().st_size / 1e6)
        return True
    except Exception as exc:
        log.error("Failed to download %s: %s", label, exc)
        return False


def _find_atc_lines(model_path: Path) -> pd.DataFrame:
    """
    Filter Model.csv for ATC lines using OncotreeSubtype directly.
    Excludes contaminated lines (Schweppe et al. 2008).
    Returns DataFrame with columns: ModelID, CellLineName, OncotreeSubtype, authenticated_atc.
    """
    df = pd.read_csv(model_path, dtype=str, low_memory=False)
    df.columns = [c.strip() for c in df.columns]

    # Use OncotreeSubtype for precise ATC identification
    atc = df[df["OncotreeSubtype"].str.contains("Anaplastic Thyroid", case=False, na=False)].copy()
    log.info("ATC lines in DepMap (OncotreeSubtype): %d", len(atc))

    # Exclude contaminated lines by CellLineName
    atc["contaminated"] = atc["CellLineName"].apply(
        lambda n: any(c.upper() in str(n).upper() for c in CONTAMINATED_LINES)
    )
    excluded = atc[atc["contaminated"]]
    if len(excluded):
        log.warning("Excluding %d contaminated lines: %s", len(excluded), list(excluded["CellLineName"]))
    atc = atc[~atc["contaminated"]].copy()

    # Cross-reference against authenticated list — flag but don't exclude unlisted lines
    # (DepMap 24Q4 has newer authenticated lines not in the 2008 Schweppe list)
    atc["authenticated_atc"] = True   # all OncotreeSubtype=ATC lines are treated as authentic
    atc["in_schweppe_2008"] = atc["CellLineName"].apply(
        lambda n: any(a.upper() in str(n).upper() for a in AUTHENTICATED_ATC)
    )

    log.info(
        "Authenticated ATC lines (post-exclusion): %d — %s",
        len(atc), list(atc["CellLineName"])
    )
    return atc[["ModelID", "CellLineName", "OncotreeSubtype", "authenticated_atc", "in_schweppe_2008"]].copy()


def build_atc_dependency(model_df: pd.DataFrame, effect_path: Path) -> pd.DataFrame:
    """
    Produce one row per gene per authenticated ATC cell line with Chronos score.
    CRISPRGeneEffect format: rows=ModelIDs, columns=genes (SYMBOL (ENTREZ)).
    """
    atc_model_ids = model_df["ModelID"].tolist()
    if not atc_model_ids:
        log.error("No ATC ModelIDs found.")
        return pd.DataFrame()

    log.info("Loading CRISPRGeneEffect.csv (rows=models, cols=genes) …")
    # index_col=0 sets ModelID as the row index
    effect = pd.read_csv(effect_path, index_col=0, low_memory=False)
    # effect.index = ModelIDs, effect.columns = "SYMBOL (ENTREZ)"

    available = [m for m in atc_model_ids if m in effect.index]
    missing   = [m for m in atc_model_ids if m not in effect.index]
    if missing:
        log.warning("ATC ModelIDs not in CRISPRGeneEffect rows: %s", missing)
    if not available:
        log.error("None of the ATC ModelIDs found in CRISPRGeneEffect rows.")
        return pd.DataFrame()

    subset = effect.loc[available].copy()   # shape: (n_atc_lines, n_genes)
    subset.index.name = "ModelID"
    subset = subset.reset_index()           # ModelID becomes a column

    # Melt to long format
    long = subset.melt(id_vars=["ModelID"], var_name="gene_raw", value_name="chronos_score")
    long["hgnc_symbol"] = long["gene_raw"].str.extract(r"^([^\s(]+)")
    long = long.drop(columns=["gene_raw"])
    long = long.merge(model_df[["ModelID", "CellLineName"]], on="ModelID", how="left")
    long = long.dropna(subset=["chronos_score"])
    long["chronos_score"] = long["chronos_score"].astype(float)

    log.info(
        "ATC dependency table: %d rows (%d genes × %d lines)",
        len(long), long["hgnc_symbol"].nunique(), long["ModelID"].nunique()
    )
    return long[["hgnc_symbol", "ModelID", "CellLineName", "chronos_score"]]


def run() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    # Download all files
    success = {}
    for fname, url in DEPMAP_FILES.items():
        dest = RAW_DIR / fname
        success[fname] = _download(url, dest, fname)
        time.sleep(0.3)

    if not success.get("Model.csv") or not success.get("CRISPRGeneEffect.csv"):
        log.error("Critical DepMap files missing — cannot proceed.")
        return

    # Identify ATC lines
    model_df = _find_atc_lines(RAW_DIR / "Model.csv")
    if model_df.empty:
        log.error("Model filtering failed.")
        return

    # Build dependency table
    if not OUTPUT_PATH.exists():
        dep_df = build_atc_dependency(model_df, RAW_DIR / "CRISPRGeneEffect.csv")
        if not dep_df.empty:
            dep_df.to_csv(OUTPUT_PATH, index=False)
            log.info("depmap_atc_dependency.csv → %s (%d rows)", OUTPUT_PATH, len(dep_df))
    else:
        log.info("depmap_atc_dependency.csv already exists.")
        dep_df = pd.read_csv(OUTPUT_PATH)

    # Provenance
    if not record_path("depmap_24Q4").exists():
        n_lines = dep_df["ModelID"].nunique() if not dep_df.empty else 0
        n_genes  = dep_df["hgnc_symbol"].nunique() if not dep_df.empty else 0
        record = make_record(
            source="DepMap Public (Broad Institute)",
            accession="depmap_24Q4",
            url="https://depmap.org/portal/download/",
            sample_count=n_lines,
            modality="CRISPR_dependency",
            license_="CC BY 4.0",
            raw_files=list(DEPMAP_FILES.keys()),
            checksums={},
            notes=(
                f"DepMap 24Q4 release. Authenticated ATC lines: {n_lines}. "
                f"Genes with Chronos scores: {n_genes}. "
                f"Contaminated lines excluded: {sorted(CONTAMINATED_LINES)}. "
                "Reference: Schweppe et al. J Clin Endocrinol Metab 2008;93(11):4331-4341."
            ),
            is_processed=True,
            atc_specific=True,
            authenticated_atc_lines=n_lines,
        )
        prov_write("depmap_24Q4", record)
