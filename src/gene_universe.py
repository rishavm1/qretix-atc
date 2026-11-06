"""
gene_universe.py
────────────────
Central gene ID mapping layer. Must be built ONCE before any downstream module
runs. Every join across datasets goes through this table — no ad-hoc remapping.

Primary key : Ensembl Gene ID (ENSG...)
Mapped IDs  : HGNC Symbol, HGNC ID, UniProt (Swiss-Prot), Entrez Gene ID
Scope       : protein-coding genes only (~19,000–20,000 rows)

Output      : data/processed/gene_universe.parquet
"""

from __future__ import annotations

import hashlib
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

log = logging.getLogger(__name__)

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
PROVENANCE_DIR = ROOT / "data" / "provenance"

HGNC_URL = (
    "https://www.genenames.org/cgi-bin/download/custom"
    "?col=gd_hgnc_id&col=gd_app_sym&col=gd_pub_ensembl_id"
    "&col=gd_pub_eg_id&col=md_prot_id"
    "&status=Approved&locus_type=gene+with+protein+product"
    "&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
)

OUTPUT_PATH = PROCESSED_DIR / "gene_universe.parquet"
PROVENANCE_PATH = PROVENANCE_DIR / "gene_universe.json"
RAW_PATH = RAW_DIR / "hgnc_protein_coding.tsv"

REQUIRED_COLS = {
    "HGNC ID": "hgnc_id",
    "Approved symbol": "hgnc_symbol",
    "Ensembl gene ID": "ensembl_gene_id",
    "NCBI Gene ID": "entrez_gene_id",
    "UniProt ID(supplied by UniProt)": "uniprot_id",
}

MIN_GENE_COUNT = 19_000


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def download_hgnc(force: bool = False) -> Path:
    """Download HGNC protein-coding gene table to data/raw/."""
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    if RAW_PATH.exists() and not force:
        log.info("HGNC raw file already present, skipping download.")
        return RAW_PATH

    log.info("Downloading HGNC protein-coding gene set …")
    resp = requests.get(HGNC_URL, timeout=60)
    resp.raise_for_status()
    RAW_PATH.write_bytes(resp.content)
    log.info("Saved → %s", RAW_PATH)
    return RAW_PATH


def build(force: bool = False) -> pd.DataFrame:
    """
    Build gene_universe.parquet.

    Returns the DataFrame. Raises ValueError if row count < MIN_GENE_COUNT
    (acceptance criterion: ≥19,000 protein-coding genes).
    """
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    PROVENANCE_DIR.mkdir(parents=True, exist_ok=True)

    if OUTPUT_PATH.exists() and not force:
        log.info("gene_universe.parquet already exists, loading.")
        return pd.read_parquet(OUTPUT_PATH)

    raw_path = download_hgnc(force=force)

    df = pd.read_csv(raw_path, sep="\t", dtype=str)

    # Rename to canonical column names
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"HGNC download missing expected columns: {missing}")

    df = df[list(REQUIRED_COLS.keys())].rename(columns=REQUIRED_COLS)

    # Drop rows with no Ensembl ID — these cannot be used as primary keys
    before = len(df)
    df = df.dropna(subset=["ensembl_gene_id"]).reset_index(drop=True)
    dropped = before - len(df)
    if dropped:
        log.warning("Dropped %d rows with missing Ensembl Gene ID.", dropped)

    # Deduplicate on Ensembl ID (keep first — HGNC is already curated)
    df = df.drop_duplicates(subset=["ensembl_gene_id"]).reset_index(drop=True)

    if len(df) < MIN_GENE_COUNT:
        raise ValueError(
            f"Gene universe has only {len(df)} rows — expected ≥{MIN_GENE_COUNT}. "
            "Check HGNC download."
        )

    df.to_parquet(OUTPUT_PATH, index=False)
    log.info("gene_universe.parquet written: %d genes → %s", len(df), OUTPUT_PATH)

    # ── Provenance record ──────────────────────────────────────────────────────
    prov = {
        "source": "HGNC (HUGO Gene Nomenclature Committee)",
        "accession": "hgnc_protein_coding_complete_set",
        "url": HGNC_URL,
        "date_pulled": datetime.now(timezone.utc).isoformat(),
        "sample_count": len(df),
        "modality": "gene_id_mapping",
        "license": "Creative Commons Attribution 4.0 (CC BY 4.0)",
        "raw_files": [str(RAW_PATH.name)],
        "checksums": {"hgnc_protein_coding.tsv": f"sha256:{_sha256(raw_path)}"},
        "notes": (
            "Protein-coding genes only (locus_type=gene with protein product). "
            f"Rows with missing Ensembl ID dropped: {dropped}."
        ),
        "qc_thresholds_applied": None,
        "is_processed": True,
    }
    PROVENANCE_PATH.write_text(json.dumps(prov, indent=2))
    log.info("Provenance written → %s", PROVENANCE_PATH)

    return df


def load() -> pd.DataFrame:
    """Load the pre-built gene universe. Raises FileNotFoundError if not built."""
    if not OUTPUT_PATH.exists():
        raise FileNotFoundError(
            f"{OUTPUT_PATH} not found. Run gene_universe.build() first."
        )
    return pd.read_parquet(OUTPUT_PATH)


def map_symbols_to_ensembl(symbols: list[str]) -> dict[str, str]:
    """Return {hgnc_symbol: ensembl_gene_id} for the given symbol list."""
    df = load()
    subset = df[df["hgnc_symbol"].isin(symbols)][["hgnc_symbol", "ensembl_gene_id"]]
    return dict(zip(subset["hgnc_symbol"], subset["ensembl_gene_id"]))


def validate() -> bool:
    """Phase 0 acceptance check: parquet exists with ≥19,000 rows and 5 columns."""
    if not OUTPUT_PATH.exists():
        log.error("FAIL: gene_universe.parquet missing.")
        return False
    df = pd.read_parquet(OUTPUT_PATH)
    if len(df) < MIN_GENE_COUNT:
        log.error("FAIL: only %d genes (need ≥%d).", len(df), MIN_GENE_COUNT)
        return False
    expected_cols = set(REQUIRED_COLS.values())
    missing = expected_cols - set(df.columns)
    if missing:
        log.error("FAIL: missing columns %s.", missing)
        return False
    log.info("PASS: gene_universe.parquet — %d genes, columns OK.", len(df))
    return True


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    build()
    validate()
