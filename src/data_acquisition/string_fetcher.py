"""
string_fetcher.py
─────────────────
Phase 1.5: STRING PPI network acquisition.

Downloads:
  - 9606.protein.links.v12.0.txt.gz  (full human PPI, combined score)
  - 9606.protein.info.v12.0.txt.gz   (protein metadata: preferred name, annotation)

Filters to combined_score >= 400 during Phase 3 feature engineering.
Protein IDs are Ensembl protein IDs (ENSP) — mapped to gene IDs via gene_universe.

Output: data/raw/string/
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import requests

from src.provenance import make_record, write as prov_write, record_path

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw" / "string"

STRING_FILES = {
    "9606.protein.links.v12.0.txt.gz":
        "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz",
    "9606.protein.info.v12.0.txt.gz":
        "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz",
}


def _download(url: str, dest: Path, label: str) -> bool:
    if dest.exists():
        log.info("%s already downloaded.", label)
        return True
    log.info("Downloading %s (may be large) …", label)
    try:
        resp = requests.get(url, timeout=600, stream=True)
        resp.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                f.write(chunk)
        log.info("Saved %s → %.1f MB", dest.name, dest.stat().st_size / 1e6)
        return True
    except Exception as exc:
        log.error("Failed to download %s: %s", label, exc)
        return False


def run() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    success = {}
    for fname, url in STRING_FILES.items():
        success[fname] = _download(url, RAW_DIR / fname, fname)
        time.sleep(0.5)

    if not record_path("string_v12").exists():
        record = make_record(
            source="STRING v12.0",
            accession="string_v12",
            url="https://string-db.org",
            sample_count=0,
            modality="protein_protein_interaction",
            license_="CC BY 4.0",
            raw_files=list(STRING_FILES.keys()),
            checksums={},
            notes=(
                "Full human PPI network (taxon 9606). "
                "combined_score >= 400 filter applied in Phase 3. "
                "Protein IDs are ENSP — mapped to ENSG via gene_universe."
            ),
            is_processed=False,
        )
        prov_write("string_v12", record)
