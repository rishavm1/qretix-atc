"""
reactome_fetcher.py
───────────────────
Phase 1.5: Reactome pathway acquisition.

Downloads:
  - ReactomePathways.gmt          (gene sets per pathway)
  - ReactomePathways.txt          (pathway ID → name → species)
  - Ensembl2Reactome_All_Levels.txt (gene → pathway mapping, all levels)

Output: data/raw/reactome/
"""

from __future__ import annotations

import logging
from pathlib import Path

import requests

from src.provenance import make_record, write as prov_write, record_path

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw" / "reactome"

REACTOME_FILES = {
    "ReactomePathways.gmt":
        "https://reactome.org/download/current/ReactomePathways.gmt.zip",
    "ReactomePathways.txt":
        "https://reactome.org/download/current/ReactomePathways.txt",
    "Ensembl2Reactome_All_Levels.txt":
        "https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt",
}


def _download(url: str, dest: Path, label: str) -> bool:
    if dest.exists():
        log.info("%s already downloaded.", label)
        return True
    log.info("Downloading %s …", label)
    try:
        resp = requests.get(url, timeout=120, stream=True)
        resp.raise_for_status()
        dest.write_bytes(resp.content)
        log.info("Saved %s → %.2f MB", dest.name, dest.stat().st_size / 1e6)
        return True
    except Exception as exc:
        log.error("Failed to download %s: %s", label, exc)
        return False


def run() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    for fname, url in REACTOME_FILES.items():
        _download(url, RAW_DIR / fname, fname)

    if not record_path("reactome_current").exists():
        record = make_record(
            source="Reactome",
            accession="reactome_current",
            url="https://reactome.org/download/current/",
            sample_count=0,
            modality="pathway_annotation",
            license_="CC BY 4.0",
            raw_files=list(REACTOME_FILES.keys()),
            checksums={},
            notes=(
                "Reactome current release. Human pathways. "
                "Ensembl2Reactome_All_Levels.txt used for gene→pathway mapping in Phase 3. "
                "Mapped to gene_universe via Ensembl Gene ID."
            ),
            is_processed=False,
        )
        prov_write("reactome_current", record)
