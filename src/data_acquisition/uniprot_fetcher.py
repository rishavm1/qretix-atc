"""
uniprot_fetcher.py
──────────────────
Phase 1.5: UniProt protein feature acquisition.

Fetches for all protein-coding genes (from gene_universe):
  - Sequence length, domain count (Pfam), subcellular localization,
    PTMs, existing PDB structures, Swiss-Prot reviewed status.

Uses UniProt REST API streaming endpoint for bulk retrieval.
Output: data/raw/uniprot/uniprot_human_proteome.tsv
"""

from __future__ import annotations

import logging
from pathlib import Path

import requests

from src.provenance import make_record, write as prov_write, record_path

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw" / "uniprot"

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
OUTPUT_FILE = RAW_DIR / "uniprot_human_proteome.tsv"

FIELDS = ",".join([
    "accession", "gene_names", "gene_primary",
    "length", "organism_name",
    "cc_subcellular_location", "ft_domain",
    "xref_pdb", "reviewed",
    "protein_existence", "annotation_score",
    "xref_ensembl",
])

QUERY = "organism_id:9606 AND reviewed:true"


def run() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    if OUTPUT_FILE.exists():
        log.info("uniprot_human_proteome.tsv already exists.")
    else:
        log.info("Downloading UniProt human Swiss-Prot proteome …")
        params = {
            "query": QUERY,
            "format": "tsv",
            "fields": FIELDS,
            "compressed": "false",
        }
        try:
            resp = requests.get(UNIPROT_STREAM, params=params, timeout=300, stream=True)
            resp.raise_for_status()
            with open(OUTPUT_FILE, "wb") as f:
                for chunk in resp.iter_content(chunk_size=1024 * 512):
                    f.write(chunk)
            log.info("UniProt proteome saved → %s (%.1f MB)",
                     OUTPUT_FILE, OUTPUT_FILE.stat().st_size / 1e6)
        except Exception as exc:
            log.error("UniProt download failed: %s", exc)
            return

    if not record_path("uniprot_human_swissprot").exists():
        record = make_record(
            source="UniProt Swiss-Prot",
            accession="uniprot_human_swissprot",
            url="https://rest.uniprot.org/uniprotkb/stream",
            sample_count=0,
            modality="protein_annotation",
            license_="CC BY 4.0",
            raw_files=[OUTPUT_FILE.name],
            checksums={},
            notes=(
                "Human Swiss-Prot reviewed entries (organism_id:9606, reviewed:true). "
                f"Fields: {FIELDS}. "
                "Mapped to gene_universe via xref_ensembl in Phase 3."
            ),
            is_processed=False,
        )
        prov_write("uniprot_human_swissprot", record)
