"""
encode_fetcher.py
─────────────────
Phase 1.3: ENCODE epigenomics acquisition.

Searches for thyroid and anaplastic thyroid biosamples.
Expected finding: no ATC biosamples — normal thyroid only.
Sparsity is logged explicitly, not treated as null result.

Output: data/raw/encode/ + data/provenance/encode_thyroid.json
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import requests

from src.provenance import make_record, write as prov_write, record_path

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw" / "encode"

ENCODE_API = "https://www.encodeproject.org"
SEARCH_PARAMS_BASE = {
    "type": "Experiment",
    "organism.scientific_name": "Homo sapiens",
    "format": "json",
    "limit": 200,
    "field": "accession,biosample_ontology.term_name,assay_title,status,files.href",
}


def _search(biosample_term: str) -> list[dict]:
    params = {**SEARCH_PARAMS_BASE, "biosample_ontology.term_name": biosample_term}
    try:
        resp = requests.get(
            f"{ENCODE_API}/search/",
            params=params,
            headers={"Accept": "application/json"},
            timeout=60,
        )
        resp.raise_for_status()
        return resp.json().get("@graph", [])
    except Exception as exc:
        log.error("ENCODE search for '%s' failed: %s", biosample_term, exc)
        return []


def run() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    if record_path("encode_thyroid").exists():
        log.info("ENCODE provenance already exists.")
        return

    atc_results    = _search("anaplastic thyroid carcinoma")
    thyroid_results = _search("thyroid gland")

    log.info("ENCODE: ATC biosamples found: %d", len(atc_results))
    log.info("ENCODE: Normal thyroid biosamples found: %d", len(thyroid_results))

    summary = {
        "atc_biosample_count": len(atc_results),
        "thyroid_normal_count": len(thyroid_results),
        "atc_experiments": [e.get("accession") for e in atc_results],
        "note": (
            "ATC-specific ENCODE biosamples: NONE FOUND (as expected). "
            "Normal thyroid data available as reference only. "
            "This sparsity is logged explicitly — absence of ATC ENCODE data is a "
            "documented gap, not a null result."
        ),
    }
    (RAW_DIR / "encode_search_summary.json").write_text(json.dumps(summary, indent=2))

    record = make_record(
        source="ENCODE",
        accession="encode_thyroid",
        url="https://www.encodeproject.org/help/rest-api/",
        sample_count=len(thyroid_results),
        modality="epigenomics_atac_h3k27ac",
        license_="CC BY 4.0",
        raw_files=["encode_search_summary.json"],
        checksums={},
        notes=(
            f"ATC biosamples: {len(atc_results)} (NONE). "
            f"Normal thyroid biosamples: {len(thyroid_results)}. "
            "ENCODE does not have ATC-specific epigenomic data. "
            "Normal thyroid used as reference baseline only."
        ),
        is_processed=False,
        atc_specific=False,
        data_gap="No ATC-specific ENCODE biosamples available.",
    )
    prov_write("encode_thyroid", record)
