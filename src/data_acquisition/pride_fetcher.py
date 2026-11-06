"""
pride_fetcher.py
────────────────
Phase 1.3: PRIDE proteomics acquisition.

Searches for thyroid and anaplastic thyroid proteomics datasets.
Expected finding: sparse — a small number of ATC proteomics datasets.
Sparsity is logged explicitly, not treated as null result.

Output: data/raw/pride/ + data/provenance/pride_thyroid.json
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import requests

from src.provenance import make_record, write as prov_write, record_path

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw" / "pride"

PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v2"


def _search(keyword: str, page_size: int = 100) -> list[dict]:
    params = {
        "keyword": keyword,
        "pageSize": page_size,
        "page": 0,
        "sortDirection": "DESC",
        "sortConditions": "submissionDate",
    }
    try:
        resp = requests.get(
            f"{PRIDE_API}/projects",
            params=params,
            timeout=60,
        )
        resp.raise_for_status()
        return resp.json() if isinstance(resp.json(), list) else resp.json().get("_embedded", {}).get("projects", [])
    except Exception as exc:
        log.error("PRIDE search for '%s' failed: %s", keyword, exc)
        return []


def run() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    if record_path("pride_thyroid").exists():
        log.info("PRIDE provenance already exists.")
        return

    atc_results     = _search("anaplastic thyroid")
    thyroid_results = _search("thyroid cancer proteomics")

    log.info("PRIDE: ATC proteomics datasets found: %d", len(atc_results))
    log.info("PRIDE: Thyroid cancer proteomics datasets found: %d", len(thyroid_results))

    def _summarise(results: list[dict]) -> list[dict]:
        out = []
        for r in results:
            out.append({
                "accession": r.get("accession", r.get("projectAccession", "")),
                "title": r.get("title", ""),
                "submission_date": r.get("submissionDate", ""),
                "species": r.get("organisms", []),
            })
        return out

    summary = {
        "atc_datasets": _summarise(atc_results),
        "thyroid_datasets": _summarise(thyroid_results),
        "note": (
            f"ATC-specific PRIDE datasets: {len(atc_results)}. "
            f"Thyroid cancer proteomics datasets: {len(thyroid_results)}. "
            "Proteomics data for ATC is sparse — this is a documented gap. "
            "Known ATC proteomics datasets: PXD030304, PXD063991 (if present above)."
        ),
    }
    (RAW_DIR / "pride_search_summary.json").write_text(json.dumps(summary, indent=2))

    record = make_record(
        source="PRIDE Archive",
        accession="pride_thyroid",
        url="https://www.ebi.ac.uk/pride/ws/archive/v2/",
        sample_count=len(atc_results) + len(thyroid_results),
        modality="proteomics_phosphoproteomics",
        license_="Public (PRIDE submission agreement)",
        raw_files=["pride_search_summary.json"],
        checksums={},
        notes=(
            f"ATC proteomics datasets: {len(atc_results)}. "
            f"Thyroid cancer proteomics: {len(thyroid_results)}. "
            "Sparsity explicitly documented. "
            "Quantitation tables to be downloaded in Phase 2 if datasets are available."
        ),
        is_processed=False,
        atc_specific=len(atc_results) > 0,
        data_gap=f"ATC proteomics is sparse ({len(atc_results)} datasets found).",
    )
    prov_write("pride_thyroid", record)
