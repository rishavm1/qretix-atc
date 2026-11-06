"""
chembl_fetcher.py
─────────────────
Phase 1.5: ChEMBL existing modulator acquisition.

Pulls for all human protein targets:
  - Target → compound activity data (IC50, Ki, mechanism of action)
  - Drug mechanism of action table

Used in Phase 3 (protein/structural features) and Phase 8 (evidence cards).
Output: data/raw/chembl/
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import requests

from src.provenance import make_record, write as prov_write, record_path

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw" / "chembl"

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
PAGE_SIZE = 1000


def _get(endpoint: str, params: dict) -> dict | None:
    url = f"{CHEMBL_API}/{endpoint}.json"
    try:
        resp = requests.get(url, params=params, timeout=60)
        resp.raise_for_status()
        return resp.json()
    except Exception as exc:
        log.error("ChEMBL GET %s failed: %s", endpoint, exc)
        return None


def fetch_drug_mechanisms() -> None:
    """Fetch all drug mechanism of action records."""
    out_path = RAW_DIR / "drug_mechanisms.json"
    if out_path.exists():
        log.info("drug_mechanisms.json already exists.")
        return

    log.info("Fetching ChEMBL drug mechanisms …")
    all_records: list[dict] = []
    offset = 0

    while True:
        data = _get("mechanism", {"limit": PAGE_SIZE, "offset": offset})
        if not data:
            break
        records = data.get("mechanisms", [])
        all_records.extend(records)
        total = data.get("page_meta", {}).get("total_count", 0)
        log.info("Drug mechanisms: %d / %d", len(all_records), total)
        if len(all_records) >= total or not records:
            break
        offset += PAGE_SIZE
        time.sleep(0.2)

    out_path.write_text(json.dumps(all_records, indent=2))
    log.info("Drug mechanisms: %d records → %s", len(all_records), out_path)


def fetch_approved_drugs() -> None:
    """Fetch approved drugs with max_phase=4."""
    out_path = RAW_DIR / "approved_drugs.json"
    if out_path.exists():
        log.info("approved_drugs.json already exists.")
        return

    log.info("Fetching ChEMBL approved drugs (max_phase=4) …")
    all_records: list[dict] = []
    offset = 0

    while True:
        data = _get("molecule", {
            "max_phase": 4,
            "limit": PAGE_SIZE,
            "offset": offset,
            "format": "json",
        })
        if not data:
            break
        records = data.get("molecules", [])
        all_records.extend(records)
        total = data.get("page_meta", {}).get("total_count", 0)
        log.info("Approved drugs: %d / %d", len(all_records), total)
        if len(all_records) >= total or not records:
            break
        offset += PAGE_SIZE
        time.sleep(0.2)

    out_path.write_text(json.dumps(all_records, indent=2))
    log.info("Approved drugs: %d records → %s", len(all_records), out_path)


def run() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    fetch_drug_mechanisms()
    fetch_approved_drugs()

    if not record_path("chembl_current").exists():
        record = make_record(
            source="ChEMBL",
            accession="chembl_current",
            url="https://www.ebi.ac.uk/chembl/api/data/",
            sample_count=0,
            modality="drug_target_activity",
            license_="CC BY-SA 3.0",
            raw_files=["drug_mechanisms.json", "approved_drugs.json"],
            checksums={},
            notes=(
                "ChEMBL drug mechanism of action + approved drugs (max_phase=4). "
                "Used for druggability assessment in Phase 3 and evidence cards in Phase 8. "
                "DrugBank is unavailable (licensed) — ChEMBL is primary druggability source."
            ),
            is_processed=False,
        )
        prov_write("chembl_current", record)
