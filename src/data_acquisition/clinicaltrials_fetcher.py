"""
clinicaltrials_fetcher.py
─────────────────────────
Phase 1.6: ClinicalTrials.gov acquisition.

Query: Anaplastic Thyroid Cancer
API:   ClinicalTrials.gov v2  https://clinicaltrials.gov/api/v2/studies

Every record stores: NCT ID, phase, interventions, status, study_start_date.
study_start_date is MANDATORY — required for Phase 7 temporal splits.

Output: literature/atc_clinical_trials.json
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
LIT_DIR = ROOT / "literature"

CT_API = "https://clinicaltrials.gov/api/v2/studies"
OUTPUT_PATH = LIT_DIR / "atc_clinical_trials.json"
PAGE_SIZE = 100


def _fetch_page(query: str, page_token: str | None) -> dict:
    params: dict = {
        "query.cond": query,
        "pageSize": PAGE_SIZE,
        "format": "json",
        "fields": (
            "NCTId,BriefTitle,OverallStatus,Phase,StudyType,"
            "InterventionName,InterventionType,"
            "StartDate,PrimaryCompletionDate,CompletionDate,"
            "LeadSponsorName,Condition"
        ),
    }
    if page_token:
        params["pageToken"] = page_token
    resp = requests.get(CT_API, params=params, timeout=60)
    resp.raise_for_status()
    return resp.json()


def _extract(study: dict) -> dict:
    proto = study.get("protocolSection", {})
    id_mod     = proto.get("identificationModule", {})
    status_mod = proto.get("statusModule", {})
    design_mod = proto.get("designModule", {})
    arms_mod   = proto.get("armsInterventionsModule", {})
    sponsor_mod= proto.get("sponsorCollaboratorsModule", {})
    cond_mod   = proto.get("conditionsModule", {})

    interventions = [
        {"name": iv.get("interventionName", ""), "type": iv.get("interventionType", "")}
        for iv in arms_mod.get("interventions", [])
    ]

    # study_start_date is mandatory for temporal splits
    start_date = (
        status_mod.get("startDateStruct", {}).get("date", "") or
        status_mod.get("startDate", "")
    )

    return {
        "nct_id":           id_mod.get("nctId", ""),
        "title":            id_mod.get("briefTitle", ""),
        "status":           status_mod.get("overallStatus", ""),
        "phase":            design_mod.get("phases", []),
        "study_type":       design_mod.get("studyType", ""),
        "interventions":    interventions,
        "conditions":       cond_mod.get("conditions", []),
        "lead_sponsor":     sponsor_mod.get("leadSponsor", {}).get("name", ""),
        "study_start_date": start_date,
        "primary_completion": status_mod.get("primaryCompletionDateStruct", {}).get("date", ""),
    }


def run() -> None:
    LIT_DIR.mkdir(parents=True, exist_ok=True)

    if OUTPUT_PATH.exists():
        existing = json.loads(OUTPUT_PATH.read_text())
        log.info("atc_clinical_trials.json already exists (%d trials).", len(existing))
        return

    query = "Anaplastic Thyroid Cancer"
    log.info("Fetching ClinicalTrials.gov: '%s'", query)

    trials: list[dict] = []
    page_token: str | None = None

    while True:
        try:
            data = _fetch_page(query, page_token)
        except Exception as exc:
            log.error("ClinicalTrials.gov fetch failed: %s", exc)
            break

        studies = data.get("studies", [])
        for s in studies:
            trials.append(_extract(s))

        page_token = data.get("nextPageToken")
        log.info("Fetched %d trials so far …", len(trials))

        if not page_token:
            break
        time.sleep(0.5)

    OUTPUT_PATH.write_text(json.dumps(trials, indent=2))
    log.info("Clinical trials written: %d records → %s", len(trials), OUTPUT_PATH)

    # Validate: every record must have study_start_date
    missing = sum(1 for t in trials if not t.get("study_start_date"))
    if missing:
        log.warning("%d trials missing study_start_date.", missing)

    # Provenance
    if not record_path("clinicaltrials_atc").exists():
        record = make_record(
            source="ClinicalTrials.gov",
            accession="clinicaltrials_atc",
            url="https://clinicaltrials.gov/api/v2/studies",
            sample_count=len(trials),
            modality="clinical_trials",
            license_="Public domain (US government)",
            raw_files=["atc_clinical_trials.json"],
            checksums={},
            notes=f"Query: '{query}'. {len(trials)} trials retrieved. "
                  "study_start_date stored per record for Phase 7 temporal splits.",
            is_processed=False,
        )
        prov_write("clinicaltrials_atc", record)
