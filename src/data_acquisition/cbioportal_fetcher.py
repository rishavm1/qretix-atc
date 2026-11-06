"""
cbioportal_fetcher.py
─────────────────────
Phase 1.2: cBioPortal genomics acquisition.

Studies pulled:
  - thca_gatci_cellrep_2024  : 190 ATC tumors (CRITICAL — largest ATC cohort)
  - thca_mskcc_2016          : 33 ATC + 84 PDTC

Data pulled per study: mutations, CNA, structural variants, clinical data.

TCGA-THCA is NOT pulled here — it is a PTC reference, handled separately.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import pandas as pd
import requests

from src.provenance import make_record, write as prov_write, record_path

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw" / "cbioportal"

BASE_URL = "https://www.cbioportal.org/api"
HEADERS = {"Accept": "application/json", "Content-Type": "application/json"}

STUDIES = {
    "thca_gatci_cellrep_2024": {
        "notes": "GATCI 2024. 190 ATC tumors. Cell Reports 2024. Largest unified ATC genomic cohort.",
        "atc_specific": True,
        "expected_samples": 190,
    },
    "thca_mskcc_2016": {
        "notes": "Landa et al. JCI 2016. 33 ATC + 84 PDTC. Mutations + expression z-scores.",
        "atc_specific": True,
        "expected_samples": 117,
    },
}

MOLECULAR_PROFILES = {
    "mutations":   "_mutations",
    "cna":         "_gistic",
    "expression":  "_mrna_median_Zscores",
    "structural":  "_structural_variants",
}


def _get(endpoint: str, params: dict | None = None) -> dict | list | None:
    url = f"{BASE_URL}/{endpoint}"
    try:
        resp = requests.get(url, headers=HEADERS, params=params, timeout=60)
        resp.raise_for_status()
        return resp.json()
    except Exception as exc:
        log.error("cBioPortal GET %s failed: %s", endpoint, exc)
        return None


def _post(endpoint: str, payload: dict) -> dict | list | None:
    url = f"{BASE_URL}/{endpoint}"
    try:
        resp = requests.post(url, headers=HEADERS, json=payload, timeout=120)
        resp.raise_for_status()
        return resp.json()
    except Exception as exc:
        log.error("cBioPortal POST %s failed: %s", endpoint, exc)
        return None


def fetch_study_metadata(study_id: str) -> dict | None:
    return _get(f"studies/{study_id}")


def fetch_samples(study_id: str) -> list:
    result = _get(f"studies/{study_id}/samples", params={"pageSize": 10000})
    return result if isinstance(result, list) else []


def fetch_mutations(study_id: str, sample_ids: list[str]) -> pd.DataFrame:
    profile_id = f"{study_id}{MOLECULAR_PROFILES['mutations']}"
    payload = {
        "sampleIds": sample_ids,
        "molecularProfileId": profile_id,
    }
    data = _post("mutations/fetch", payload)
    if not data:
        return pd.DataFrame()
    df = pd.DataFrame(data)
    log.info("%s mutations: %d records", study_id, len(df))
    return df


def fetch_cna(study_id: str, sample_ids: list[str]) -> pd.DataFrame:
    profile_id = f"{study_id}{MOLECULAR_PROFILES['cna']}"
    payload = {
        "sampleIds": sample_ids,
        "molecularProfileId": profile_id,
    }
    data = _post("discrete-copy-number/fetch", payload)
    if not data:
        return pd.DataFrame()
    df = pd.DataFrame(data)
    log.info("%s CNA: %d records", study_id, len(df))
    return df


def fetch_clinical(study_id: str) -> pd.DataFrame:
    data = _get(f"studies/{study_id}/clinical-data", params={"clinicalDataType": "SAMPLE", "pageSize": 100000})
    if not data:
        return pd.DataFrame()
    df = pd.DataFrame(data)
    log.info("%s clinical data: %d records", study_id, len(df))
    return df


def fetch_structural_variants(study_id: str, sample_ids: list[str]) -> pd.DataFrame:
    profile_id = f"{study_id}{MOLECULAR_PROFILES['structural']}"
    payload = {
        "sampleIds": sample_ids,
        "molecularProfileId": profile_id,
    }
    data = _post("structural-variant/fetch", payload)
    if not data:
        return pd.DataFrame()
    df = pd.DataFrame(data)
    log.info("%s structural variants: %d records", study_id, len(df))
    return df


def run_study(study_id: str, meta: dict) -> None:
    study_dir = RAW_DIR / study_id
    study_dir.mkdir(parents=True, exist_ok=True)

    if record_path(study_id).exists():
        log.info("%s provenance already exists, skipping.", study_id)
        return

    log.info("── Fetching cBioPortal study: %s ──", study_id)

    # Verify study exists
    study_meta = fetch_study_metadata(study_id)
    if not study_meta:
        log.error("Study %s not found on cBioPortal.", study_id)
        _write_gap(study_id, meta, "Study not found on cBioPortal")
        return

    # Samples
    samples = fetch_samples(study_id)
    sample_ids = [s["sampleId"] for s in samples]
    log.info("%s: %d samples", study_id, len(sample_ids))

    if samples:
        pd.DataFrame(samples).to_csv(study_dir / "samples.csv", index=False)

    # Mutations
    time.sleep(0.3)
    mut_df = fetch_mutations(study_id, sample_ids)
    if not mut_df.empty:
        mut_df.to_csv(study_dir / "mutations.csv", index=False)

    # CNA
    time.sleep(0.3)
    cna_df = fetch_cna(study_id, sample_ids)
    if not cna_df.empty:
        cna_df.to_csv(study_dir / "cna.csv", index=False)

    # Clinical
    time.sleep(0.3)
    clin_df = fetch_clinical(study_id)
    if not clin_df.empty:
        clin_df.to_csv(study_dir / "clinical.csv", index=False)

    # Structural variants
    time.sleep(0.3)
    sv_df = fetch_structural_variants(study_id, sample_ids)
    if not sv_df.empty:
        sv_df.to_csv(study_dir / "structural_variants.csv", index=False)

    raw_files = [p.name for p in study_dir.iterdir() if p.is_file()]
    record = make_record(
        source="cBioPortal",
        accession=study_id,
        url=f"https://www.cbioportal.org/study/summary?id={study_id}",
        sample_count=len(sample_ids),
        modality="genomics_mutations_cna",
        license_="Public (cBioPortal open access)",
        raw_files=raw_files,
        checksums={},
        notes=meta["notes"] + " NOTE: pre-called mutations from cBioPortal, not raw VCFs.",
        is_processed=True,   # pre-called, not raw VCFs
        atc_specific=meta["atc_specific"],
    )
    prov_write(study_id, record)
    log.info("%s complete: %d samples, provenance written.", study_id, len(sample_ids))


def _write_gap(study_id: str, meta: dict, reason: str) -> None:
    record = make_record(
        source="cBioPortal",
        accession=study_id,
        url=f"https://www.cbioportal.org/study/summary?id={study_id}",
        sample_count=0,
        modality="genomics_mutations_cna",
        license_="Public (cBioPortal open access)",
        raw_files=[],
        checksums={},
        notes=f"ACQUISITION GAP: {reason}. {meta.get('notes','')}",
        is_processed=False,
        atc_specific=meta.get("atc_specific", False),
    )
    prov_write(study_id, record)


def run() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    for study_id, meta in STUDIES.items():
        run_study(study_id, meta)
        time.sleep(1.0)
