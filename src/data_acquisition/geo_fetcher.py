"""
geo_fetcher.py
──────────────
Phase 1: GEO transcriptomics acquisition.
  - GSE193581  : scRNA-seq (10 ATC, 7 PTC, 6 normal)
  - GSE76039   : bulk microarray PDTC/ATC
  - GSM7030310 : verify if standalone or part of GSE193581
  - GEO search : enumerate all anaplastic thyroid cancer datasets
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import GEOparse
import requests

from src.provenance import make_record, write as prov_write, record_path

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw"
PROV_DIR = ROOT / "data" / "provenance"

NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
GEO_ACC_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"

DATASETS = {
    "GSE193581": {
        "modality": "scRNA-seq",
        "notes": "Lu et al. JCI 2023. 10 ATC, 7 PTC, 6 adjacent normal thyroid + ATC cell lines.",
        "license": "Public domain (GEO submission agreement)",
    },
    "GSE76039": {
        "modality": "bulk_microarray",
        "notes": "Landa et al. JCI 2016. 37 PDTC/ATC genomic/transcriptomic characterization.",
        "license": "Public domain (GEO submission agreement)",
    },
}


def _load_keys() -> tuple[str, str]:
    cfg = json.loads((ROOT / "configs" / "api_keys.json").read_text())
    return cfg["ncbi_api_key"], cfg["ncbi_email"]


def _geo_dir(accession: str) -> Path:
    d = RAW_DIR / accession
    d.mkdir(parents=True, exist_ok=True)
    return d


def fetch_gse(accession: str, meta: dict) -> None:
    """Download a GSE series using GEOparse and write provenance."""
    if record_path(accession).exists():
        log.info("%s provenance already exists, skipping.", accession)
        return

    dest = _geo_dir(accession)
    log.info("Fetching %s → %s", accession, dest)

    try:
        gse = GEOparse.get_GEO(geo=accession, destdir=str(dest), silent=True)
    except Exception as exc:
        log.error("Failed to fetch %s: %s", accession, exc)
        _write_gap_provenance(accession, meta, str(exc))
        return

    gsm_count = len(gse.gsms)
    raw_files = [p.name for p in dest.iterdir() if p.is_file()]

    record = make_record(
        source="NCBI GEO",
        accession=accession,
        url=GEO_ACC_URL.format(acc=accession),
        sample_count=gsm_count,
        modality=meta["modality"],
        license_=meta["license"],
        raw_files=raw_files,
        checksums={},          # checksums computed in Phase 2 after full download
        notes=meta["notes"],
        is_processed=False,
        platform=_get_platform(gse),
    )
    prov_write(accession, record)
    log.info("%s: %d samples, provenance written.", accession, gsm_count)


def _get_platform(gse: GEOparse.GEOTypes.GSE) -> str:
    try:
        platforms = list(gse.gpls.keys())
        return ", ".join(platforms) if platforms else "unknown"
    except Exception:
        return "unknown"


def _write_gap_provenance(accession: str, meta: dict, error: str) -> None:
    record = make_record(
        source="NCBI GEO",
        accession=accession,
        url=GEO_ACC_URL.format(acc=accession),
        sample_count=0,
        modality=meta.get("modality", "unknown"),
        license_=meta.get("license", "unknown"),
        raw_files=[],
        checksums={},
        notes=f"ACQUISITION FAILED: {error}",
        is_processed=False,
    )
    prov_write(accession, record)


def verify_gsm7030310() -> None:
    """Verify whether GSM7030310 is standalone or part of GSE193581."""
    out_path = PROV_DIR / "GSM7030310_verification.json"
    if out_path.exists():
        log.info("GSM7030310 verification already done.")
        return

    api_key, email = _load_keys()
    params = {
        "db": "gds", "id": "GSM7030310", "retmode": "json",
        "api_key": api_key, "email": email,
    }
    try:
        resp = requests.get(NCBI_ESUMMARY, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        result = data.get("result", {})
        uids = result.get("uids", [])
        series_found = []
        for uid in uids:
            entry = result.get(uid, {})
            series_found.append(entry.get("gse", "unknown"))

        finding = {
            "accession": "GSM7030310",
            "query_date": __import__("datetime").datetime.utcnow().isoformat(),
            "is_standalone": len(uids) > 0 and "GSE193581" not in str(series_found),
            "parent_series": series_found,
            "raw_response_uids": uids,
            "conclusion": (
                "GSM7030310 is a sample within GSE193581, not an independent dataset."
                if "GSE193581" in str(series_found)
                else f"GSM7030310 parent series: {series_found}. Verify manually."
            ),
        }
    except Exception as exc:
        finding = {
            "accession": "GSM7030310",
            "query_date": __import__("datetime").datetime.utcnow().isoformat(),
            "error": str(exc),
            "conclusion": "Verification failed — manual check required.",
        }

    out_path.write_text(json.dumps(finding, indent=2))
    log.info("GSM7030310 verification → %s: %s", out_path, finding["conclusion"])


def geo_search_enumerate() -> None:
    """
    Run GEO esearch for anaplastic thyroid cancer, enumerate ALL results,
    log accession + metadata. Does NOT download — inventory only.
    """
    out_path = PROV_DIR / "geo_search_results.json"
    if out_path.exists():
        log.info("GEO search results already enumerated.")
        return

    api_key, email = _load_keys()
    query = "anaplastic thyroid cancer[All Fields]"
    log.info("GEO search: '%s'", query)

    # Step 1: get total count + all IDs
    search_params = {
        "db": "gds", "term": query, "retmax": 500,
        "retmode": "json", "usehistory": "y",
        "api_key": api_key, "email": email,
    }
    try:
        resp = requests.get(NCBI_ESEARCH, params=search_params, timeout=30)
        resp.raise_for_status()
        search_data = resp.json()
    except Exception as exc:
        log.error("GEO search failed: %s", exc)
        out_path.write_text(json.dumps({"error": str(exc)}, indent=2))
        return

    id_list = search_data.get("esearchresult", {}).get("idlist", [])
    total = search_data.get("esearchresult", {}).get("count", "0")
    log.info("GEO search: %s total results, fetching summaries for %d IDs.", total, len(id_list))

    # Step 2: fetch summaries in batches of 50
    results = []
    batch_size = 50
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i + batch_size]
        sum_params = {
            "db": "gds", "id": ",".join(batch),
            "retmode": "json", "api_key": api_key, "email": email,
        }
        try:
            sr = requests.get(NCBI_ESUMMARY, params=sum_params, timeout=30)
            sr.raise_for_status()
            sdata = sr.json().get("result", {})
            for uid in sdata.get("uids", []):
                entry = sdata[uid]
                results.append({
                    "uid": uid,
                    "accession": entry.get("accession", ""),
                    "title": entry.get("title", ""),
                    "gdstype": entry.get("gdstype", ""),
                    "platform": entry.get("platform", ""),
                    "n_samples": entry.get("n_samples", ""),
                    "pubdate": entry.get("pubdate", ""),
                    "summary": entry.get("summary", "")[:300],
                })
        except Exception as exc:
            log.warning("Batch %d summary fetch failed: %s", i, exc)
        time.sleep(0.11)  # 10 req/s with API key

    output = {
        "query": query,
        "query_date": __import__("datetime").datetime.utcnow().isoformat(),
        "total_reported": total,
        "retrieved": len(results),
        "results": results,
    }
    out_path.write_text(json.dumps(output, indent=2))
    log.info("GEO search enumeration complete: %d datasets logged → %s", len(results), out_path)


def run() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PROV_DIR.mkdir(parents=True, exist_ok=True)

    for accession, meta in DATASETS.items():
        fetch_gse(accession, meta)
        time.sleep(0.5)

    verify_gsm7030310()
    geo_search_enumerate()
