"""
pubmed_fetcher.py
─────────────────
Phase 1.6: PubMed literature corpus acquisition.

Corpora:
  1. ATC-specific : ("anaplastic thyroid carcinoma"[tiab] OR
                     "anaplastic thyroid cancer"[tiab])
     → literature/atc_pubmed_corpus.json

  2. Broad thyroid : ("thyroid cancer"[tiab] OR "thyroid carcinoma"[tiab])
     → literature/thyroid_pubmed_corpus.json

Every record stores: PMID, title, abstract, publication_date, journal, MeSH terms.
publication_date is MANDATORY — required for Phase 7 temporal splits.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any

import requests
from Bio import Entrez

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
LIT_DIR = ROOT / "literature"

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

QUERIES = {
    "atc": (
        '("anaplastic thyroid carcinoma"[tiab] OR "anaplastic thyroid cancer"[tiab])',
        "atc_pubmed_corpus.json",
    ),
    "thyroid": (
        '("thyroid cancer"[tiab] OR "thyroid carcinoma"[tiab])',
        "thyroid_pubmed_corpus.json",
    ),
}

BATCH_SIZE = 200   # efetch XML batch size
SLEEP = 0.11       # 10 req/s with API key


def _load_keys() -> tuple[str, str]:
    cfg = json.loads((ROOT / "configs" / "api_keys.json").read_text())
    return cfg["ncbi_api_key"], cfg["ncbi_email"]


def _search_pmids(query: str, api_key: str, email: str) -> list[str]:
    """Return all PMIDs matching query."""
    params = {
        "db": "pubmed", "term": query, "retmax": 100000,
        "retmode": "json", "usehistory": "y",
        "api_key": api_key, "email": email,
    }
    resp = requests.get(ESEARCH_URL, params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()
    pmids = data["esearchresult"]["idlist"]
    total = data["esearchresult"]["count"]
    log.info("Query '%s': %s total, retrieved %d PMIDs.", query[:60], total, len(pmids))
    return pmids


def _parse_pubdate(article: Any) -> str:
    """Extract publication date from a parsed Medline record."""
    try:
        pub = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]
        year  = str(pub.get("Year",  pub.get("MedlineDate", "")[:4]))
        month = str(pub.get("Month", "01"))
        day   = str(pub.get("Day",   "01"))
        # Normalise month abbreviations
        months = {"Jan":"01","Feb":"02","Mar":"03","Apr":"04","May":"05","Jun":"06",
                  "Jul":"07","Aug":"08","Sep":"09","Oct":"10","Nov":"11","Dec":"12"}
        month = months.get(month, month.zfill(2))
        day   = day.zfill(2)
        return f"{year}-{month}-{day}"
    except Exception:
        return ""


def _parse_mesh(article: Any) -> list[str]:
    try:
        mesh_list = article["MedlineCitation"].get("MeshHeadingList", [])
        return [str(m["DescriptorName"]) for m in mesh_list]
    except Exception:
        return []


def fetch_corpus(query: str, out_path: Path, api_key: str, email: str) -> None:
    """Fetch all abstracts for a query and write to out_path as JSON."""
    if out_path.exists():
        existing = json.loads(out_path.read_text())
        log.info("%s already exists (%d records), skipping.", out_path.name, len(existing))
        return

    Entrez.api_key = api_key
    Entrez.email   = email

    pmids = _search_pmids(query, api_key, email)
    if not pmids:
        log.warning("No PMIDs found for query: %s", query)
        out_path.write_text(json.dumps([]))
        return

    records: list[dict] = []
    total = len(pmids)

    for i in range(0, total, BATCH_SIZE):
        batch = pmids[i:i + BATCH_SIZE]
        log.info("Fetching abstracts %d–%d / %d …", i + 1, min(i + BATCH_SIZE, total), total)

        try:
            handle = Entrez.efetch(
                db="pubmed", id=",".join(batch),
                rettype="medline", retmode="xml",
            )
            parsed = Entrez.read(handle)
            handle.close()
        except Exception as exc:
            log.warning("Batch %d efetch failed: %s — retrying once.", i, exc)
            time.sleep(2)
            try:
                handle = Entrez.efetch(
                    db="pubmed", id=",".join(batch),
                    rettype="medline", retmode="xml",
                )
                parsed = Entrez.read(handle)
                handle.close()
            except Exception as exc2:
                log.error("Batch %d failed after retry: %s", i, exc2)
                continue

        articles = parsed.get("PubmedArticle", [])
        for art in articles:
            try:
                medline = art["MedlineCitation"]
                article = medline["Article"]
                pmid = str(medline["PMID"])
                title = str(article.get("ArticleTitle", ""))
                abstract_parts = article.get("Abstract", {}).get("AbstractText", [])
                abstract = " ".join(str(p) for p in abstract_parts) if abstract_parts else ""
                journal = str(article["Journal"].get("Title", ""))
                pub_date = _parse_pubdate(art)
                mesh = _parse_mesh(art)

                records.append({
                    "pmid": pmid,
                    "title": title,
                    "abstract": abstract,
                    "journal": journal,
                    "publication_date": pub_date,
                    "mesh_terms": mesh,
                })
            except Exception as exc:
                log.debug("Record parse error: %s", exc)

        # Checkpoint every 1000 records
        if len(records) % 1000 < BATCH_SIZE:
            out_path.write_text(json.dumps(records, indent=2))

        time.sleep(SLEEP)

    out_path.write_text(json.dumps(records, indent=2))
    log.info("Corpus written: %d records → %s", len(records), out_path)

    # Validate: every record must have publication_date
    missing_dates = sum(1 for r in records if not r.get("publication_date"))
    if missing_dates:
        log.warning("%d records missing publication_date (%.1f%%).",
                    missing_dates, 100 * missing_dates / len(records))


def run() -> None:
    LIT_DIR.mkdir(parents=True, exist_ok=True)
    api_key, email = _load_keys()

    for corpus_name, (query, filename) in QUERIES.items():
        out_path = LIT_DIR / filename
        log.info("── PubMed corpus: %s ──", corpus_name)
        fetch_corpus(query, out_path, api_key, email)
        time.sleep(1.0)
