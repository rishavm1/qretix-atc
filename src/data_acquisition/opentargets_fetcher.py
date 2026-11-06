"""
opentargets_fetcher.py
──────────────────────
Phase 1.5: Open Targets Platform acquisition.

Pulls:
  - ATC disease-target associations (MONDO:0006468)
  - Tractability scores per target (small molecule / antibody)

Used as:
  - Cross-check for label construction (Phase 4)
  - Benchmark baseline (Phase 9)
  - Novelty validation (Phase 6 eight-point checklist)

Output:
  data/raw/opentargets/atc_associations.json
  data/raw/opentargets/tractability.json
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
RAW_DIR = ROOT / "data" / "raw" / "opentargets"

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
ATC_DISEASE_ID = "MONDO_0006468"   # Anaplastic thyroid carcinoma

ASSOC_QUERY = """
query ATCAssociations($diseaseId: String!, $page: Int!) {
  disease(efoId: $diseaseId) {
    associatedTargets(page: {index: $page, size: 500}) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
          biotype
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"""

TRACTABILITY_QUERY = """
query Tractability($targetId: String!) {
  target(ensemblId: $targetId) {
    id
    approvedSymbol
    tractability {
      label
      modality
      value
    }
  }
}
"""


def _gql(query: str, variables: dict) -> dict | None:
    try:
        resp = requests.post(
            OT_GRAPHQL,
            json={"query": query, "variables": variables},
            timeout=60,
        )
        resp.raise_for_status()
        return resp.json().get("data")
    except Exception as exc:
        log.error("Open Targets GraphQL failed: %s", exc)
        return None


def fetch_atc_associations() -> list[dict]:
    out_path = RAW_DIR / "atc_associations.json"
    if out_path.exists():
        log.info("atc_associations.json already exists.")
        return json.loads(out_path.read_text())

    log.info("Fetching Open Targets ATC associations (MONDO:0006468) …")
    all_rows: list[dict] = []
    page = 0

    while True:
        data = _gql(ASSOC_QUERY, {"diseaseId": ATC_DISEASE_ID, "page": page})
        if not data:
            break
        rows = data["disease"]["associatedTargets"]["rows"]
        if not rows:
            break
        all_rows.extend(rows)
        total = data["disease"]["associatedTargets"]["count"]
        log.info("Open Targets ATC: %d / %d associations fetched.", len(all_rows), total)
        if len(all_rows) >= total:
            break
        page += 1
        time.sleep(0.2)

    out_path.write_text(json.dumps(all_rows, indent=2))
    log.info("Open Targets ATC associations: %d targets → %s", len(all_rows), out_path)
    return all_rows


def fetch_tractability(ensembl_ids: list[str]) -> dict[str, list]:
    """Fetch tractability for a list of Ensembl gene IDs."""
    out_path = RAW_DIR / "tractability.json"
    if out_path.exists():
        log.info("tractability.json already exists.")
        return json.loads(out_path.read_text())

    log.info("Fetching tractability for %d targets …", len(ensembl_ids))
    results: dict[str, list] = {}

    for i, eid in enumerate(ensembl_ids):
        data = _gql(TRACTABILITY_QUERY, {"targetId": eid})
        if data and data.get("target"):
            t = data["target"]
            results[eid] = {
                "symbol": t.get("approvedSymbol", ""),
                "tractability": t.get("tractability", []),
            }
        if i % 100 == 0:
            log.info("Tractability: %d / %d done.", i, len(ensembl_ids))
        time.sleep(0.05)

    out_path.write_text(json.dumps(results, indent=2))
    log.info("Tractability written for %d targets → %s", len(results), out_path)
    return results


def run() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    associations = fetch_atc_associations()

    # Extract Ensembl IDs for tractability fetch
    ensembl_ids = [r["target"]["id"] for r in associations if r.get("target")]
    if ensembl_ids:
        fetch_tractability(ensembl_ids)

    if not record_path("opentargets_atc").exists():
        record = make_record(
            source="Open Targets Platform",
            accession="opentargets_atc",
            url="https://api.platform.opentargets.org/api/v4/graphql",
            sample_count=len(associations),
            modality="disease_target_associations",
            license_="CC BY 4.0",
            raw_files=["atc_associations.json", "tractability.json"],
            checksums={},
            notes=(
                f"ATC disease ID: MONDO:0006468. {len(associations)} target associations. "
                "Used as cross-check for labels (Phase 4) and benchmark baseline (Phase 9). "
                "NOT used as primary label source — cross-check only."
            ),
            is_processed=False,
            atc_specific=True,
        )
        prov_write("opentargets_atc", record)
