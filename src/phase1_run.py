"""
phase1_run.py
─────────────
Phase 1 orchestrator. Runs all data acquisition fetchers in priority order:
  1. Gene Universe (prerequisite for all joins)
  2. DepMap (high priority)
  3. GEO transcriptomics
  4. cBioPortal genomics
  5. Open Targets
  6. PubMed literature (slow — large corpora)
  7. ClinicalTrials.gov
  8. STRING PPI
  9. UniProt
  10. Reactome
  11. ChEMBL
  12. ENCODE (expected: no ATC data)
  13. PRIDE (expected: sparse)
  14. Generate dataset_inventory.md

Run: python -m src.phase1_run
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

log = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]

def _step(name: str, fn) -> bool:
    log.info("════ %s ════", name)
    try:
        fn()
        log.info("✓ %s complete.", name)
        return True
    except Exception as exc:
        log.error("✗ %s FAILED: %s", name, exc)
        _log_decision(f"Phase 1 step '{name}' failed: {exc}")
        return False


def _log_decision(msg: str) -> None:
    log_path = ROOT / "docs" / "decisions_log.md"
    entry = (
        f"\n## {datetime.now(timezone.utc).strftime('%Y-%m-%d')} [Phase 1] — Runtime Gap\n"
        f"**Decision:** {msg}\n"
        f"**Logged by:** phase1_run.py\n\n---\n"
    )
    with open(log_path, "a") as f:
        f.write(entry)


def generate_inventory() -> None:
    """Build docs/dataset_inventory.md from all provenance JSON files."""
    prov_dir = ROOT / "data" / "provenance"
    out_path = ROOT / "docs" / "dataset_inventory.md"

    records = []
    for p in sorted(prov_dir.glob("*.json")):
        try:
            r = json.loads(p.read_text())
            records.append(r)
        except Exception:
            pass

    lines = [
        "# QRETIX-ATC v0.1 — Dataset Inventory",
        "",
        f"_Generated: {datetime.now(timezone.utc).isoformat()}_",
        "",
        "| Accession | Source | Modality | Samples | Date Pulled | License | ATC-Specific | Notes |",
        "|-----------|--------|----------|---------|-------------|---------|--------------|-------|",
    ]
    for r in records:
        lines.append(
            f"| {r.get('accession','')} "
            f"| {r.get('source','')} "
            f"| {r.get('modality','')} "
            f"| {r.get('sample_count','')} "
            f"| {str(r.get('date_pulled',''))[:10]} "
            f"| {r.get('license','')} "
            f"| {r.get('atc_specific', '—')} "
            f"| {str(r.get('notes',''))[:120]} |"
        )

    lines += [
        "",
        "## Data Gaps (Manual Approval Required)",
        "",
        "| Source | Reason | Action Required |",
        "|--------|--------|-----------------|",
        "| COSMIC | License acceptance required | https://cancer.sanger.ac.uk/cosmic |",
        "| DrugBank | Licensed database | Academic/commercial license required |",
        "| DisGeNET | Registration required | https://www.disgenet.org/api/ |",
        "| BioGRID | Free API key required | https://webservice.thebiogrid.org/ |",
    ]

    out_path.write_text("\n".join(lines), encoding="utf-8")
    log.info("dataset_inventory.md written: %d datasets → %s", len(records), out_path)


def run() -> None:
    from src import gene_universe
    from src.data_acquisition import (
        depmap_fetcher, geo_fetcher, cbioportal_fetcher,
        opentargets_fetcher, pubmed_fetcher, clinicaltrials_fetcher,
        string_fetcher, uniprot_fetcher, reactome_fetcher,
        chembl_fetcher, encode_fetcher, pride_fetcher,
    )

    steps = [
        ("1.0 Gene Universe",        gene_universe.build),
        ("1.4 DepMap",               depmap_fetcher.run),
        ("1.1 GEO Transcriptomics",  geo_fetcher.run),
        ("1.2 cBioPortal Genomics",  cbioportal_fetcher.run),
        ("1.5 Open Targets",         opentargets_fetcher.run),
        ("1.6 PubMed Literature",    pubmed_fetcher.run),
        ("1.6 ClinicalTrials.gov",   clinicaltrials_fetcher.run),
        ("1.5 STRING PPI",           string_fetcher.run),
        ("1.5 UniProt",              uniprot_fetcher.run),
        ("1.5 Reactome",             reactome_fetcher.run),
        ("1.5 ChEMBL",               chembl_fetcher.run),
        ("1.3 ENCODE",               encode_fetcher.run),
        ("1.3 PRIDE",                pride_fetcher.run),
        ("Inventory",                generate_inventory),
    ]

    results = {}
    for name, fn in steps:
        results[name] = _step(name, fn)

    log.info("\n════ Phase 1 Summary ════")
    for name, ok in results.items():
        log.info("  %s  %s", "✓" if ok else "✗", name)

    failed = [n for n, ok in results.items() if not ok]
    if failed:
        log.warning("Phase 1 completed with %d failures: %s", len(failed), failed)
        log.warning("Check docs/decisions_log.md for details.")
    else:
        log.info("Phase 1 complete — all steps succeeded.")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )
    run()
