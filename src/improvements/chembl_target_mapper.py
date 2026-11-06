"""
Improvement 1 — ChEMBL Target→Gene Mapper
==========================================
Fetches target component (gene) data from ChEMBL REST API for all
target_chembl_ids present in drug_mechanisms.json.

Produces: data/raw/chembl/target_gene_map.json
  {target_chembl_id: {gene_symbols: [...], target_name: str, target_type: str}}

Then rewrites protein_structural.py ChEMBL features with proper coverage.
"""
from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import requests

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

ROOT = Path(__file__).resolve().parents[2]
MECH_PATH = ROOT / "data" / "raw" / "chembl" / "drug_mechanisms.json"
OUT_PATH  = ROOT / "data" / "raw" / "chembl" / "target_gene_map.json"

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data/target"
BATCH_SIZE = 50
SLEEP_BETWEEN = 0.3  # seconds between requests — respect rate limit


def fetch_target_gene_map() -> dict:
    """Fetch gene symbols for all unique target_chembl_ids from ChEMBL REST API."""
    import concurrent.futures

    with open(MECH_PATH) as f:
        mechs = json.load(f)

    target_ids = list({m.get("target_chembl_id", "") for m in mechs if m.get("target_chembl_id")})
    log.info("Fetching gene data for %d unique ChEMBL target IDs...", len(target_ids))

    def fetch_one(tid: str) -> tuple[str, dict]:
        url = f"{CHEMBL_API}/{tid}.json"
        try:
            resp = requests.get(url, timeout=15)
            resp.raise_for_status()
            tgt = resp.json()
            gene_syms = []
            for comp in tgt.get("target_components", []):
                for syn in comp.get("target_component_synonyms", []):
                    if syn.get("syn_type") in ("GENE_SYMBOL", "HGNC_SYMBOL"):
                        sym = syn.get("component_synonym", "").strip().upper()
                        if sym:
                            gene_syms.append(sym)
            return tid, {
                "gene_symbols": list(set(gene_syms)),
                "target_name": tgt.get("pref_name", ""),
                "target_type": tgt.get("target_type", ""),
            }
        except Exception:
            return tid, {"gene_symbols": [], "target_name": "", "target_type": ""}

    result = {}
    # 10 parallel workers — stays within rate limits
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as pool:
        futures = {pool.submit(fetch_one, tid): tid for tid in target_ids}
        done = 0
        for future in concurrent.futures.as_completed(futures):
            tid, data = future.result()
            result[tid] = data
            done += 1
            if done % 200 == 0:
                with_genes = sum(1 for v in result.values() if v["gene_symbols"])
                log.info("  %d / %d done (%d with genes)", done, len(target_ids), with_genes)

    OUT_PATH.write_text(json.dumps(result, indent=2))
    with_genes = sum(1 for v in result.values() if v["gene_symbols"])
    log.info("Wrote %s — %d targets, %d with gene symbols", OUT_PATH, len(result), with_genes)
    return result

    OUT_PATH.write_text(json.dumps(result, indent=2))
    log.info("Wrote %s with %d target mappings", OUT_PATH, len(result))

    # Report coverage
    with_genes = sum(1 for v in result.values() if v["gene_symbols"])
    log.info("Targets with >=1 gene symbol: %d / %d", with_genes, len(result))
    return result


def build_chembl_features(gene_universe_path: Path) -> dict:
    """
    Build gene→{modulator_count, approved_count, clinical_count} using target map.
    Returns dict suitable for adding to feature table.
    """
    import pandas as pd

    if not OUT_PATH.exists():
        log.info("target_gene_map.json not found — running fetch...")
        fetch_target_gene_map()

    with open(OUT_PATH) as f:
        target_map = json.load(f)
    with open(MECH_PATH) as f:
        mechs = json.load(f)

    gu = pd.read_parquet(gene_universe_path)
    valid_symbols = set(gu["hgnc_symbol"].str.upper())

    modulator_sets: dict[str, set] = {}
    approved_sets: dict[str, set] = {}
    clinical_sets: dict[str, set] = {}

    for mech in mechs:
        tid = mech.get("target_chembl_id", "")
        mol = mech.get("molecule_chembl_id", "")
        phase = int(mech.get("max_phase") or 0)

        syms = target_map.get(tid, {}).get("gene_symbols", [])
        for sym in syms:
            if sym in valid_symbols:
                modulator_sets.setdefault(sym, set()).add(mol)
                if phase >= 1:
                    clinical_sets.setdefault(sym, set()).add(mol)
                if phase >= 4:
                    approved_sets.setdefault(sym, set()).add(mol)

    # Normalize back to original case
    sym_case_map = {s.upper(): s for s in gu["hgnc_symbol"]}
    result = {}
    for sym_upper in modulator_sets:
        orig = sym_case_map.get(sym_upper, sym_upper)
        result[orig] = {
            "chembl_modulator_count": len(modulator_sets[sym_upper]),
            "chembl_clinical_count": len(clinical_sets.get(sym_upper, set())),
            "chembl_approved_count": len(approved_sets.get(sym_upper, set())),
        }

    log.info("ChEMBL feature coverage: %d genes with >=1 modulator", len(result))
    return result


if __name__ == "__main__":
    fetch_target_gene_map()
