"""
Phase 3 — Feature: Protein / Structural Features
Sources:
  - UniProt human proteome TSV  (data/raw/uniprot/uniprot_human_proteome.tsv)
  - ChEMBL drug mechanisms JSON (data/raw/chembl/drug_mechanisms.json)
  - Open Targets tractability   (data/raw/opentargets/tractability.json)
    (tractability tier already computed in disease_relevance.py — imported here
     so protein_structural is self-contained; no cross-module import needed)

Outputs (indexed by hgnc_symbol):
    uniprot_length          — protein sequence length (aa)
    uniprot_domain_count    — number of annotated functional domains
    uniprot_has_structure   — 1 if at least one PDB entry listed in UniProt
    uniprot_is_membrane     — 1 if subcellular location contains 'Membrane' or 'membrane'
    uniprot_is_secreted     — 1 if subcellular location contains 'Secreted'
    uniprot_is_reviewed     — 1 if SwissProt reviewed entry
    chembl_modulator_count  — number of distinct ChEMBL compounds targeting this gene
                              (from drug_mechanisms.json; any max_phase)
    chembl_approved_count   — number with max_phase == 4 (approved)
    chembl_clinical_count   — number with max_phase >= 1 (clinical)

Notes:
  - UniProt is joined via hgnc_symbol (Gene Names primary column).
  - ChEMBL drug_mechanisms.json lacks direct gene symbol — joined via target_chembl_id
    mapped through approved_drugs.json target field (best-effort; some gaps expected).
  - Domain count parsed from 'Domain [FT]' column (count semicolon-separated entries).
  - Missing values: length=NaN (gene not in UniProt), counts=0 (no ChEMBL data).
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

ROOT = Path(__file__).parents[2]
RAW_UNIPROT = ROOT / "data" / "raw" / "uniprot" / "uniprot_human_proteome.tsv"
RAW_CHEMBL_MECH = ROOT / "data" / "raw" / "chembl" / "drug_mechanisms.json"
RAW_CHEMBL_DRUGS = ROOT / "data" / "raw" / "chembl" / "approved_drugs.json"


def _load_uniprot() -> pd.DataFrame:
    """Parse UniProt TSV → DataFrame keyed by hgnc_symbol (primary gene name)."""
    df = pd.read_csv(RAW_UNIPROT, sep="\t", low_memory=False)
    df = df.rename(columns={
        "Gene Names (primary)": "hgnc_symbol",
        "Length": "uniprot_length",
        "Subcellular location [CC]": "subcell",
        "Domain [FT]": "domain_ft",
        "PDB": "pdb",
        "Reviewed": "is_reviewed",
    })

    # Drop rows with no primary gene symbol
    df = df.dropna(subset=["hgnc_symbol"])
    df["hgnc_symbol"] = df["hgnc_symbol"].str.strip()

    # is_reviewed
    df["uniprot_is_reviewed"] = (df["is_reviewed"] == "reviewed").astype(int)

    # has_structure: PDB column is non-null and non-empty
    df["uniprot_has_structure"] = (
        df["pdb"].notna() & (df["pdb"].astype(str).str.strip() != "")
    ).astype(int)

    # Domain count: count entries in Domain [FT] column
    def _count_domains(val):
        if pd.isna(val) or str(val).strip() in ("", "nan"):
            return 0
        # Domains are listed like "DOMAIN 1..50; /note=...;"
        return len(re.findall(r"DOMAIN\s+\d+\.\.\d+", str(val)))

    df["uniprot_domain_count"] = df["domain_ft"].apply(_count_domains)

    # Subcellular location flags
    def _has_loc(val, term: str) -> int:
        if pd.isna(val):
            return 0
        return int(term.lower() in str(val).lower())

    df["uniprot_is_membrane"] = df["subcell"].apply(lambda v: _has_loc(v, "membrane"))
    df["uniprot_is_secreted"]  = df["subcell"].apply(lambda v: _has_loc(v, "secreted"))

    # Keep only relevant columns; deduplicate by taking first reviewed entry per symbol
    keep = ["hgnc_symbol", "uniprot_length", "uniprot_domain_count",
            "uniprot_has_structure", "uniprot_is_membrane", "uniprot_is_secreted",
            "uniprot_is_reviewed"]
    df = df[keep].drop_duplicates(subset="hgnc_symbol", keep="first")
    df = df.set_index("hgnc_symbol")

    log.info("Loaded UniProt: %d unique gene symbols", len(df))
    return df


def _load_chembl_modulators(gene_universe: pd.DataFrame) -> pd.DataFrame:
    """
    Parse drug_mechanisms.json to get modulator counts per gene symbol.
    Gene symbols are extracted from mechanism_of_action text patterns.
    Returns DataFrame indexed by hgnc_symbol.
    """
    with open(RAW_CHEMBL_MECH, encoding="utf-8") as f:
        mechanisms = json.load(f)

    valid_symbols = set(gene_universe["hgnc_symbol"].str.upper())

    # Pattern: "GENESYM inhibitor/activator/antagonist/agonist/modulator/blocker"
    MOA_GENE_RE = re.compile(
        r"^([A-Z][A-Z0-9]{1,14})\s+(?:inhibitor|activator|antagonist|agonist|"
        r"modulator|blocker|opener|stabilizer|binder|ligand|substrate)",
        re.IGNORECASE,
    )

    modulator_sets: dict[str, set[str]] = {}
    approved_sets: dict[str, set[str]] = {}
    clinical_sets: dict[str, set[str]] = {}

    for mech in mechanisms:
        mol_id = mech.get("molecule_chembl_id", "")
        max_phase = mech.get("max_phase") or 0
        try:
            max_phase = int(max_phase)
        except (ValueError, TypeError):
            max_phase = 0

        moa = (mech.get("mechanism_of_action") or "").strip()
        m = MOA_GENE_RE.match(moa)
        if not m:
            continue
        sym = m.group(1).upper()
        if sym not in valid_symbols:
            continue

        modulator_sets.setdefault(sym, set()).add(mol_id)
        if max_phase >= 1:
            clinical_sets.setdefault(sym, set()).add(mol_id)
        if max_phase == 4:
            approved_sets.setdefault(sym, set()).add(mol_id)

    all_syms = sorted(set(modulator_sets))
    # Normalize back to original case from gene universe
    sym_orig_map = {s.upper(): s for s in gene_universe["hgnc_symbol"]}
    rows = []
    for sym in all_syms:
        orig = sym_orig_map.get(sym, sym)
        rows.append({
            "hgnc_symbol": orig,
            "chembl_modulator_count": len(modulator_sets.get(sym, set())),
            "chembl_clinical_count": len(clinical_sets.get(sym, set())),
            "chembl_approved_count": len(approved_sets.get(sym, set())),
        })

    df = pd.DataFrame(rows).set_index("hgnc_symbol") if rows else pd.DataFrame(
        columns=["chembl_modulator_count", "chembl_clinical_count", "chembl_approved_count"]
    )
    log.info("ChEMBL modulators: %d genes with >=1 modulator", len(df))
    return df


def compute(gene_universe: pd.DataFrame) -> pd.DataFrame:
    """Return feature DataFrame indexed by hgnc_symbol."""
    symbols = gene_universe["hgnc_symbol"].tolist()

    uniprot_df = _load_uniprot()
    chembl_df  = _load_chembl_modulators(gene_universe)

    # Normalize chembl index to match gene universe case
    chembl_df.index = chembl_df.index.str.upper()
    sym_upper_map = {s.upper(): s for s in symbols}
    chembl_df.index = [sym_upper_map.get(s, s) for s in chembl_df.index]

    result = (
        pd.DataFrame(index=pd.Index(symbols, name="hgnc_symbol"))
        .join(uniprot_df, how="left")
        .join(chembl_df, how="left")
    )

    # Fill integer columns with 0
    int_cols = ["uniprot_domain_count", "uniprot_has_structure",
                "uniprot_is_membrane", "uniprot_is_secreted",
                "uniprot_is_reviewed",
                "chembl_modulator_count", "chembl_clinical_count", "chembl_approved_count"]
    for col in int_cols:
        if col in result.columns:
            result[col] = result[col].fillna(0).astype(int)

    # uniprot_length: keep as NaN (gene not in UniProt proteome)

    log.info("Protein structural: %d genes with UniProt data, %d with ChEMBL data",
             result["uniprot_length"].notna().sum(),
             (result.get("chembl_modulator_count", pd.Series(0)) > 0).sum())

    return result
