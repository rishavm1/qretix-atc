"""
Phase 3 — Feature: Disease Relevance
Sources:
  - Open Targets ATC association scores  (data/raw/opentargets/atc_associations.json)
  - Open Targets tractability            (data/raw/opentargets/tractability.json)

Status: IMPLEMENTED — real data parsed from local raw files.

Outputs (indexed by hgnc_symbol):
    ot_atc_score         — Open Targets overall ATC association score [0,1]
    ot_literature_score  — OT literature evidence sub-score
    ot_somatic_mut_score — OT somatic mutation evidence sub-score
    ot_sm_approved       — 1 if any approved small-molecule drug exists (tractability)
    ot_ab_approved       — 1 if any approved antibody exists (tractability)
    ot_druggability_tier — A=approved drug, B=clinical, C=tractable, D=novel/undrugged

Notes:
  - Genes not in OT ATC associations get score=0 (genuinely no OT association).
  - Tractability tier follows Open Targets SM/AB pipeline labels.
  - All scores traceable to data/raw/opentargets/*.json (provenance in decisions_log).
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

ROOT = Path(__file__).parents[2]
RAW_OT = ROOT / "data" / "raw" / "opentargets"
PROC_OT = ROOT / "data" / "processed" / "opentargets_atc_targets.parquet"


def _load_atc_associations() -> pd.DataFrame:
    """Parse atc_associations.json → DataFrame keyed by ensembl_gene_id."""
    path = RAW_OT / "atc_associations.json"
    with open(path, encoding="utf-8") as f:
        records = json.load(f)

    rows = []
    for rec in records:
        tgt = rec.get("target", {})
        ensembl = tgt.get("id", "")
        symbol  = tgt.get("approvedSymbol", "")
        overall = float(rec.get("score", 0.0))

        # Parse datatype sub-scores
        lit_score = 0.0
        som_score = 0.0
        for dt in rec.get("datatypeScores", []):
            if dt["id"] == "literature":
                lit_score = float(dt["score"])
            elif dt["id"] == "somatic_mutation":
                som_score = float(dt["score"])

        rows.append({
            "ensembl_gene_id": ensembl,
            "hgnc_symbol": symbol,
            "ot_atc_score": overall,
            "ot_literature_score": lit_score,
            "ot_somatic_mut_score": som_score,
        })

    df = pd.DataFrame(rows)
    log.info("Loaded %d OT ATC associations", len(df))
    return df


def _load_tractability() -> pd.DataFrame:
    """
    Parse tractability.json → DataFrame keyed by ensembl_gene_id.
    Derives:
      ot_sm_approved  — has approved small-molecule drug
      ot_ab_approved  — has approved antibody
      ot_druggability_tier — A/B/C/D
    """
    path = RAW_OT / "tractability.json"
    with open(path, encoding="utf-8") as f:
        data = json.load(f)  # dict: ensembl_id → {symbol, tractability: [...]}

    rows = []
    for ensembl, info in data.items():
        sm_approved = False
        ab_approved = False
        sm_clinical = False
        sm_tractable = False

        for item in info.get("tractability", []):
            label    = item.get("label", "")
            modality = item.get("modality", "")
            value    = bool(item.get("value", False))

            if modality == "SM":
                if label == "Approved Drug" and value:
                    sm_approved = True
                if label in ("Advanced Clinical", "Phase 1 Clinical") and value:
                    sm_clinical = True
                if label in ("Structure with Ligand", "High-Quality Ligand",
                             "High-Quality Pocket", "Druggable Family") and value:
                    sm_tractable = True
            elif modality == "AB":
                if label == "Approved Drug" and value:
                    ab_approved = True

        # Tier assignment (best tier wins)
        if sm_approved or ab_approved:
            tier = "A"
        elif sm_clinical:
            tier = "B"
        elif sm_tractable:
            tier = "C"
        else:
            tier = "D"

        rows.append({
            "ensembl_gene_id": ensembl,
            "ot_sm_approved": int(sm_approved),
            "ot_ab_approved": int(ab_approved),
            "ot_druggability_tier": tier,
        })

    df = pd.DataFrame(rows).set_index("ensembl_gene_id")
    log.info("Loaded tractability for %d genes", len(df))
    return df


def compute(gene_universe: pd.DataFrame) -> pd.DataFrame:
    """Return feature DataFrame indexed by hgnc_symbol."""
    sym_to_ensembl = (
        gene_universe.set_index("hgnc_symbol")["ensembl_gene_id"]
        .dropna()
        .to_dict()
    )

    # --- Load OT associations ---
    # Prefer the richer processed parquet (has rna_expression + pathway scores)
    if PROC_OT.exists():
        ot_proc = pd.read_parquet(PROC_OT)
        # Map ensembl → symbol for join
        ensembl_to_sym = {v: k for k, v in sym_to_ensembl.items()}
        ot_proc["hgnc_symbol"] = ot_proc["target_ensembl_id"].map(ensembl_to_sym).fillna(
            ot_proc.get("approved_symbol", pd.Series(dtype=str))
        )
        ot_proc = ot_proc.dropna(subset=["hgnc_symbol"]).set_index("hgnc_symbol")
        # Map processed column names to expected output names
        ot_assoc = pd.DataFrame(index=ot_proc.index)
        ot_assoc["ot_atc_score"] = ot_proc.get("overall_association_score", 0.0)
        ot_assoc["ot_literature_score"] = ot_proc.get("score_literature", 0.0)
        ot_assoc["ot_somatic_mut_score"] = ot_proc.get("score_somatic_mutation", 0.0)
        log.info("Loaded OT from processed parquet: %d targets", len(ot_assoc))
    else:
        ot_assoc = _load_atc_associations().set_index("hgnc_symbol")

    # --- Load tractability (still from raw JSON — richer tier info) ---
    tract_df = _load_tractability()
    ensembl_to_sym = {v: k for k, v in sym_to_ensembl.items()}
    tract_by_sym = tract_df.copy()
    tract_by_sym.index = [ensembl_to_sym.get(e, e) for e in tract_by_sym.index]

    # --- Join onto gene universe ---
    symbols = gene_universe["hgnc_symbol"].tolist()
    base = pd.DataFrame(index=pd.Index(symbols, name="hgnc_symbol"))

    df = (
        base
        .join(ot_assoc[["ot_atc_score", "ot_literature_score", "ot_somatic_mut_score"]], how="left")
        .join(tract_by_sym[["ot_sm_approved", "ot_ab_approved", "ot_druggability_tier"]], how="left")
    )

    for col in ["ot_atc_score", "ot_literature_score", "ot_somatic_mut_score"]:
        df[col] = df[col].fillna(0.0)
    for col in ["ot_sm_approved", "ot_ab_approved"]:
        df[col] = df[col].fillna(0).astype(int)
    df["ot_druggability_tier"] = df["ot_druggability_tier"].fillna("D")

    n_assoc = (df["ot_atc_score"] > 0).sum()
    n_tier_a = (df["ot_druggability_tier"] == "A").sum()
    log.info("Disease relevance: %d genes with OT ATC score>0, %d Tier-A druggable", n_assoc, n_tier_a)
    return df
