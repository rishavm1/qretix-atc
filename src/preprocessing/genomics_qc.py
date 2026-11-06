"""
genomics_qc.py
──────────────
Phase 2.3: Genomics QC.

cBioPortal GATCI 2024 and MSKCC 2016 returned 404 — gap logged.
Using DepMap OmicsSomaticMutations + OmicsCNGene for ATC cell lines as proxy.

Outputs:
  data/processed/atc_mutations.parquet
  data/processed/atc_cna.parquet
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from src import gene_universe

log = logging.getLogger(__name__)

ROOT       = Path(__file__).resolve().parents[2]
PROC_DIR   = ROOT / "data" / "processed"
RAW_DEPMAP = ROOT / "data" / "raw" / "depmap"

MUT_OUT = PROC_DIR / "atc_mutations.parquet"
CNA_OUT = PROC_DIR / "atc_cna.parquet"

PROTEIN_ALTERING = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
    "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site",
}


def _gu() -> pd.DataFrame:
    return gene_universe.load()[["hgnc_symbol", "ensembl_gene_id"]].dropna()


def _atc_models() -> list[str]:
    dep = PROC_DIR / "depmap_atc_dependency.csv"
    if not dep.exists():
        return []
    return pd.read_csv(dep)["ModelID"].unique().tolist()


def process_mutations() -> pd.DataFrame:
    path = RAW_DEPMAP / "OmicsSomaticMutations.csv"
    if not path.exists():
        log.error("OmicsSomaticMutations.csv missing.")
        return pd.DataFrame()

    models = _atc_models()
    if not models:
        log.error("No ATC ModelIDs available.")
        return pd.DataFrame()

    log.info("Loading OmicsSomaticMutations.csv …")
    mut = pd.read_csv(path, low_memory=False)

    model_col = next((c for c in mut.columns
                      if c.lower() in {"modelid", "depmap_id", "model_id"}), None)
    hugo_col  = next((c for c in mut.columns
                      if "hugo" in c.lower() or c.lower() == "genename"), None)
    class_col = next((c for c in mut.columns
                      if "variantclassification" in c.lower().replace("_", "")), None)

    if not model_col or not hugo_col:
        log.error("Cannot identify key columns. Columns: %s", list(mut.columns[:15]))
        return pd.DataFrame()

    atc = mut[mut[model_col].isin(models)].copy()
    log.info("ATC mutations raw: %d", len(atc))

    atc = atc.rename(columns={model_col: "ModelID", hugo_col: "hgnc_symbol"})
    if class_col:
        atc = atc.rename(columns={class_col: "variant_classification"})
        atc = atc[atc["variant_classification"].isin(PROTEIN_ALTERING)]
        log.info("After protein-altering filter: %d", len(atc))

    atc = atc.merge(_gu(), on="hgnc_symbol", how="left")
    atc["platform"] = "DepMap_precalled_somatic"
    atc["batch"]    = "depmap_24Q4"
    atc["source_note"] = (
        "DepMap ATC cell line somatic mutations — proxy for patient tumor mutations. "
        "cBioPortal ATC cohorts unavailable (404). See decisions_log.md."
    )

    keep = ["ensembl_gene_id", "hgnc_symbol", "ModelID",
            "platform", "batch", "source_note"]
    if "variant_classification" in atc.columns:
        keep.append("variant_classification")

    return atc[[c for c in keep if c in atc.columns]]


def process_cna() -> pd.DataFrame:
    path = RAW_DEPMAP / "OmicsCNGene.csv"
    if not path.exists():
        log.error("OmicsCNGene.csv missing.")
        return pd.DataFrame()

    models = _atc_models()
    if not models:
        return pd.DataFrame()

    log.info("Loading OmicsCNGene.csv …")
    cna = pd.read_csv(path, index_col=0, low_memory=False)

    available = [m for m in models if m in cna.index]
    missing   = [m for m in models if m not in cna.index]
    if missing:
        log.warning("ATC models not in CNA: %s", missing)
    if not available:
        log.error("No ATC models in CNA rows.")
        return pd.DataFrame()

    subset = cna.loc[available].reset_index()
    subset = subset.rename(columns={"index": "ModelID",
                                    subset.columns[0]: "ModelID"})
    # ensure first col is ModelID
    subset.columns = ["ModelID"] + list(subset.columns[1:])

    long = subset.melt(id_vars=["ModelID"],
                       var_name="gene_raw",
                       value_name="log2_copy_ratio")
    long["hgnc_symbol"] = long["gene_raw"].str.extract(r"^([^\s(]+)")
    long = long.drop(columns=["gene_raw"])
    long = long.merge(_gu(), on="hgnc_symbol", how="left")
    long = long.dropna(subset=["log2_copy_ratio"])
    long["platform"] = "DepMap_precalled_CNA"
    long["batch"]    = "depmap_24Q4"

    log.info("ATC CNA: %d rows (%d genes × %d lines)",
             len(long), long["hgnc_symbol"].nunique(), long["ModelID"].nunique())
    return long[["ensembl_gene_id", "hgnc_symbol", "ModelID",
                 "log2_copy_ratio", "platform", "batch"]]


def run() -> None:
    PROC_DIR.mkdir(parents=True, exist_ok=True)

    if not MUT_OUT.exists():
        df = process_mutations()
        if df.empty:
            df = pd.DataFrame(columns=["ensembl_gene_id", "hgnc_symbol",
                                        "ModelID", "platform", "batch", "source_note"])
        df.to_parquet(MUT_OUT, index=False)
        log.info("atc_mutations.parquet → %d rows", len(df))
    else:
        log.info("atc_mutations.parquet already exists.")

    if not CNA_OUT.exists():
        df = process_cna()
        if df.empty:
            df = pd.DataFrame(columns=["ensembl_gene_id", "hgnc_symbol",
                                        "ModelID", "log2_copy_ratio",
                                        "platform", "batch"])
        df.to_parquet(CNA_OUT, index=False)
        log.info("atc_cna.parquet → %d rows", len(df))
    else:
        log.info("atc_cna.parquet already exists.")
