"""
Phase 3 — Feature: Evidence Independence
Sources:
  - data/processed/atc_mutations.parquet     (genomics cohort evidence)
  - data/processed/atc_cna.parquet           (genomics cohort evidence)
  - data/processed/depmap_atc_dependency.csv (functional evidence)
  - data/processed/GSE193581_annotated.h5ad  (scRNA-seq evidence)
  - data/processed/GSE76039_expression.parquet (bulk RNA evidence)
  - literature/atc_pubmed_corpus.json        (literature evidence)

Purpose:
  Measures how many INDEPENDENT lines of evidence support each gene as
  an ATC candidate. This is a key anti-leakage metric: genes with broad,
  cross-modality, cross-lab support are more credible than genes prominent
  in a single study.

Outputs (indexed by hgnc_symbol):
    ev_n_modalities         — number of distinct evidence modalities (0-5):
                              genomics, transcriptomics_bulk, transcriptomics_sc,
                              functional_dependency, literature
    ev_genomics_flag        — 1 if any protein-altering mutation or focal CNA
    ev_bulk_rna_flag        — 1 if |bulk_atc_vs_pdtc_fc| > 1.0 (2-fold change)
    ev_scrna_flag           — 1 if scrna_specificity > 0.5 in malignant cells
    ev_depmap_flag          — 1 if median Chronos score < -0.3 across ATC lines
    ev_literature_flag      — 1 if ≥2 ATC-specific PubMed mentions
    ev_cross_modality       — binary: gene supported by ≥2 distinct modalities
    ev_lit_n_atc_mentions   — raw ATC paper count (int, for reference)

Notes:
  - All thresholds documented in feature_spec.md.
  - This module intentionally does NOT import feature values from other modules
    to avoid circular dependencies. It re-derives flags from the processed data
    files directly (independent code path, auditable).
  - The feature table join in phase3_run.py will merge this with other features.
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

ROOT = Path(__file__).parents[2]
DATA = ROOT / "data" / "processed"
LIT  = ROOT / "literature"

# Thresholds (documented)
BULK_FC_THRESHOLD = 1.0         # |log2FC| > 1.0 (ATC vs PDTC)
SCRNA_SPEC_THRESHOLD = 0.5      # malignant specificity score
DEPMAP_THRESHOLD = -0.3         # Chronos score (moderate dependency)
LIT_MIN_MENTIONS = 2            # minimum ATC PubMed mentions


def _genomics_flag(symbols: list[str]) -> pd.Series:
    """1 if gene has ≥1 somatic mutation OR focal CNA in any ATC line."""
    mut = pd.read_parquet(DATA / "atc_mutations.parquet")
    cna = pd.read_parquet(DATA / "atc_cna.parquet")

    mut_genes = set(mut["hgnc_symbol"].dropna())
    # Focal CNA: amp or del
    focal_amp = set(cna.loc[cna["log2_copy_ratio"] > 1.0, "hgnc_symbol"].dropna())
    focal_del = set(cna.loc[cna["log2_copy_ratio"] < -1.0, "hgnc_symbol"].dropna())
    genomic_genes = mut_genes | focal_amp | focal_del

    return pd.Series(
        [int(s in genomic_genes) for s in symbols],
        index=pd.Index(symbols, name="hgnc_symbol"),
        name="ev_genomics_flag",
    )


def _bulk_rna_flag(symbols: list[str]) -> pd.Series:
    """1 if |ATC - PDTC mean expression| > threshold (log-space fold-change)."""
    expr = pd.read_parquet(DATA / "GSE76039_expression.parquet")
    meta = pd.read_parquet(DATA / "GSE76039_metadata.parquet")

    atc_cols  = meta.loc[meta["tumor_type"] == "ATC",  "sample_id"].values
    pdtc_cols = meta.loc[meta["tumor_type"] == "PDTC", "sample_id"].values

    atc_cols  = [c for c in atc_cols  if c in expr.columns]
    pdtc_cols = [c for c in pdtc_cols if c in expr.columns]

    fc = (expr[atc_cols].mean(axis=1) - expr[pdtc_cols].mean(axis=1)).abs()
    flag = (fc > BULK_FC_THRESHOLD).astype(int)
    flag.index.name = "hgnc_symbol"

    result = pd.Series(
        pd.DataFrame(index=pd.Index(symbols, name="hgnc_symbol"))
        .join(flag.rename("ev_bulk_rna_flag"), how="left")["ev_bulk_rna_flag"]
        .fillna(0)
        .astype(int)
    )
    result.index = pd.Index(symbols, name="hgnc_symbol")
    result.name = "ev_bulk_rna_flag"
    return result


def _scrna_flag(symbols: list[str]) -> pd.Series:
    """1 if malignant specificity score > threshold."""
    try:
        import anndata as ad
        h = ad.read_h5ad(DATA / "GSE193581_annotated.h5ad")
        gene_mask = h.var_names.isin(set(symbols))
        h = h[:, gene_mask]

        import scipy.sparse as sp
        mal_mask = (h.obs["malignant_status"] == "malignant").values
        X = h.X
        if sp.issparse(X):
            mal_mean  = np.asarray(X[mal_mask].mean(axis=0)).flatten()
            nonm_mean = np.asarray(X[~mal_mask].mean(axis=0)).flatten()
        else:
            mal_mean  = X[mal_mask].mean(axis=0)
            nonm_mean = X[~mal_mask].mean(axis=0)

        spec = np.log2(mal_mean + 1) - np.log2(nonm_mean + 1)
        spec_series = pd.Series(spec, index=h.var_names, name="spec")

        result = (
            pd.Series(index=pd.Index(symbols, name="hgnc_symbol"), dtype=float, name="ev_scrna_flag")
            .combine_first(spec_series.rename("ev_scrna_flag"))
        )
        result = (result > SCRNA_SPEC_THRESHOLD).fillna(False).astype(int)
        result = result.reindex(symbols, fill_value=0)
        result.index.name = "hgnc_symbol"
        return result
    except Exception as e:
        log.warning("scRNA flag computation failed: %s — defaulting to 0", e)
        return pd.Series(0, index=pd.Index(symbols, name="hgnc_symbol"), name="ev_scrna_flag")


def _depmap_flag(symbols: list[str]) -> pd.Series:
    """1 if median Chronos score across ATC lines < threshold."""
    dep = pd.read_csv(DATA / "depmap_atc_dependency.csv")
    dep.columns = dep.columns.str.strip()
    median_chronos = dep.groupby("hgnc_symbol")["chronos_score"].median()

    flag = (median_chronos < DEPMAP_THRESHOLD).astype(int)
    flag.index.name = "hgnc_symbol"

    result = (
        pd.DataFrame(index=pd.Index(symbols, name="hgnc_symbol"))
        .join(flag.rename("ev_depmap_flag"), how="left")["ev_depmap_flag"]
        .fillna(0)
        .astype(int)
    )
    result.index = pd.Index(symbols, name="hgnc_symbol")
    result.name = "ev_depmap_flag"
    return result


def _literature_flag(symbols: list[str]) -> tuple[pd.Series, pd.Series]:
    """Return (flag, n_mentions): 1 if >= LIT_MIN_MENTIONS ATC papers cite the gene."""
    with open(LIT / "atc_pubmed_corpus.json", encoding="utf-8") as f:
        records = json.load(f)

    TOKEN_RE = re.compile(r'\b[A-Z][A-Z0-9]{2,14}\b')
    symbol_set = {s for s in symbols if len(s) >= 3 and s == s.upper()}
    counts = {s: 0 for s in symbols}

    for rec in records:
        text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
        hits = set(TOKEN_RE.findall(text)) & symbol_set
        for sym in hits:
            counts[sym] += 1

    n_series = pd.Series(counts, name="ev_lit_n_atc_mentions")
    n_series.index.name = "hgnc_symbol"
    flag_series = (n_series >= LIT_MIN_MENTIONS).astype(int).rename("ev_literature_flag")
    return flag_series, n_series


def compute(gene_universe: pd.DataFrame) -> pd.DataFrame:
    """Return feature DataFrame indexed by hgnc_symbol."""
    symbols = gene_universe["hgnc_symbol"].tolist()

    log.info("Evidence independence: computing per-modality flags...")
    g_flag  = _genomics_flag(symbols)
    b_flag  = _bulk_rna_flag(symbols)
    s_flag  = _scrna_flag(symbols)
    d_flag  = _depmap_flag(symbols)
    l_flag, l_count = _literature_flag(symbols)

    df = pd.DataFrame({
        "ev_genomics_flag":     g_flag,
        "ev_bulk_rna_flag":     b_flag,
        "ev_scrna_flag":        s_flag,
        "ev_depmap_flag":       d_flag,
        "ev_literature_flag":   l_flag,
        "ev_lit_n_atc_mentions": l_count,
    }, index=pd.Index(symbols, name="hgnc_symbol"))

    flag_cols = ["ev_genomics_flag", "ev_bulk_rna_flag", "ev_scrna_flag",
                 "ev_depmap_flag", "ev_literature_flag"]
    df["ev_n_modalities"] = df[flag_cols].sum(axis=1).astype(int)
    df["ev_cross_modality"] = (df["ev_n_modalities"] >= 2).astype(int)

    log.info(
        "Evidence independence done. Cross-modality support (≥2 modalities): %d genes",
        df["ev_cross_modality"].sum()
    )
    return df
