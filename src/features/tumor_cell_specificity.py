"""
Phase 3 — Feature: Tumor Cell Specificity
Inputs : data/processed/GSE193581_annotated.h5ad  (scRNA)
         data/processed/GSE76039_expression.parquet + GSE76039_metadata.parquet (bulk)
Outputs: DataFrame indexed by hgnc_symbol with columns:
    scrna_malignant_mean   — mean log-norm expression in malignant cells
    scrna_nonmal_mean      — mean log-norm expression in non-malignant cells
    scrna_specificity      — log2(malignant_mean + 1) - log2(nonmal_mean + 1)
    bulk_atc_mean          — mean expression in ATC samples (GSE76039, already log-scaled)
    bulk_pdtc_mean         — mean expression in PDTC samples
    bulk_atc_vs_pdtc_fc    — bulk_atc_mean - bulk_pdtc_mean (log-space fold-change)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

ROOT = Path(__file__).parents[2]
DATA = ROOT / "data" / "processed"


def compute(gene_universe: pd.DataFrame) -> pd.DataFrame:
    symbols = set(gene_universe["hgnc_symbol"])

    # --- scRNA ---
    h = ad.read_h5ad(DATA / "GSE193581_annotated.h5ad")
    # keep only genes in universe
    gene_mask = h.var_names.isin(symbols)
    h = h[:, gene_mask]

    mal_mask = (h.obs["malignant_status"] == "malignant").values
    X = h.X  # sparse or dense, log-normalised by scanpy in scrna_qc

    import scipy.sparse as sp
    if sp.issparse(X):
        mal_mean = np.asarray(X[mal_mask].mean(axis=0)).flatten()
        nonmal_mean = np.asarray(X[~mal_mask].mean(axis=0)).flatten()
    else:
        mal_mean = X[mal_mask].mean(axis=0)
        nonmal_mean = X[~mal_mask].mean(axis=0)

    scrna = pd.DataFrame({
        "scrna_malignant_mean": mal_mean,
        "scrna_nonmal_mean": nonmal_mean,
    }, index=h.var_names)
    scrna["scrna_specificity"] = (
        np.log2(scrna["scrna_malignant_mean"] + 1) -
        np.log2(scrna["scrna_nonmal_mean"] + 1)
    )

    # --- bulk ---
    expr = pd.read_parquet(DATA / "GSE76039_expression.parquet")   # genes × samples
    meta = pd.read_parquet(DATA / "GSE76039_metadata.parquet")

    atc_samples = meta.loc[meta["tumor_type"] == "ATC", "sample_id"].values
    pdtc_samples = meta.loc[meta["tumor_type"] == "PDTC", "sample_id"].values

    atc_cols = [c for c in atc_samples if c in expr.columns]
    pdtc_cols = [c for c in pdtc_samples if c in expr.columns]

    bulk = pd.DataFrame({
        "bulk_atc_mean": expr[atc_cols].mean(axis=1),
        "bulk_pdtc_mean": expr[pdtc_cols].mean(axis=1),
    })
    bulk["bulk_atc_vs_pdtc_fc"] = bulk["bulk_atc_mean"] - bulk["bulk_pdtc_mean"]
    bulk.index.name = "hgnc_symbol"

    # --- join onto gene universe ---
    df = (
        pd.DataFrame(index=pd.Index(sorted(symbols), name="hgnc_symbol"))
        .join(scrna, how="left")
        .join(bulk, how="left")
    )
    df[["scrna_malignant_mean", "scrna_nonmal_mean", "scrna_specificity"]] = \
        df[["scrna_malignant_mean", "scrna_nonmal_mean", "scrna_specificity"]].fillna(0.0)
    return df
