"""
Improvement 2 — Biology-informed NaN Imputation
=================================================
Problem: 57% of genes have NaN for DepMap/bulk RNA features because
those datasets only cover ~19k protein-coding genes. Filling with
column median unfairly penalizes genes not in the dataset vs genes
that genuinely have low scores.

Strategy:
  - DepMap NaN → fill with 0.0 (no dependency signal, not missing)
  - Bulk RNA NaN → fill with 0.0 (no differential expression, not missing)
  - UniProt NaN → fill with median of protein-coding genes only
  - ppi_driver_proximity NaN → fill with max_observed + 1 (far from drivers)

Also adds a missingness indicator flag per dataset group so the model
can learn that "not in DepMap" is itself a biological signal.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
FEATURES_PATH = ROOT / "features" / "feature_table.parquet"


def improve_imputation() -> pd.DataFrame:
    ft = pd.read_parquet(FEATURES_PATH)
    log.info("Feature table before imputation: %s", ft.shape)

    # --- DepMap features: NaN = not in DepMap = no dependency (fill 0) ---
    dep_cols = ["dep_mean_chronos", "dep_min_chronos", "dep_essential_frac"]
    dep_n_col = "dep_n_lines"
    for col in dep_cols:
        if col in ft.columns:
            ft[col] = ft[col].fillna(0.0)
    if dep_n_col in ft.columns:
        ft[dep_n_col] = ft[dep_n_col].fillna(0).astype(int)

    # Add missingness flag: was this gene in DepMap at all?
    ft["dep_in_depmap_flag"] = (ft[dep_n_col] > 0).astype(int)

    # --- Bulk RNA: NaN = gene not in microarray probe set = fill 0 ---
    bulk_cols = ["bulk_atc_mean", "bulk_pdtc_mean", "bulk_atc_vs_pdtc_fc"]
    for col in bulk_cols:
        if col in ft.columns:
            ft[col] = ft[col].fillna(0.0)

    # Add missingness flag
    ft["bulk_in_array_flag"] = (ft["bulk_atc_mean"] != 0.0).astype(int)

    # --- UniProt: NaN = not in reviewed proteome ---
    # Fill length with 0 (unknown protein length), flags already 0
    if "uniprot_length" in ft.columns:
        ft["uniprot_in_swissprot_flag"] = ft["uniprot_length"].notna().astype(int)
        ft["uniprot_length"] = ft["uniprot_length"].fillna(0.0)

    # --- ppi_driver_proximity: NaN = isolated from drivers ---
    # Fill with max observed distance + 2 (penalize isolation)
    if "ppi_driver_proximity" in ft.columns:
        max_prox = ft["ppi_driver_proximity"].dropna().max()
        ft["ppi_driver_proximity"] = ft["ppi_driver_proximity"].fillna(max_prox + 2)
        # Add flag: was gene reachable from drivers?
        ft["ppi_driver_reachable_flag"] = (ft["ppi_driver_proximity"] <= max_prox).astype(int)

    # --- Remaining NaNs: fill with 0 ---
    ft = ft.fillna(0)

    log.info("Feature table after imputation: %s  (NaN remaining: %d)",
             ft.shape, ft.isna().sum().sum())

    # Save
    ft.to_parquet(FEATURES_PATH)
    log.info("Saved improved feature table → %s", FEATURES_PATH)
    return ft


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    improve_imputation()
