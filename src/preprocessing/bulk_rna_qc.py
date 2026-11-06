"""
bulk_rna_qc.py
──────────────
Phase 2.2: Bulk microarray QC for GSE76039.

Platform: Affymetrix HG-U133 Plus 2.0 (GPL570)
37 tumors: 17 PDTC + 20 ATC

Steps:
  1. Parse GSE76039_family.soft.gz — extract expression matrix + sample metadata
  2. Probe-to-gene mapping via GPL570 platform annotation in SOFT
  3. Collapse multiple probes per gene: max-mean strategy
  4. Log2 transform if not already log-scaled
  5. Tag each sample with platform, batch, tissue_type (ATC vs PDTC)
  6. Output: data/processed/GSE76039_expression.parquet
             data/processed/GSE76039_metadata.parquet
"""

from __future__ import annotations

import gzip
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from src.provenance import add_qc_thresholds

log = logging.getLogger(__name__)

ROOT     = Path(__file__).resolve().parents[2]
RAW_DIR  = ROOT / "data" / "raw" / "GSE76039"
PROC_DIR = ROOT / "data" / "processed"

SOFT_PATH   = RAW_DIR / "GSE76039_family.soft.gz"
EXPR_OUT    = PROC_DIR / "GSE76039_expression.parquet"
META_OUT    = PROC_DIR / "GSE76039_metadata.parquet"


def _parse_soft() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Parse SOFT file. Returns:
      expr_df   : probes × samples (float)
      meta_df   : samples × metadata fields
      probe_map : {probe_id: gene_symbol}
    """
    log.info("Parsing %s …", SOFT_PATH.name)

    samples: dict[str, dict] = {}
    current_sample = None
    expr_rows: dict[str, dict] = {}   # probe_id -> {sample_id: value}
    probe_map: dict[str, str] = {}    # probe_id -> gene_symbol
    in_platform_table = False
    in_sample_table   = False
    platform_done     = False

    opener = gzip.open(SOFT_PATH, "rt", encoding="utf-8", errors="replace")

    for line in opener:
        line = line.rstrip()
        ll   = line.lower()

        # ── Platform annotation ───────────────────────────────────────────────
        if ll.startswith("!platform_table_begin"):
            in_platform_table = True
            continue
        if ll.startswith("!platform_table_end"):
            in_platform_table = False
            platform_done = True
            continue
        if in_platform_table:
            parts = line.split("\t")
            if len(parts) > 10 and parts[0] not in ("ID", ""):
                gene_sym = parts[10].strip()
                if gene_sym and gene_sym != "---":
                    primary = gene_sym.split("///")[0].strip()
                    if primary:
                        probe_map[parts[0]] = primary
            continue

        # ── Sample blocks ─────────────────────────────────────────────────────
        if line.startswith("^SAMPLE"):
            current_sample = line.split("=")[1].strip()
            samples[current_sample] = {"sample_id": current_sample}
            in_sample_table = False
            continue

        if ll.startswith("!sample_table_begin"):
            in_sample_table = True
            continue
        if ll.startswith("!sample_table_end"):
            in_sample_table = False
            continue

        if in_sample_table and current_sample:
            parts = line.split("\t")
            if len(parts) >= 2 and parts[0] not in ("ID_REF", ""):
                try:
                    val = float(parts[1])
                    if parts[0] not in expr_rows:
                        expr_rows[parts[0]] = {}
                    expr_rows[parts[0]][current_sample] = val
                except ValueError:
                    pass
            continue

        if current_sample and not in_sample_table:
            if line.startswith("!Sample_title"):
                samples[current_sample]["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1"):
                val = line.split("=", 1)[1].strip()
                if any(x in val.lower() for x in ["tissue", "type", "histology"]):
                    samples[current_sample]["tissue_type"] = val
            elif line.startswith("!Sample_source_name"):
                samples[current_sample]["source"] = line.split("=", 1)[1].strip()

    opener.close()

    log.info("Parsed: %d samples, %d probes, %d probe→gene mappings",
             len(samples), len(expr_rows), len(probe_map))

    expr_df = pd.DataFrame(expr_rows).T   # probes × samples
    meta_df = pd.DataFrame(samples).T

    return expr_df, meta_df, probe_map


def _collapse_probes(expr_df: pd.DataFrame, probe_map: dict) -> pd.DataFrame:
    """
    Map probes to gene symbols and collapse by max-mean strategy.
    Max-mean: for each gene, keep the probe with the highest mean expression.
    """
    expr_df = expr_df.copy()
    expr_df.index.name = "probe_id"
    expr_df = expr_df.reset_index()
    expr_df["gene_symbol"] = expr_df["probe_id"].map(probe_map)
    expr_df = expr_df.dropna(subset=["gene_symbol"])
    expr_df = expr_df[expr_df["gene_symbol"] != ""]

    sample_cols = [c for c in expr_df.columns if c not in ("probe_id", "gene_symbol")]
    expr_df[sample_cols] = expr_df[sample_cols].apply(pd.to_numeric, errors="coerce")

    # Max-mean: keep probe with highest mean across samples per gene
    expr_df["probe_mean"] = expr_df[sample_cols].mean(axis=1)
    expr_df = expr_df.sort_values("probe_mean", ascending=False)
    expr_df = expr_df.drop_duplicates(subset=["gene_symbol"], keep="first")
    expr_df = expr_df.set_index("gene_symbol")[sample_cols]

    log.info("After probe collapse: %d genes × %d samples", len(expr_df), len(sample_cols))
    return expr_df


def _check_log_scale(expr_df: pd.DataFrame) -> pd.DataFrame:
    """Log2-transform if data appears to be on linear scale (max > 100)."""
    max_val = expr_df.values.max()
    if max_val > 100:
        log.info("Data appears linear (max=%.1f) — applying log2(x+1).", max_val)
        expr_df = np.log2(expr_df + 1)
    else:
        log.info("Data appears log-scaled (max=%.2f) — no transform applied.", max_val)
    return expr_df


def _tag_metadata(meta_df: pd.DataFrame) -> pd.DataFrame:
    """Add platform, batch, and tissue_type columns."""
    meta_df = meta_df.copy()
    meta_df["platform"] = "Affymetrix_HGU133Plus2_GPL570"
    meta_df["batch"] = "GSE76039"

    # Infer ATC vs PDTC from title/source/tissue_type
    def _classify(row):
        text = " ".join(str(v) for v in row.values).lower()
        if "anaplastic" in text or " atc" in text:
            return "ATC"
        elif "poorly" in text or "pdtc" in text:
            return "PDTC"
        return "unknown"

    meta_df["tumor_type"] = meta_df.apply(_classify, axis=1)
    log.info("Tumor type distribution: %s", meta_df["tumor_type"].value_counts().to_dict())
    return meta_df


def run() -> None:
    PROC_DIR.mkdir(parents=True, exist_ok=True)

    if EXPR_OUT.exists() and META_OUT.exists():
        log.info("GSE76039 processed files already exist, skipping.")
        return

    if not SOFT_PATH.exists():
        log.error("GSE76039 SOFT file missing: %s", SOFT_PATH)
        return

    expr_df, meta_df, probe_map = _parse_soft()

    if expr_df.empty:
        log.error("Expression matrix empty after SOFT parse.")
        return

    # QC stats before collapse
    n_probes_before = len(expr_df)
    n_samples = len(meta_df)

    expr_df = _collapse_probes(expr_df, probe_map)
    expr_df = _check_log_scale(expr_df)
    meta_df = _tag_metadata(meta_df)

    # Save
    expr_df.to_parquet(EXPR_OUT)
    meta_df.to_parquet(META_OUT)
    log.info("GSE76039 expression → %s (%d genes × %d samples)", EXPR_OUT, len(expr_df), len(expr_df.columns))
    log.info("GSE76039 metadata   → %s", META_OUT)

    # Update provenance
    add_qc_thresholds("GSE76039", {
        "n_probes_raw": n_probes_before,
        "n_genes_after_collapse": len(expr_df),
        "n_samples": n_samples,
        "probe_collapse_method": "max_mean (keep probe with highest mean expression per gene)",
        "log_transform": "log2(x+1) if max>100 else none",
        "platform": "Affymetrix_HGU133Plus2_GPL570",
        "batch": "GSE76039",
        "tumor_type_distribution": meta_df["tumor_type"].value_counts().to_dict(),
        "note": "Pre-called expression values from GEO series matrix — not raw CEL files.",
    })
