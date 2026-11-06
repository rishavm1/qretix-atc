"""
scrna_qc.py
───────────
Phase 2.1: scRNA-seq QC and cell-type annotation for GSE193581.

Steps:
  1. Download GSE193581_RAW.tar + celltype_annotation from GEO FTP
  2. Extract per-sample 10x matrices
  3. Per-sample QC: min_genes, max_genes, pct_mito thresholds
  4. Doublet detection via scanpy score_genes heuristic (scrublet unavailable)
  5. Cell-type annotation: use published GSE193581_celltype_annotation.txt.gz
     as primary source; marker-gene fallback for unannotated cells
  6. Assign malignant_status via PAX8+ + EMT markers + published annotation
  7. Output: data/processed/GSE193581_annotated.h5ad

QC thresholds are logged per-sample in provenance JSON.
"""

from __future__ import annotations

import ftplib
import gzip
import io
import logging
import os
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from src.provenance import add_qc_thresholds, load as prov_load, record_path

log = logging.getLogger(__name__)

ROOT       = Path(__file__).resolve().parents[2]
RAW_DIR    = ROOT / "data" / "raw" / "GSE193581"
PROC_DIR   = ROOT / "data" / "processed"
OUTPUT_H5AD = PROC_DIR / "GSE193581_annotated.h5ad"

GEO_FTP_BASE = "/geo/series/GSE193nnn/GSE193581/suppl"
FTP_FILES = {
    "GSE193581_RAW.tar":                    "GSE193581_RAW.tar",
    "GSE193581_celltype_annotation.txt.gz": "GSE193581_celltype_annotation.txt.gz",
    "GSE193581_ATC_cell_line_bulk_count.txt.gz": "GSE193581_ATC_cell_line_bulk_count.txt.gz",
}

# ── Marker genes for fallback annotation ─────────────────────────────────────
MARKERS = {
    "Thyrocyte_differentiated": ["TG", "TPO", "SLC5A5", "TSHR", "PAX8", "NKX2-1"],
    "ATC_malignant_mesenchymal": ["VIM", "FN1", "CDH2", "ZEB1", "TWIST1", "SPP1"],
    "ATC_malignant_inflammatory": ["CXCL8", "IL6", "NFKBIA", "ICAM1"],
    "TAM_M2":    ["CD68", "CD163", "MRC1", "APOE", "C1QA", "TREM2"],
    "TAM_M1":    ["CD68", "CD86", "IL1B", "CXCL9"],
    "T_CD8":     ["CD3D", "CD8A", "GZMA", "GZMB", "PRF1"],
    "T_CD4":     ["CD3D", "CD4", "IL7R"],
    "Treg":      ["FOXP3", "IL2RA", "CTLA4"],
    "T_exhausted": ["PDCD1", "HAVCR2", "LAG3", "TIGIT"],
    "B_cell":    ["MS4A1", "CD19", "CD79A"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "CAF":       ["ACTA2", "COL1A1", "FAP", "PDGFRA"],
    "NK":        ["NCAM1", "NKG7", "GNLY"],
}

MALIGNANT_TYPES = {
    "ATC_malignant_mesenchymal", "ATC_malignant_inflammatory",
    "Thyrocyte_differentiated",   # PTC malignant cells
}


# ── FTP download ──────────────────────────────────────────────────────────────

def _ftp_download(fname: str, dest: Path) -> bool:
    if dest.exists() and dest.stat().st_size > 1000:
        log.info("%s already downloaded.", fname)
        return True
    log.info("FTP downloading %s …", fname)
    try:
        ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login()
        ftp.cwd(GEO_FTP_BASE)
        with open(dest, "wb") as f:
            ftp.retrbinary(f"RETR {fname}", f.write)
        ftp.quit()
        log.info("Saved %s → %.1f MB", fname, dest.stat().st_size / 1e6)
        return True
    except Exception as exc:
        log.error("FTP download failed for %s: %s", fname, exc)
        return False


def download_raw() -> bool:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    ok = True
    for fname in FTP_FILES:
        ok &= _ftp_download(fname, RAW_DIR / fname)
    return ok


# ── Extract 10x matrices ──────────────────────────────────────────────────────

def extract_matrices() -> Path:
    """Extract GSE193581_RAW.tar into RAW_DIR/matrices/"""
    mat_dir = RAW_DIR / "matrices"
    if mat_dir.exists() and any(mat_dir.iterdir()):
        log.info("Matrices already extracted.")
        return mat_dir
    mat_dir.mkdir(exist_ok=True)
    tar_path = RAW_DIR / "GSE193581_RAW.tar"
    log.info("Extracting %s …", tar_path.name)
    with tarfile.open(tar_path) as tar:
        tar.extractall(mat_dir)
    log.info("Extracted to %s", mat_dir)
    return mat_dir


# ── Load published cell-type annotations ─────────────────────────────────────

def load_published_annotations() -> pd.DataFrame | None:
    ann_path = RAW_DIR / "GSE193581_celltype_annotation.txt.gz"
    if not ann_path.exists():
        return None
    try:
        with gzip.open(ann_path, "rt") as f:
            df = pd.read_csv(f, sep="\t")
        log.info("Published annotations: %d cells, columns: %s", len(df), list(df.columns))
        return df
    except Exception as exc:
        log.warning("Could not load published annotations: %s", exc)
        return None


# ── Per-sample QC ─────────────────────────────────────────────────────────────

def _qc_sample(adata: sc.AnnData, sample_id: str) -> tuple[sc.AnnData, dict]:
    """Apply per-sample QC. Returns filtered AnnData + threshold dict."""
    n_before = adata.n_obs

    # Mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Adaptive thresholds based on distribution
    med_genes = float(np.median(adata.obs["n_genes_by_counts"]))
    mad_genes = float(np.median(np.abs(adata.obs["n_genes_by_counts"] - med_genes)))
    min_genes = max(200, int(med_genes - 3 * mad_genes))
    max_genes = min(8000, int(med_genes + 5 * mad_genes))

    med_mito = float(np.median(adata.obs["pct_counts_mt"]))
    max_mito = min(40.0, float(med_mito + 3 * np.median(
        np.abs(adata.obs["pct_counts_mt"] - med_mito)
    )) + 5.0)
    max_mito = max(max_mito, 20.0)  # floor at 20% — ATC has elevated mito

    # Filter
    mask = (
        (adata.obs["n_genes_by_counts"] >= min_genes) &
        (adata.obs["n_genes_by_counts"] <= max_genes) &
        (adata.obs["pct_counts_mt"] <= max_mito)
    )
    adata = adata[mask].copy()
    n_after = adata.n_obs

    # Doublet score via simple heuristic: cells with unusually high total counts
    # relative to their gene count (scrublet unavailable — MSVC build tools absent)
    count_per_gene = adata.obs["total_counts"] / (adata.obs["n_genes_by_counts"] + 1)
    threshold_99 = float(np.percentile(count_per_gene, 99))
    adata.obs["doublet_score"] = (count_per_gene / threshold_99).clip(0, 2)
    adata.obs["predicted_doublet"] = count_per_gene > threshold_99

    n_doublets = int(adata.obs["predicted_doublet"].sum())
    adata = adata[~adata.obs["predicted_doublet"]].copy()

    thresholds = {
        "sample_id": sample_id,
        "n_cells_before_qc": n_before,
        "n_cells_after_qc": adata.n_obs,
        "pct_removed": round(100 * (1 - adata.n_obs / n_before), 1),
        "min_genes": min_genes,
        "max_genes": max_genes,
        "max_pct_mito": round(max_mito, 1),
        "doublets_removed": n_doublets,
        "doublet_method": "count_per_gene_99pct_heuristic (scrublet unavailable: MSVC absent)",
    }
    log.info(
        "%s: %d → %d cells (%.1f%% removed, %d doublets)",
        sample_id, n_before, adata.n_obs,
        thresholds["pct_removed"], n_doublets,
    )
    return adata, thresholds


# ── Cell-type annotation ──────────────────────────────────────────────────────

def _annotate_celltypes(adata: sc.AnnData, published: pd.DataFrame | None) -> sc.AnnData:
    """
    Assign cell_type and malignant_status.
    Primary: published GSE193581 annotations matched by cell barcode.
    Fallback: marker-gene scoring for unannotated cells.
    """
    adata.obs["cell_type"] = "Unknown"
    adata.obs["malignant_status"] = "Unknown"
    adata.obs["annotation_source"] = "unassigned"

    # Try to match published annotations by barcode
    if published is not None:
        # Detect barcode and celltype columns
        bc_col = next((c for c in published.columns
                       if any(x in c.lower() for x in ["barcode", "cell", "id"])), None)
        ct_col = next((c for c in published.columns
                       if any(x in c.lower() for x in ["celltype", "cell_type", "type", "cluster", "annotation"])), None)

        if bc_col and ct_col:
            pub_map = dict(zip(published[bc_col].astype(str), published[ct_col].astype(str)))
            matched = 0
            for bc in adata.obs_names:
                # Try exact match, then suffix-stripped match
                ct = pub_map.get(bc) or pub_map.get(bc.split("-")[0])
                if ct:
                    adata.obs.loc[bc, "cell_type"] = ct
                    adata.obs.loc[bc, "annotation_source"] = "published_GSE193581"
                    matched += 1
            log.info("Published annotation matched %d / %d cells.", matched, adata.n_obs)

    # Marker-gene scoring for unmatched cells
    unmatched = adata.obs["cell_type"] == "Unknown"
    if unmatched.sum() > 0:
        log.info("Running marker-gene scoring for %d unannotated cells.", unmatched.sum())
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        for ct, markers in MARKERS.items():
            present = [g for g in markers if g in adata.var_names]
            if present:
                sc.tl.score_genes(adata, present, score_name=f"score_{ct}", use_raw=False)

        score_cols = [f"score_{ct}" for ct in MARKERS if f"score_{ct}" in adata.obs.columns]
        if score_cols:
            scores = adata.obs[score_cols]
            best = scores.idxmax(axis=1).str.replace("score_", "")
            mask_unmatched = adata.obs["cell_type"] == "Unknown"
            adata.obs.loc[mask_unmatched, "cell_type"] = best[mask_unmatched]
            adata.obs.loc[mask_unmatched, "annotation_source"] = "marker_gene_scoring"

    # Assign malignant_status
    adata.obs["malignant_status"] = adata.obs["cell_type"].apply(
        lambda ct: "malignant" if any(m in str(ct) for m in
                   ["ATC", "malignant", "Thyrocyte", "tumor", "cancer"]) else "non-malignant"
    )

    log.info(
        "Cell-type distribution:\n%s",
        adata.obs["cell_type"].value_counts().to_string()
    )
    log.info(
        "Malignant status: %s",
        adata.obs["malignant_status"].value_counts().to_dict()
    )
    return adata


# ── Main run ──────────────────────────────────────────────────────────────────

def run() -> None:
    PROC_DIR.mkdir(parents=True, exist_ok=True)

    if OUTPUT_H5AD.exists():
        log.info("GSE193581_annotated.h5ad already exists, skipping.")
        return

    # Download
    if not download_raw():
        log.error("GSE193581 download incomplete — cannot proceed with QC.")
        return

    # Extract
    mat_dir = extract_matrices()

    # Load published annotations
    published = load_published_annotations()

    # GSE193581 uses dense UMI count matrices: GSM*_UMI.txt.gz
    # Format: genes x cells tab-separated, first col = gene symbol
    umi_files = sorted(mat_dir.glob("GSM*_UMI.txt.gz"))
    log.info("Found %d UMI count files.", len(umi_files))

    all_adatas = []
    all_thresholds = []

    for umi_path in umi_files:
        # Parse sample ID and tissue type from filename
        # e.g. GSM5814583_ATC08_UMI.txt.gz -> sample_id=ATC08
        parts = umi_path.stem.replace(".txt", "").split("_")
        sample_id = parts[1] if len(parts) > 1 else umi_path.stem
        tissue = "ATC" if "ATC" in sample_id else ("PTC" if "PTC" in sample_id else "NORM")
        log.info("Processing %s (tissue=%s) …", sample_id, tissue)
        try:
            import scipy.sparse as sp
            df = pd.read_csv(umi_path, sep="\t", index_col=0, compression="gzip")
            # df: genes x cells
            adata = sc.AnnData(
                X=sp.csr_matrix(df.values.T),   # cells x genes
                obs=pd.DataFrame(index=df.columns),
                var=pd.DataFrame(index=df.index),
            )
            adata.obs["sample_id"] = sample_id
            adata.obs["tissue_type"] = tissue
            adata.obs["batch"] = sample_id
            adata.obs["platform"] = "10x_Genomics_3prime_v3"
            adata, thresholds = _qc_sample(adata, sample_id)
            all_adatas.append(adata)
            all_thresholds.append(thresholds)
        except Exception as exc:
            log.error("Failed to process %s: %s", sample_id, exc)

    if not all_adatas:
        log.error("No samples processed successfully.")
        return

    # Concatenate
    log.info("Concatenating %d samples …", len(all_adatas))
    combined = sc.concat(all_adatas, label="sample", keys=[a.obs["sample_id"].iloc[0] for a in all_adatas])

    # Annotate
    combined = _annotate_celltypes(combined, published)

    # Save
    combined.write_h5ad(OUTPUT_H5AD)
    log.info("Saved → %s (%d cells, %d genes)", OUTPUT_H5AD, combined.n_obs, combined.n_vars)

    # Update provenance with QC thresholds
    add_qc_thresholds("GSE193581", {
        "per_sample_thresholds": all_thresholds,
        "total_cells_before": sum(t["n_cells_before_qc"] for t in all_thresholds),
        "total_cells_after": combined.n_obs,
        "doublet_method": "count_per_gene_99pct_heuristic",
        "annotation_method": "published_GSE193581 + marker_gene_scoring fallback",
    })
