"""
merge_master_dataset.py
───────────────────────
Merges YOUR datasets + HF datasets into a single master dataset.

Strategy per source:
  - gene_universe    : HF wins (richer cols) but yours has more genes (41k vs 19k)
                       → merge both, keep all genes, fill extra HF cols where available
  - depmap           : keep BOTH, yours has more rows (124k vs 106k) different cell lines
                       → union on (hgnc_symbol + model_id), normalize column names
  - cbioportal       : HF only (you had 404) → copy as-is into data/processed/
  - opentargets      : HF has cleaner parquet (1293 targets, more score cols)
                       → replace your raw JSON-derived data with HF parquet
  - literature/trials: yours is flat & more usable → keep as-is

Output files (all go into data/processed/):
  - gene_universe.parquet          (merged, enriched)
  - depmap_atc_dependency.parquet  (merged, normalized)
  - cbioportal_atc_mutations.parquet (HF, copied)
  - cbioportal_atc_clinical.parquet  (HF, copied)
  - opentargets_atc_targets.parquet  (HF, copied)

HF root parquet files are deleted after merge.
"""

from __future__ import annotations
import logging
import shutil
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")

ROOT       = Path("C:/qretix-atc")
HF_DIR     = ROOT                          # HF parquets dropped in root
PROCESSED  = ROOT / "data" / "processed"
PROCESSED.mkdir(parents=True, exist_ok=True)

HF_FILES = [
    "atc_clinical_trials.parquet",
    "atc_pubmed_corpus.parquet",
    "cbioportal_atc_clinical.parquet",
    "cbioportal_atc_mutations.parquet",
    "depmap_atc_dependency.parquet",
    "depmap_atc_models.parquet",
    "gene_universe.parquet",
    "opentargets_atc_targets.parquet",
    "provenance_manifest.parquet",
]


# ── 1. Gene Universe ─────────────────────────────────────────────────────────
def merge_gene_universe():
    hf  = pd.read_parquet(HF_DIR / "gene_universe.parquet")
    my  = pd.read_parquet(PROCESSED / "gene_universe.parquet")

    # Normalize column names
    if "hgnc_symbol" in my.columns and "symbol" not in my.columns:
        my = my.rename(columns={"hgnc_symbol": "symbol", "entrez_gene_id": "entrez_id"})

    # Yours has 41k genes (includes non-protein-coding), HF has 19k (protein-coding only)
    # Start from yours (broader), then enrich with HF extra columns
    hf_extra = hf[["ensembl_gene_id", "name", "aliases", "previous_symbols", "uniprot_ids"]].copy()
    merged = my.merge(hf_extra, on="ensembl_gene_id", how="left")

    # Standardize column names
    if "symbol" not in merged.columns and "hgnc_symbol" in merged.columns:
        merged = merged.rename(columns={"hgnc_symbol": "symbol"})

    out = PROCESSED / "gene_universe.parquet"
    merged.to_parquet(out, index=False)
    log.info("gene_universe: %d genes, %d cols → %s", len(merged), len(merged.columns), out)
    return merged


# ── 2. DepMap ─────────────────────────────────────────────────────────────────
def merge_depmap():
    hf = pd.read_parquet(HF_DIR / "depmap_atc_dependency.parquet")
    my = pd.read_csv(PROCESSED / "depmap_atc_dependency.csv")

    # Normalize your columns to match HF schema
    my_norm = my.rename(columns={
        "ModelID": "model_id",
        "CellLineName": "cell_line",
        "chronos_score": "chronos_gene_effect",
    })
    # Add missing cols
    for col in ["ensembl_gene_id", "entrez_id", "platform", "batch"]:
        if col not in my_norm.columns:
            my_norm[col] = None

    # HF has ensembl IDs and entrez — fill into yours via gene universe
    gu = pd.read_parquet(PROCESSED / "gene_universe.parquet")
    sym_map = gu.set_index("symbol")[["ensembl_gene_id", "entrez_id"]].dropna()

    my_norm["ensembl_gene_id"] = my_norm["hgnc_symbol"].map(
        sym_map["ensembl_gene_id"]
    )
    my_norm["entrez_id"] = my_norm["hgnc_symbol"].map(sym_map["entrez_id"])
    my_norm["platform"] = "CRISPR-Cas9 (Chronos)"
    my_norm["batch"] = "DepMap release"

    # Normalize entrez_id to string in both to avoid type conflict
    hf["entrez_id"] = hf["entrez_id"].astype(str)
    my_norm["entrez_id"] = my_norm["entrez_id"].astype(str)

    # Union — combine both, deduplicate on model_id + hgnc_symbol
    combined = pd.concat([hf, my_norm[hf.columns]], ignore_index=True)
    before = len(combined)
    combined = combined.drop_duplicates(subset=["model_id", "hgnc_symbol"], keep="first")
    log.info("DepMap: %d rows → after dedup %d rows", before, len(combined))

    out = PROCESSED / "depmap_atc_dependency.parquet"
    combined.to_parquet(out, index=False)
    log.info("depmap_atc_dependency: %d rows, %d cols → %s", len(combined), len(combined.columns), out)

    # Also save models
    models_src = HF_DIR / "depmap_atc_models.parquet"
    models_dst = PROCESSED / "depmap_atc_models.parquet"
    shutil.copy2(models_src, models_dst)
    log.info("depmap_atc_models copied → %s", models_dst)

    return combined


# ── 3. cBioPortal (HF only — you had 404) ────────────────────────────────────
def copy_cbioportal():
    for fname in ["cbioportal_atc_mutations.parquet", "cbioportal_atc_clinical.parquet"]:
        src = HF_DIR / fname
        dst = PROCESSED / fname
        shutil.copy2(src, dst)
        df = pd.read_parquet(dst)
        log.info("%s: %d rows → %s", fname, len(df), dst)


# ── 4. Open Targets (HF is cleaner) ──────────────────────────────────────────
def copy_opentargets():
    src = HF_DIR / "opentargets_atc_targets.parquet"
    dst = PROCESSED / "opentargets_atc_targets.parquet"
    shutil.copy2(src, dst)
    df = pd.read_parquet(dst)
    log.info("opentargets_atc_targets: %d rows → %s", len(df), dst)


# ── 5. Provenance manifest ────────────────────────────────────────────────────
def copy_provenance_manifest():
    src = HF_DIR / "provenance_manifest.parquet"
    dst = PROCESSED / "provenance_manifest.parquet"
    shutil.copy2(src, dst)
    log.info("provenance_manifest copied → %s", dst)


# ── 6. Delete HF root files ───────────────────────────────────────────────────
def cleanup_hf_root():
    for fname in HF_FILES:
        p = HF_DIR / fname
        if p.exists():
            p.unlink()
            log.info("Deleted HF root file: %s", fname)


# ── Main ──────────────────────────────────────────────────────────────────────
def run():
    log.info("═══ Step 1: Gene Universe ═══")
    merge_gene_universe()

    log.info("═══ Step 2: DepMap ═══")
    merge_depmap()

    log.info("═══ Step 3: cBioPortal ═══")
    copy_cbioportal()

    log.info("═══ Step 4: Open Targets ═══")
    copy_opentargets()

    log.info("═══ Step 5: Provenance Manifest ═══")
    copy_provenance_manifest()

    log.info("═══ Step 6: Cleanup HF root files ═══")
    cleanup_hf_root()

    log.info("═══ DONE — Master dataset in data/processed/ ═══")

    # Print final summary
    print("\n=== MASTER DATASET SUMMARY ===")
    for f in sorted(PROCESSED.glob("*.parquet")):
        df = pd.read_parquet(f)
        print(f"  {f.name}: {len(df):,} rows × {len(df.columns)} cols")
    for f in sorted(PROCESSED.glob("*.csv")):
        df = pd.read_csv(f)
        print(f"  {f.name}: {len(df):,} rows × {len(df.columns)} cols")
    for f in sorted(PROCESSED.glob("*.h5ad")):
        print(f"  {f.name}: (h5ad — scRNA-seq, use scanpy to inspect)")


if __name__ == "__main__":
    run()
