"""
phase2_run.py
─────────────
Phase 2 orchestrator: QC, preprocessing, batch tagging, QC report.

Run: python -m src.phase2_run
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

log = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]


def _step(name: str, fn) -> bool:
    log.info("════ %s ════", name)
    try:
        fn()
        log.info("✓ %s complete.", name)
        return True
    except Exception as exc:
        log.error("✗ %s FAILED: %s", name, exc)
        _log_decision(f"Phase 2 step '{name}' failed: {exc}")
        return False


def _log_decision(msg: str) -> None:
    log_path = ROOT / "docs" / "decisions_log.md"
    entry = (
        f"\n## {datetime.now(timezone.utc).strftime('%Y-%m-%d')} [Phase 2] — Runtime Gap\n"
        f"**Decision:** {msg}\n"
        f"**Logged by:** phase2_run.py\n\n---\n"
    )
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(entry)


def generate_qc_report() -> None:
    """Build docs/qc_report.md from provenance QC thresholds + processed file stats."""
    import pandas as pd

    prov_dir = ROOT / "data" / "provenance"
    proc_dir = ROOT / "data" / "processed"
    out_path = ROOT / "docs" / "qc_report.md"

    lines = [
        "# QRETIX-ATC v0.1 — QC Report",
        "",
        f"_Generated: {datetime.now(timezone.utc).isoformat()}_",
        "",
        "## Summary",
        "",
    ]

    # ── GSE193581 scRNA-seq ───────────────────────────────────────────────────
    lines += ["## GSE193581 (scRNA-seq)", ""]
    h5ad = proc_dir / "GSE193581_annotated.h5ad"
    if h5ad.exists():
        try:
            import scanpy as sc
            adata = sc.read_h5ad(h5ad)
            lines += [
                f"- Cells after QC: **{adata.n_obs}**",
                f"- Genes: {adata.n_vars}",
                f"- Cell types: {adata.obs['cell_type'].nunique()}",
                f"- Malignant cells: {(adata.obs['malignant_status']=='malignant').sum()}",
                f"- Non-malignant cells: {(adata.obs['malignant_status']=='non-malignant').sum()}",
                "",
                "### Per-sample QC thresholds",
                "",
                "| Sample | Cells in | Cells out | % removed | min_genes | max_genes | max_mito% | doublets |",
                "|--------|----------|-----------|-----------|-----------|-----------|-----------|----------|",
            ]
            prov = json.loads((prov_dir / "GSE193581.json").read_text())
            thresholds = prov.get("qc_thresholds_applied", {})
            for t in thresholds.get("per_sample_thresholds", []):
                lines.append(
                    f"| {t.get('sample_id','')} "
                    f"| {t.get('n_cells_before_qc','')} "
                    f"| {t.get('n_cells_after_qc','')} "
                    f"| {t.get('pct_removed','')}% "
                    f"| {t.get('min_genes','')} "
                    f"| {t.get('max_genes','')} "
                    f"| {t.get('max_pct_mito','')} "
                    f"| {t.get('doublets_removed','')} |"
                )
        except Exception as exc:
            lines += [f"- h5ad exists but could not read: {exc}", ""]
    else:
        lines += ["- **NOT YET PROCESSED** — GSE193581_annotated.h5ad missing.", ""]

    # ── GSE76039 bulk ─────────────────────────────────────────────────────────
    lines += ["", "## GSE76039 (Bulk Microarray)", ""]
    expr = proc_dir / "GSE76039_expression.parquet"
    meta = proc_dir / "GSE76039_metadata.parquet"
    if expr.exists() and meta.exists():
        try:
            e = pd.read_parquet(expr)
            m = pd.read_parquet(meta)
            prov = json.loads((prov_dir / "GSE76039.json").read_text())
            thresholds = prov.get("qc_thresholds_applied", {})
            lines += [
                f"- Samples in: **{thresholds.get('n_samples', len(m.columns))}**",
                f"- Samples out: **{len(m)}** (no sample-level exclusions for bulk array)",
                f"- % removed: 0%",
                f"- Probes raw: {thresholds.get('n_probes_raw', '?')}",
                f"- Genes after probe collapse: **{len(e)}**",
                f"- Probe collapse method: {thresholds.get('probe_collapse_method', 'max_mean')}",
                f"- Log transform: {thresholds.get('log_transform', '?')}",
                f"- Tumor type distribution: {thresholds.get('tumor_type_distribution', '?')}",
                f"- Platform: {thresholds.get('platform', 'Affymetrix_HGU133Plus2_GPL570')}",
                f"- Note: {thresholds.get('note', '')}",
            ]
        except Exception as exc:
            lines += [f"- Error reading processed files: {exc}"]
    else:
        lines += ["- **NOT YET PROCESSED** — GSE76039 parquet files missing."]

    # ── Genomics ──────────────────────────────────────────────────────────────
    lines += ["", "## Genomics (DepMap ATC Cell Lines)", ""]
    for fname, label in [("atc_mutations.parquet", "Mutations"),
                          ("atc_cna.parquet", "CNA")]:
        fpath = proc_dir / fname
        if fpath.exists():
            df = pd.read_parquet(fpath)
            lines += [
                f"### {label}",
                f"- Rows: {len(df)}",
                f"- Genes: {df['hgnc_symbol'].nunique() if 'hgnc_symbol' in df.columns else '?'}",
                f"- Source: DepMap precalled (not raw VCFs) — cBioPortal ATC cohorts unavailable (404)",
                "",
            ]
        else:
            lines += [f"### {label}", "- NOT YET PROCESSED", ""]

    # ── Batch tagging ─────────────────────────────────────────────────────────
    lines += ["", "## Batch Tagging", ""]
    lines += ["All processed files audited for platform/batch columns by batch_tagger.py.", ""]

    # ── Data gaps ─────────────────────────────────────────────────────────────
    lines += [
        "## Documented Data Gaps",
        "",
        "| Gap | Reason | Impact |",
        "|-----|--------|--------|",
        "| cBioPortal GATCI 2024 | Study ID thca_gatci_cellrep_2024 returned 404 | Using DepMap somatic mutations as proxy |",
        "| cBioPortal MSKCC 2016 | Study ID thca_mskcc_2016 returned 404 | Same proxy |",
        "| scrublet doublet detection | annoy build requires MSVC (absent) | Using count_per_gene_99pct heuristic |",
        "| ENCODE ATC epigenomics | No ATC biosamples in ENCODE | Normal thyroid reference only |",
        "",
        "---",
        "_All thresholds logged per-dataset in data/provenance/*.json_",
    ]

    out_path.write_text("\n".join(lines), encoding="utf-8")
    log.info("QC report → %s", out_path)


def run() -> None:
    from src.preprocessing import (
        bulk_rna_qc, genomics_qc, batch_tagger, scrna_qc,
    )

    steps = [
        ("2.1 scRNA-seq QC (GSE193581)",  scrna_qc.run),
        ("2.2 Bulk RNA QC (GSE76039)",    bulk_rna_qc.run),
        ("2.3 Genomics QC",               genomics_qc.run),
        ("2.4 Batch Tagging",             batch_tagger.run),
        ("QC Report",                     generate_qc_report),
    ]

    results = {}
    for name, fn in steps:
        results[name] = _step(name, fn)

    log.info("\n════ Phase 2 Summary ════")
    for name, ok in results.items():
        log.info("  %s  %s", "✓" if ok else "✗", name)

    failed = [n for n, ok in results.items() if not ok]
    if failed:
        log.warning("Phase 2 completed with %d failures: %s", len(failed), failed)
    else:
        log.info("Phase 2 complete — all steps succeeded.")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )
    run()
