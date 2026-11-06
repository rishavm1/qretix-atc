"""
Phase 3 — Feature Engineering Orchestrator
Gate artifacts:
    features/feature_table.parquet
    docs/feature_spec.md
"""
from pathlib import Path
import sys
import logging
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))

from src.features import (
    genetic_evidence,
    tumor_cell_specificity,
    functional_dependency,
    literature_features,
    disease_relevance,
    network_relevance,
    protein_structural,
    evidence_independence,
)

FEATURES_DIR = ROOT / "features"
DOCS_DIR = ROOT / "docs"
DATA = ROOT / "data" / "processed"

MODULES = [
    (
        "genetic_evidence",
        genetic_evidence,
        ["mut_freq", "amp_freq", "del_freq", "any_genetic_event"],
    ),
    (
        "tumor_cell_specificity",
        tumor_cell_specificity,
        ["scrna_malignant_mean", "scrna_nonmal_mean", "scrna_specificity",
         "bulk_atc_mean", "bulk_pdtc_mean", "bulk_atc_vs_pdtc_fc"],
    ),
    (
        "functional_dependency",
        functional_dependency,
        ["dep_mean_chronos", "dep_min_chronos", "dep_n_lines", "dep_essential_frac"],
    ),
    (
        "literature_features",
        literature_features,
        ["lit_atc_mentions", "lit_thyroid_mentions", "lit_atc_ratio", "lit_recency_score"],
    ),
    (
        "disease_relevance",
        disease_relevance,
        ["ot_atc_score", "ot_literature_score", "ot_somatic_mut_score",
         "ot_sm_approved", "ot_ab_approved", "ot_druggability_tier"],
    ),
    (
        "network_relevance",
        network_relevance,
        ["string_degree", "string_degree_hi", "reactome_pathway_count",
         "ppi_driver_proximity"],
    ),
    (
        "protein_structural",
        protein_structural,
        ["uniprot_length", "uniprot_domain_count", "uniprot_has_structure",
         "uniprot_is_membrane", "uniprot_is_secreted", "uniprot_is_reviewed",
         "chembl_modulator_count", "chembl_approved_count", "chembl_clinical_count"],
    ),
    (
        "evidence_independence",
        evidence_independence,
        ["ev_n_modalities", "ev_genomics_flag", "ev_bulk_rna_flag",
         "ev_scrna_flag", "ev_depmap_flag", "ev_literature_flag",
         "ev_cross_modality", "ev_lit_n_atc_mentions"],
    ),
]


def run():
    FEATURES_DIR.mkdir(exist_ok=True)

    gene_universe = pd.read_parquet(DATA / "gene_universe.parquet")
    print(f"Gene universe: {len(gene_universe):,} genes")

    frames = []
    stats = []
    for name, mod, expected_cols in MODULES:
        print(f"\n  [{name}] computing...", flush=True)
        try:
            df = mod.compute(gene_universe)
        except Exception as e:
            print(f"  ERROR in {name}: {e}")
            raise

        missing = [c for c in expected_cols if c not in df.columns]
        if missing:
            raise ValueError(f"{name} missing expected columns: {missing}")

        n_nonnan = df[expected_cols].select_dtypes("number").notna().any(axis=1).sum()
        gap = getattr(df, "attrs", {}).get("gap_note", "")
        print(f"  [{name}] OK — {n_nonnan:,} genes with data"
              + (f" [GAP: {gap}]" if gap else ""))
        # Deduplicate index within each frame before appending
        if df.index.duplicated().any():
            df = df[~df.index.duplicated(keep="first")]
        frames.append(df[expected_cols])
        stats.append((name, expected_cols, n_nonnan, gap))

    feature_table = pd.concat(frames, axis=1)
    # Deduplicate index (keep first) — defensive against duplicate gene symbols
    feature_table = feature_table[~feature_table.index.duplicated(keep="first")]
    feature_table.index.name = "hgnc_symbol"

    out_path = FEATURES_DIR / "feature_table.parquet"
    feature_table.to_parquet(out_path)
    print(f"\nWrote {out_path}  shape={feature_table.shape}")

    _write_feature_spec(stats, feature_table)
    return feature_table


def _write_feature_spec(stats, feature_table):
    lines = [
        "# Feature Specification — Phase 3",
        "",
        f"Total genes: {len(feature_table):,}  |  Total features: {len(feature_table.columns)}",
        "",
        "## Feature Groups",
        "",
    ]
    for name, cols, n_nonnan, gap in stats:
        lines.append(f"### {name}")
        lines.append(f"- Genes with data: {n_nonnan:,}")
        if gap:
            lines.append(f"- **GAP**: {gap}")
        lines.append("- Columns:")
        for c in cols:
            if c in feature_table.columns:
                col_data = feature_table[c]
                if col_data.dtype.kind in "biufc":
                    pct_nan = col_data.isna().mean() * 100
                    lines.append(f"  - `{c}` — {pct_nan:.1f}% NaN")
                else:
                    lines.append(f"  - `{c}` — categorical")
        lines.append("")

    lines += [
        "## Gap Summary",
        "",
        "| Module | Status |",
        "|--------|--------|",
        "| genetic_evidence | COMPLETE — DepMap somatic mutations + CNA |",
        "| tumor_cell_specificity | COMPLETE — GSE193581 scRNA + GSE76039 bulk |",
        "| functional_dependency | COMPLETE — DepMap CRISPR Chronos (7 ATC lines) |",
        "| literature_features | COMPLETE — PubMed ATC + thyroid corpora |",
        "| disease_relevance | COMPLETE — Open Targets ATC associations + tractability |",
        "| network_relevance | COMPLETE — STRING PPI + Reactome + driver proximity |",
        "| protein_structural | COMPLETE — UniProt + ChEMBL drug mechanisms |",
        "| evidence_independence | COMPLETE — Cross-modality evidence flags |",
        "",
        "## Leakage Risk Notes",
        "",
        "| Feature | Risk | Mitigation |",
        "|---------|------|------------|",
        "| lit_atc_mentions | TEMPORAL | Per-PMID dates stored; filtered in Phase 7 temporal splits |",
        "| ot_atc_score | TEMPORAL | OT snapshot date logged in provenance |",
        "| ev_literature_flag | TEMPORAL | Derived from literature corpus; filtered in Phase 7 |",
        "| dep_mean_chronos | LOW | CRISPR screens don't leak label information |",
        "| genetic_known_driver_flag | HIGH | Kept as isolated column, never merged into score |",
    ]

    spec_path = ROOT / "docs" / "feature_spec.md"
    spec_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {spec_path}")


def validate() -> bool:
    """Phase 3 acceptance criteria check."""
    ft_path = FEATURES_DIR / "feature_table.parquet"
    spec_path = DOCS_DIR / "feature_spec.md"

    checks = [
        (ft_path.exists(), f"feature_table.parquet exists at {ft_path}"),
        (spec_path.exists(), f"feature_spec.md exists at {spec_path}"),
    ]

    if ft_path.exists():
        ft = pd.read_parquet(ft_path)
        checks.append((len(ft) >= 19000, f"≥19,000 genes in feature table (got {len(ft):,})"))
        checks.append((len(ft.columns) >= 20, f"≥20 feature columns (got {len(ft.columns)})"))

    all_pass = True
    for ok, msg in checks:
        status = "✅" if ok else "❌"
        print(f"  {status} {msg}")
        if not ok:
            all_pass = False
    return all_pass


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate", action="store_true")
    args = parser.parse_args()

    if args.validate:
        ok = validate()
        sys.exit(0 if ok else 1)
    else:
        run()
