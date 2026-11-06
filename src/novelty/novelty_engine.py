"""
Phase 6 — Novelty Engine
=========================
Computes FIVE independent novelty sub-scores per gene.
These are NEVER collapsed into a single novelty score.

Sub-scores:
  1. disease_novelty       — inverse of ATC-specific publication volume, recency-weighted
  2. target_novelty        — inverse of cross-cancer therapeutic investigation volume
  3. mechanistic_novelty   — is ATC mechanism previously described? (literature NLP)
  4. combination_novelty   — target+pathway combination proposed in any cancer?
  5. research_saturation   — publication trajectory slope (flat/declining = saturated)

Also implements the 8-point novelty validation checklist per candidate gene.

Output:
  - models/qretix_ranker/novelty_scores.parquet  (5 sub-scores per gene)
  - evaluation/reports/novelty_checklist.md       (8-point checklist results)
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
LIT_DIR   = ROOT / "literature"
RAW_OT    = ROOT / "data" / "raw" / "opentargets"
FEATURES_PATH = ROOT / "features" / "feature_table.parquet"
LABELS_PATH   = ROOT / "labels" / "labels.parquet"
QRETIX_DIR    = ROOT / "models" / "qretix_ranker"
EVAL_DIR      = ROOT / "evaluation" / "reports"

# Known ATC targets (Stage 3+) — used to assess "not already well-known"
ATC_APPROVED_TARGETS = frozenset([
    "BRAF", "MAP2K1", "MAP2K2", "KDR", "FLT1", "FLT4", "RET",
    "NTRK1", "NTRK2", "NTRK3", "PDCD1", "CD274", "CTLA4",
])

# Known ATC clinical targets (Stage 3)
ATC_CLINICAL_TARGETS = frozenset([
    "CDK4", "CDK6", "WEE1", "AURKA", "AURKB", "PLK1", "MTOR",
    "AKT1", "PIK3CA", "MDM2", "HDAC1", "HDAC2", "MET", "EGFR",
    "ALK", "ERBB2", "BCL2", "STAT3", "KRAS", "TERT",
])


def _load_pubmed() -> list[dict]:
    with open(LIT_DIR / "atc_pubmed_corpus.json", encoding="utf-8") as f:
        return json.load(f)


def _load_thyroid_pubmed() -> list[dict]:
    with open(LIT_DIR / "thyroid_pubmed_corpus.json", encoding="utf-8") as f:
        return json.load(f)


def _load_ot_assoc() -> dict[str, float]:
    """Return symbol → OT ATC score."""
    path = RAW_OT / "atc_associations.json"
    with open(path, encoding="utf-8") as f:
        records = json.load(f)
    return {
        r["target"]["approvedSymbol"]: float(r.get("score", 0.0))
        for r in records
        if "target" in r and "approvedSymbol" in r["target"]
    }


def _count_atc_mentions_by_year(pubmed: list[dict], symbols: list[str]) -> dict[str, dict[int, int]]:
    """Returns {symbol: {year: count}} using fast token intersection."""
    TOKEN_RE = re.compile(r'\b[A-Z][A-Z0-9]{2,14}\b')
    symbol_set = {s for s in symbols if len(s) >= 3 and s == s.upper()}
    by_year: dict[str, dict[int, int]] = {s: {} for s in symbols}

    for rec in pubmed:
        text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
        year_raw = rec.get("year") or rec.get("publication_date", "0")
        try:
            year = int(str(year_raw)[:4])
        except (ValueError, TypeError):
            year = 0
        if year < 1990 or year > 2030:
            continue

        hits = set(TOKEN_RE.findall(text)) & symbol_set
        for sym in hits:
            by_year[sym][year] = by_year[sym].get(year, 0) + 1

    return by_year


def _count_thyroid_mentions(thyroid_pubmed: list[dict], symbols: list[str]) -> dict[str, int]:
    """Count thyroid corpus mentions per gene using fast token intersection."""
    TOKEN_RE = re.compile(r'\b[A-Z][A-Z0-9]{2,14}\b')
    symbol_set = {s for s in symbols if len(s) >= 3 and s == s.upper()}
    counts = {s: 0 for s in symbols}
    for rec in thyroid_pubmed:
        text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
        hits = set(TOKEN_RE.findall(text)) & symbol_set
        for sym in hits:
            counts[sym] += 1
    return counts


def _disease_novelty_score(
    by_year: dict[str, dict[int, int]],
    symbols: list[str],
    current_year: int = 2026,
) -> pd.Series:
    """
    Disease novelty = inverse of ATC mention volume, weighted by recency.
    Novel genes have LOW publication counts in ATC corpus.
    Score ∈ [0, 1] where 1 = most novel (0 ATC papers).
    """
    scores = {}
    for sym in symbols:
        ydata = by_year.get(sym, {})
        total = sum(ydata.values())
        # Recency weight: papers in last 5 years count more
        recent = sum(v for y, v in ydata.items() if y >= current_year - 5)

        # Weighted count: older papers discounted
        weighted = total + recent  # double-count recent
        # Inverse sigmoid mapping to [0,1] (higher count = lower novelty)
        score = 1.0 / (1.0 + np.log1p(weighted))
        scores[sym] = score

    return pd.Series(scores, name="disease_novelty")


def _target_novelty_score(
    thy_counts: dict[str, int],
    ot_scores: dict[str, float],
    symbols: list[str],
) -> pd.Series:
    """
    Target novelty = inverse of cross-cancer therapeutic investigation.
    Uses thyroid corpus mention count + OT score as proxy.
    """
    scores = {}
    for sym in symbols:
        thy = thy_counts.get(sym, 0)
        ot  = ot_scores.get(sym, 0.0)
        # Combine: high OT score or high thyroid mentions = low target novelty
        combined = thy + ot * 100   # normalize OT to comparable scale
        score = 1.0 / (1.0 + np.log1p(combined))
        scores[sym] = score

    return pd.Series(scores, name="target_novelty")


def _mechanistic_novelty_score(pubmed: list[dict], symbols: list[str]) -> pd.Series:
    """Mechanistic novelty using fast token intersection."""
    TOKEN_RE = re.compile(r'\b[A-Z][A-Z0-9]{2,14}\b')
    symbol_set = {s for s in symbols if len(s) >= 3 and s == s.upper()}
    mech_terms = re.compile(
        r"\b(mechanism|pathway|signaling|crosstalk|downstream|upstream|cascade)\b",
        re.IGNORECASE,
    )
    mech_counts = {s: 0 for s in symbols}

    for rec in pubmed:
        text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
        if not mech_terms.search(text):
            continue
        hits = set(TOKEN_RE.findall(text)) & symbol_set
        for sym in hits:
            mech_counts[sym] += 1

    scores = {}
    for sym in symbols:
        n = mech_counts.get(sym, 0)
        score = max(0.0, 1.0 - min(n, 5) / 5.0)
        scores[sym] = score

    return pd.Series(scores, name="mechanistic_novelty")


def _combination_novelty_score(symbols: list[str]) -> pd.Series:
    """
    Combination novelty = whether this target + known ATC pathway has been
    proposed together. Approximated by: is the gene in a known ATC pathway
    (from feature_table reactome_pathway_count) AND has cross-modality support?
    Score ∈ [0,1].
    Note: Full combination novelty requires pathway-level literature search;
    this is a first-pass approximation using available features.
    """
    if not FEATURES_PATH.exists():
        return pd.Series(0.5, index=pd.Index(symbols, name="hgnc_symbol"),
                         name="combination_novelty")

    ft = pd.read_parquet(FEATURES_PATH)
    ft = ft.reindex(symbols)

    # Higher pathway count → more pathway context known → lower combination novelty
    path_count = ft.get("reactome_pathway_count", pd.Series(0, index=ft.index)).fillna(0)
    # Normalize: 0 pathways → most novel (score=1), many pathways → less novel
    max_count = path_count.quantile(0.99)
    scores = 1.0 - (path_count / max(max_count, 1)).clip(0, 1)
    scores.name = "combination_novelty"
    scores.index = pd.Index(symbols, name="hgnc_symbol")
    return scores


def _research_saturation_score(
    by_year: dict[str, dict[int, int]],
    symbols: list[str],
    window_start: int = 2019,
    window_end: int = 2024,
) -> pd.Series:
    """
    Research saturation = slope of publication rate over the last 5 years.
    Negative or flat slope → saturated (score = 0).
    Positive slope → growing (score = 1) → more novelty potential.
    Score ∈ [0, 1].
    """
    scores = {}
    years = list(range(window_start, window_end + 1))

    for sym in symbols:
        ydata = by_year.get(sym, {})
        counts = np.array([ydata.get(y, 0) for y in years], dtype=float)

        if counts.sum() == 0:
            # No papers → research hasn't started → very novel (high score)
            scores[sym] = 1.0
            continue

        # Linear regression slope
        x = np.arange(len(years), dtype=float)
        slope = np.polyfit(x, counts, 1)[0]

        # Positive slope → growing field → more opportunity (higher novelty value)
        # Negative slope → declining interest (but still novel from attention perspective)
        # Normalize slope to [0,1] using sigmoid
        score = 1.0 / (1.0 + np.exp(-slope))  # logistic of slope
        scores[sym] = float(score)

    return pd.Series(scores, name="research_saturation")


def _eight_point_checklist(symbol: str, ot_score: float) -> dict[str, bool]:
    """Implement 8-point novelty validation checklist for a given gene."""
    return {
        "not_in_OT_ATC_above_0.5": ot_score < 0.5,
        "not_approved_target_for_thyroid": symbol not in ATC_APPROVED_TARGETS,
        "not_ATC_clinical_target": symbol not in ATC_CLINICAL_TARGETS,
        "not_in_ATC_approved_label": symbol not in ATC_APPROVED_TARGETS,
        # Remaining checks require DisGeNET / COSMIC / patent search (marked as N/A without access)
        "disgenet_check": None,    # N/A — requires DisGeNET registration
        "cosmic_cgc_check": None,  # N/A — requires COSMIC license
        "not_primary_trial_target": None,  # N/A — requires manual trial review
        "not_in_systematic_review": None,  # N/A — requires manual literature review
    }


def compute_novelty_scores(gene_universe: pd.DataFrame) -> pd.DataFrame:
    """Main entry point — returns DataFrame with all 5 sub-scores."""
    symbols = gene_universe["hgnc_symbol"].tolist()

    log.info("Loading literature corpora for novelty analysis...")
    pubmed = _load_pubmed()
    thyroid_pubmed = _load_thyroid_pubmed()
    ot_scores = _load_ot_assoc()

    log.info("Counting ATC mentions by year (%d genes)...", len(symbols))
    by_year = _count_atc_mentions_by_year(pubmed, symbols)

    log.info("Counting thyroid corpus mentions...")
    thy_counts = _count_thyroid_mentions(thyroid_pubmed, symbols)

    log.info("Computing 5 novelty sub-scores...")
    dns = _disease_novelty_score(by_year, symbols)
    tns = _target_novelty_score(thy_counts, ot_scores, symbols)
    mns = _mechanistic_novelty_score(pubmed, symbols)
    cns = _combination_novelty_score(symbols)
    rss = _research_saturation_score(by_year, symbols)

    result = pd.DataFrame({
        "disease_novelty":      dns,
        "target_novelty":       tns,
        "mechanistic_novelty":  mns,
        "combination_novelty":  cns,
        "research_saturation":  rss,
    }, index=pd.Index(symbols, name="hgnc_symbol"))

    return result


def run():
    """Phase 6 entry point."""
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")

    gene_universe = pd.read_parquet(ROOT / "data" / "processed" / "gene_universe.parquet")
    print(f"Gene universe: {len(gene_universe):,} genes")

    novelty_df = compute_novelty_scores(gene_universe)

    # Save to both locations — pipeline validates features/novelty_scores.parquet
    QRETIX_DIR.mkdir(parents=True, exist_ok=True)
    out_path = QRETIX_DIR / "novelty_scores.parquet"
    novelty_df.to_parquet(out_path)

    features_path = ROOT / "features" / "novelty_scores.parquet"
    novelty_df.to_parquet(features_path)
    print(f"Wrote {out_path}  shape={novelty_df.shape}")
    print(f"Wrote {features_path}  (pipeline gate copy)")

    # Print summary stats per sub-score
    print("\nNovelty sub-score summaries (higher = more novel):")
    for col in novelty_df.columns:
        s = novelty_df[col]
        print(f"  {col}: mean={s.mean():.3f} | median={s.median():.3f} | "
              f"top10%>{s.quantile(0.9):.3f}")

    # 8-point checklist for top candidates (if QRETIX scores available)
    _write_checklist_report(gene_universe, novelty_df)

    return novelty_df


def _write_checklist_report(gene_universe: pd.DataFrame, novelty_df: pd.DataFrame):
    """Generate 8-point novelty checklist for top candidates."""
    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    ot_scores = _load_ot_assoc()

    # Determine top candidates
    qretix_path = QRETIX_DIR / "qretix_scores.parquet"
    if qretix_path.exists():
        qretix = pd.read_parquet(qretix_path)
        top_genes = qretix.head(50).index.tolist()
    else:
        # Fallback: top genes by average novelty
        top_genes = novelty_df.mean(axis=1).nlargest(50).index.tolist()

    lines = [
        "# Novelty Checklist Report — Phase 6",
        "",
        "8-point novelty validation for top 50 candidates.",
        "N/A = requires external registration (COSMIC/DisGeNET/manual review).",
        "",
        "| Gene | OT<0.5 | Not Approved | Not Clinical | DisGeNET | COSMIC | Trial | SysRev |",
        "|------|--------|-------------|-------------|----------|--------|-------|--------|",
    ]

    for sym in top_genes:
        ot = ot_scores.get(sym, 0.0)
        chk = _eight_point_checklist(sym, ot)
        def _fmt(v):
            if v is None: return "N/A"
            return "✅" if v else "❌"
        lines.append(
            f"| {sym} | {_fmt(chk['not_in_OT_ATC_above_0.5'])} | "
            f"{_fmt(chk['not_approved_target_for_thyroid'])} | "
            f"{_fmt(chk['not_ATC_clinical_target'])} | "
            f"{_fmt(chk['disgenet_check'])} | "
            f"{_fmt(chk['cosmic_cgc_check'])} | "
            f"{_fmt(chk['not_primary_trial_target'])} | "
            f"{_fmt(chk['not_in_systematic_review'])} |"
        )

    path = EVAL_DIR / "novelty_checklist.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {path}")


if __name__ == "__main__":
    run()
