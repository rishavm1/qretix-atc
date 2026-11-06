"""
Phase 3 — Feature: Literature Features
Inputs : literature/atc_pubmed_corpus.json, literature/thyroid_pubmed_corpus.json
Outputs: DataFrame indexed by hgnc_symbol with columns:
    lit_atc_mentions      — number of ATC-corpus abstracts mentioning the gene symbol
    lit_thyroid_mentions  — number of thyroid-corpus abstracts mentioning the gene symbol
    lit_atc_ratio         — lit_atc_mentions / (lit_thyroid_mentions + 1)
    lit_recency_score     — fraction of ATC mentions in papers published >= 2019

Performance: uses word-tokenization instead of per-gene regex to achieve
O(papers × words) instead of O(papers × genes).
"""
from pathlib import Path
import json
import re
import pandas as pd

ROOT = Path(__file__).parents[2]
LIT = ROOT / "literature"

RECENCY_YEAR = 2019
# Minimum symbol length to avoid noise (e.g. 'OR', 'IF', 'GO')
MIN_SYM_LEN = 3


def _count_mentions_fast(records: list[dict], symbol_set: set[str]) -> tuple[dict, dict]:
    """
    Scan each record once: extract all uppercase word tokens,
    intersect with symbol_set. O(papers × avg_tokens) instead of O(papers × genes).
    """
    counts: dict[str, int] = {s: 0 for s in symbol_set}
    recent: dict[str, int] = {s: 0 for s in symbol_set}

    TOKEN_RE = re.compile(r'\b[A-Z][A-Z0-9]{2,14}\b')

    for rec in records:
        text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
        # Extract all uppercase tokens that could be gene symbols
        tokens_in_paper = set(TOKEN_RE.findall(text))
        # Intersect with our gene symbol set
        hits = tokens_in_paper & symbol_set
        year_raw = rec.get("year") or rec.get("publication_date", "0")
        try:
            year = int(str(year_raw)[:4])
        except (ValueError, TypeError):
            year = 0

        for sym in hits:
            counts[sym] += 1
            if year >= RECENCY_YEAR:
                recent[sym] += 1

    return counts, recent


def compute(gene_universe: pd.DataFrame) -> pd.DataFrame:
    symbols = gene_universe["hgnc_symbol"].tolist()
    # Only try to match symbols that are uppercase and >= MIN_SYM_LEN
    symbol_set = {s for s in symbols if len(s) >= MIN_SYM_LEN and s == s.upper()}

    with open(LIT / "atc_pubmed_corpus.json", encoding="utf-8") as f:
        atc_records = json.load(f)
    with open(LIT / "thyroid_pubmed_corpus.json", encoding="utf-8") as f:
        thy_records = json.load(f)

    atc_counts, atc_recent = _count_mentions_fast(atc_records, symbol_set)
    thy_counts, _ = _count_mentions_fast(thy_records, symbol_set)

    # Fill zeros for symbols that were excluded (lowercase/short)
    for s in symbols:
        atc_counts.setdefault(s, 0)
        atc_recent.setdefault(s, 0)
        thy_counts.setdefault(s, 0)

    df = pd.DataFrame({
        "lit_atc_mentions":     pd.Series(atc_counts),
        "lit_thyroid_mentions": pd.Series(thy_counts),
    })
    df.index.name = "hgnc_symbol"
    df["lit_atc_ratio"] = df["lit_atc_mentions"] / (df["lit_thyroid_mentions"] + 1)
    df["lit_recency_score"] = (
        pd.Series(atc_recent) / df["lit_atc_mentions"].replace(0, 1)
    )

    result = (
        pd.DataFrame(index=pd.Index(symbols, name="hgnc_symbol"))
        .join(df, how="left")
        .fillna(0.0)
    )
    return result
