"""
Phase 3 — Feature: Functional Dependency
Input : data/processed/depmap_atc_dependency.csv  (long: index=gene, cols=ModelID/CellLineName/chronos_score)
Outputs: DataFrame indexed by hgnc_symbol with columns:
    dep_mean_chronos   — mean Chronos score across ATC lines (lower = more essential)
    dep_min_chronos    — min Chronos score (most dependent line)
    dep_n_lines        — number of ATC lines screened for this gene
    dep_essential_frac — fraction of lines with Chronos < -0.5
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).parents[2]
DATA = ROOT / "data" / "processed"

ESSENTIAL_THRESH = -0.5


def compute(gene_universe: pd.DataFrame) -> pd.DataFrame:
    dep = pd.read_csv(DATA / "depmap_atc_dependency.csv", index_col=0)
    # index = hgnc_symbol (repeated per line), columns: ModelID, CellLineName, chronos_score
    dep.index.name = "hgnc_symbol"

    grp = dep.groupby("hgnc_symbol")["chronos_score"]
    result = pd.DataFrame({
        "dep_mean_chronos": grp.mean(),
        "dep_min_chronos": grp.min(),
        "dep_n_lines": grp.count(),
        "dep_essential_frac": grp.apply(lambda s: (s < ESSENTIAL_THRESH).mean()),
    })

    symbols = gene_universe["hgnc_symbol"]
    df = pd.DataFrame(index=pd.Index(symbols, name="hgnc_symbol")).join(result, how="left")
    # genes with no CRISPR data get NaN — intentional (not zero)
    return df
