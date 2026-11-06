"""
Phase 3 — Feature: Genetic Evidence
Inputs : data/processed/atc_mutations.parquet, data/processed/atc_cna.parquet
Outputs: DataFrame indexed by hgnc_symbol with columns:
    mut_freq          — fraction of ATC lines with ≥1 somatic mutation
    amp_freq          — fraction of ATC lines with log2_copy_ratio > 1.0 (amplification)
    del_freq          — fraction of ATC lines with log2_copy_ratio < -1.0 (deep deletion)
    any_genetic_event — mut_freq + amp_freq + del_freq (capped at 1.0)
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).parents[2]
DATA = ROOT / "data" / "processed"

AMP_THRESH = 1.0   # log2 copy ratio
DEL_THRESH = -1.0


def compute(gene_universe: pd.DataFrame) -> pd.DataFrame:
    symbols = gene_universe["hgnc_symbol"]

    mut = pd.read_parquet(DATA / "atc_mutations.parquet")
    cna = pd.read_parquet(DATA / "atc_cna.parquet")

    n_mut_lines = mut["ModelID"].nunique()
    n_cna_lines = cna["ModelID"].nunique()

    # mutation frequency
    mut_freq = (
        mut.groupby("hgnc_symbol")["ModelID"].nunique() / n_mut_lines
    ).rename("mut_freq")

    # CNA frequencies
    amp = cna[cna["log2_copy_ratio"] > AMP_THRESH]
    del_ = cna[cna["log2_copy_ratio"] < DEL_THRESH]
    amp_freq = (amp.groupby("hgnc_symbol")["ModelID"].nunique() / n_cna_lines).rename("amp_freq")
    del_freq = (del_.groupby("hgnc_symbol")["ModelID"].nunique() / n_cna_lines).rename("del_freq")

    df = (
        pd.DataFrame(index=symbols)
        .join(mut_freq, how="left")
        .join(amp_freq, how="left")
        .join(del_freq, how="left")
        .fillna(0.0)
    )
    df["any_genetic_event"] = (df["mut_freq"] + df["amp_freq"] + df["del_freq"]).clip(upper=1.0)
    df.index.name = "hgnc_symbol"
    return df
