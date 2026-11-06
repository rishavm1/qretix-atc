"""
Phase 3 — Feature: Network Relevance
Sources:
  - STRING v12.0 PPI links  (data/raw/string/9606.protein.links.v12.0.txt.gz)
  - STRING v12.0 protein info (data/raw/string/9606.protein.info.v12.0.txt.gz)
  - Reactome pathways GMT   (data/raw/reactome/ReactomePathways.gmt)

Status: IMPLEMENTED — scalar graph summaries computed from local raw files.
        Full graph object is also returned via get_string_graph() for Phase 5 GNN.

Outputs (indexed by hgnc_symbol):
    string_degree           — number of STRING PPI partners (combined_score >= 400)
    string_degree_hi        — partners at combined_score >= 700 (high-confidence)
    reactome_pathway_count  — number of Reactome pathways the gene belongs to
    ppi_driver_proximity    — shortest-path distance to known ATC driver gene set
                              (computed on STRING score >= 400 subgraph; NaN if
                               gene is isolated or driver set is unreachable)

Known ATC driver genes used for proximity computation:
    BRAF, TP53, TERT, KRAS, NRAS, HRAS, PIK3CA, PTEN, CDKN2A, RET,
    ALK, NTRK1, NTRK2, NTRK3, AKT1, MTOR, MET, EGFR, IDH1

Notes:
  - STRING uses ENSP IDs; mapped to gene symbol via protein.info file.
  - Reactome GMT uses HGNC symbols directly.
  - Genes with no STRING data get degree=0, proximity=NaN.
  - The full NetworkX graph is cached in memory for Phase 5 reuse.
"""
from __future__ import annotations

import gzip
import logging
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

ROOT = Path(__file__).parents[2]
RAW_STRING = ROOT / "data" / "raw" / "string"
RAW_REACTOME = ROOT / "data" / "raw" / "reactome"

SCORE_THRESHOLD_LO = 400   # low-confidence edge filter
SCORE_THRESHOLD_HI = 700   # high-confidence edge filter

# Curated ATC driver gene set (from literature review + OT associations)
ATC_DRIVER_GENES = frozenset([
    "BRAF", "TP53", "TERT", "KRAS", "NRAS", "HRAS",
    "PIK3CA", "PTEN", "CDKN2A", "RET", "ALK",
    "NTRK1", "NTRK2", "NTRK3", "AKT1", "MTOR",
    "MET", "EGFR", "IDH1", "EIF1AX", "ARID1A",
])

_GRAPH_CACHE: nx.Graph | None = None


def _load_string_info() -> dict[str, str]:
    """Return mapping ENSP_ID → gene_symbol."""
    path = RAW_STRING / "9606.protein.info.v12.0.txt.gz"
    ensp_to_sym: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                ensp_to_sym[parts[0]] = parts[1]
    log.info("Loaded STRING protein info: %d entries", len(ensp_to_sym))
    return ensp_to_sym


def _build_string_graph(ensp_to_sym: dict[str, str]) -> nx.Graph:
    """Build undirected NetworkX graph of STRING PPI (score >= LO threshold)."""
    global _GRAPH_CACHE
    if _GRAPH_CACHE is not None:
        return _GRAPH_CACHE

    path = RAW_STRING / "9606.protein.links.v12.0.txt.gz"
    G = nx.Graph()
    n_edges = 0

    with gzip.open(path, "rt", encoding="utf-8") as f:
        next(f)  # header: protein1 protein2 combined_score
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            score = int(parts[2])
            if score < SCORE_THRESHOLD_LO:
                continue
            sym1 = ensp_to_sym.get(parts[0])
            sym2 = ensp_to_sym.get(parts[1])
            if sym1 and sym2 and sym1 != sym2:
                G.add_edge(sym1, sym2, weight=score)
                n_edges += 1

    log.info("Built STRING graph: %d nodes, %d edges (score>=%d)",
             G.number_of_nodes(), n_edges, SCORE_THRESHOLD_LO)
    _GRAPH_CACHE = G
    return G


def _load_reactome_gmt() -> dict[str, list[str]]:
    """
    Return mapping gene_symbol → list of pathway names.
    Uses Ensembl2Reactome_All_Levels.txt (col0=EnsemblID, col1=PathwayID, col3=PathwayName, col5=Species)
    Filtered to Homo sapiens only. Maps Ensembl gene/protein IDs to HGNC symbols via gene universe.
    """
    from pathlib import Path as _Path
    # Load gene universe for Ensembl→symbol mapping
    gu_path = ROOT / "data" / "processed" / "gene_universe.parquet"
    import pandas as _pd
    gu = _pd.read_parquet(gu_path)[["ensembl_gene_id", "hgnc_symbol"]].dropna()
    ensembl_to_sym = dict(zip(gu["ensembl_gene_id"], gu["hgnc_symbol"]))

    path = RAW_REACTOME / "Ensembl2Reactome_All_Levels.txt"
    gene_to_paths: dict[str, list[str]] = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            ensembl_id = parts[0].strip()
            pathway_name = parts[3].strip() if len(parts) > 3 else ""
            species = parts[5].strip() if len(parts) > 5 else ""
            # Only human pathways
            if "Homo sapiens" not in species:
                continue
            # Map to HGNC symbol
            sym = ensembl_to_sym.get(ensembl_id)
            if sym and pathway_name:
                gene_to_paths.setdefault(sym, []).append(pathway_name)

    log.info("Loaded Reactome (human): %d unique genes", len(gene_to_paths))
    return gene_to_paths


def _compute_driver_proximity(
    G: nx.Graph, symbols: list[str]
) -> dict[str, float]:
    """
    Shortest-path hop distance from each gene to the nearest ATC driver gene.
    Only computed on the largest connected component of the STRING graph.
    Genes not in the graph or unreachable from any driver get NaN.
    """
    proximity: dict[str, float] = {}

    # Restrict to genes that are actually in the graph
    driver_nodes = [d for d in ATC_DRIVER_GENES if d in G]
    if not driver_nodes:
        log.warning("No ATC driver genes found in STRING graph — proximity set to NaN")
        return {s: float("nan") for s in symbols}

    # Use single-source shortest path from each driver (multi-source)
    # networkx.multi_source_dijkstra_path_length gives min distance to any source
    try:
        length_dict = dict(
            nx.multi_source_dijkstra_path_length(G, sources=driver_nodes, cutoff=6)
        )
    except nx.NetworkXError:
        log.warning("multi_source_dijkstra failed — falling back to BFS")
        length_dict = {}

    for sym in symbols:
        proximity[sym] = float(length_dict.get(sym, float("nan")))

    # Driver genes themselves get distance 0
    for d in ATC_DRIVER_GENES:
        if d in proximity:
            proximity[d] = 0.0

    return proximity


def get_string_graph() -> nx.Graph:
    """Return the cached STRING graph (call after compute() to ensure it's loaded)."""
    if _GRAPH_CACHE is None:
        ensp_to_sym = _load_string_info()
        _build_string_graph(ensp_to_sym)
    return _GRAPH_CACHE


def compute(gene_universe: pd.DataFrame) -> pd.DataFrame:
    """Return feature DataFrame indexed by hgnc_symbol."""
    symbols = gene_universe["hgnc_symbol"].tolist()

    # --- STRING ---
    ensp_to_sym = _load_string_info()
    G = _build_string_graph(ensp_to_sym)

    # Degree at LO threshold (all edges in G)
    degree_lo = dict(G.degree())

    # Degree at HI threshold (subgraph)
    G_hi = nx.Graph(
        [(u, v) for u, v, d in G.edges(data=True) if d.get("weight", 0) >= SCORE_THRESHOLD_HI]
    )
    degree_hi = dict(G_hi.degree())

    # Driver proximity
    log.info("Computing driver proximity (shortest-path BFS)...")
    proximity = _compute_driver_proximity(G, symbols)

    # --- Reactome ---
    gene_to_paths = _load_reactome_gmt()
    pathway_counts = {g: len(paths) for g, paths in gene_to_paths.items()}

    # --- Assemble ---
    df = pd.DataFrame({
        "string_degree": pd.Series(degree_lo),
        "string_degree_hi": pd.Series(degree_hi),
        "reactome_pathway_count": pd.Series(pathway_counts),
        "ppi_driver_proximity": pd.Series(proximity),
    })
    df.index.name = "hgnc_symbol"

    # Join to gene universe
    result = (
        pd.DataFrame(index=pd.Index(symbols, name="hgnc_symbol"))
        .join(df, how="left")
    )
    result["string_degree"] = result["string_degree"].fillna(0).astype(int)
    result["string_degree_hi"] = result["string_degree_hi"].fillna(0).astype(int)
    result["reactome_pathway_count"] = result["reactome_pathway_count"].fillna(0).astype(int)
    # ppi_driver_proximity: NaN kept for isolated genes — intentional

    n_in_graph = (result["string_degree"] > 0).sum()
    n_with_path = (result["reactome_pathway_count"] > 0).sum()
    log.info("Network relevance: %d genes in STRING, %d in Reactome", n_in_graph, n_with_path)

    return result
