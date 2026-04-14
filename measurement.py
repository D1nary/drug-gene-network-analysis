from __future__ import annotations

from pathlib import Path

import networkx as nx
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CBP_INPUT_DIR = PROJECT_ROOT / "data_for_cbp"
DEFAULT_CBP_OUTPUT_DIR = PROJECT_ROOT / "results" / "cbp_measurement"

DEFAULT_EDGE_LIST = DEFAULT_CBP_INPUT_DIR / "edge_list.csv"

REQUIRED_EDGE_COLUMNS = {"drug_1", "drug_2"}


def _validate_columns(df: pd.DataFrame, required: set[str], file_path: Path) -> None:
    missing = required - set(df.columns)
    if missing:
        missing_str = ", ".join(sorted(missing))
        raise ValueError(f"{file_path} is missing required columns: {missing_str}")


def _build_graph(edge_df: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()
    has_weight = "similarity" in edge_df.columns
    for row in edge_df.itertuples(index=False):
        if has_weight:
            graph.add_edge(str(row.drug_1), str(row.drug_2), weight=float(row.similarity))
        else:
            graph.add_edge(str(row.drug_1), str(row.drug_2))
    return graph


def compute_betweenness_centrality(
    graph: nx.Graph,
    output_path: Path,
    weight: str | None = "weight",
) -> pd.DataFrame:
    """Compute betweenness centrality for all nodes and save results to CSV."""
    scores = nx.betweenness_centrality(graph, weight=weight, normalized=True)
    df = (
        pd.DataFrame(scores.items(), columns=["node", "betweenness_centrality"])
        .sort_values("betweenness_centrality", ascending=False)
        .reset_index(drop=True)
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    return df


def compute_closeness_centrality(
    graph: nx.Graph,
    output_path: Path,
) -> pd.DataFrame:
    """Compute closeness centrality for all nodes and save results to CSV."""
    scores = nx.closeness_centrality(graph)
    df = (
        pd.DataFrame(scores.items(), columns=["node", "closeness_centrality"])
        .sort_values("closeness_centrality", ascending=False)
        .reset_index(drop=True)
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    return df


def compute_pagerank(
    graph: nx.Graph,
    output_path: Path,
    weight: str | None = "weight",
    alpha: float = 0.85,
) -> pd.DataFrame:
    """Compute PageRank for all nodes and save results to CSV."""
    scores = nx.pagerank(graph, alpha=alpha, weight=weight)
    df = (
        pd.DataFrame(scores.items(), columns=["node", "pagerank"])
        .sort_values("pagerank", ascending=False)
        .reset_index(drop=True)
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    return df


def _extract_giant_component(graph: nx.Graph) -> nx.Graph:
    """Return the largest connected component of *graph* as a subgraph view."""
    largest_cc = max(nx.connected_components(graph), key=len)
    return graph.subgraph(largest_cc).copy()


def compute_cbp_measurements(
    edge_list_path: Path = DEFAULT_EDGE_LIST,
    output_dir: Path = DEFAULT_CBP_OUTPUT_DIR,
    giant_component_only: bool = True,
) -> dict[str, pd.DataFrame]:
    """Compute Betweenness Centrality, Closeness Centrality and PageRank.

    Reads the graph from *edge_list_path* and writes three CSV files into
    *output_dir*:
      - betweenness_centrality.csv
      - closeness_centrality.csv
      - pagerank.csv

    If *giant_component_only* is True (default), metrics are computed only on
    the largest connected component of the graph.

    Returns a dict mapping metric name to the corresponding DataFrame.
    """
    if not edge_list_path.exists():
        raise FileNotFoundError(f"Edge list not found: {edge_list_path}")

    edge_df = pd.read_csv(edge_list_path)
    _validate_columns(edge_df, REQUIRED_EDGE_COLUMNS, edge_list_path)

    graph = _build_graph(edge_df)

    if giant_component_only:
        graph = _extract_giant_component(graph)
        print(
            f"Giant component: {graph.number_of_nodes()} nodes, "
            f"{graph.number_of_edges()} edges"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    betweenness_df = compute_betweenness_centrality(
        graph, output_dir / "betweenness_centrality.csv"
    )
    closeness_df = compute_closeness_centrality(
        graph, output_dir / "closeness_centrality.csv"
    )
    pagerank_df = compute_pagerank(graph, output_dir / "pagerank.csv")

    return {
        "betweenness_centrality": betweenness_df,
        "closeness_centrality": closeness_df,
        "pagerank": pagerank_df,
    }
