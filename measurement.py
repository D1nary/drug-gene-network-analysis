from __future__ import annotations

import argparse
from pathlib import Path

import networkx as nx
import pandas as pd


# ---------------------------------------------------------------------------
# Community measurement defaults
# ---------------------------------------------------------------------------
DEFAULT_SIMILARITY_BASE_DIR = Path("results") / "similarity"
DEFAULT_COMMUNITY_NODE_METRICS = DEFAULT_SIMILARITY_BASE_DIR / "node_metrics.csv"
DEFAULT_COMMUNITY_EDGE_LIST = DEFAULT_SIMILARITY_BASE_DIR / "edge_list.csv"
DEFAULT_COMMUNITY_OUTPUT = DEFAULT_SIMILARITY_BASE_DIR / "community_measurement.csv"

REQUIRED_NODE_COLUMNS: set[str] = {"drug_id", "community_id"}
REQUIRED_COMMUNITY_EDGE_COLUMNS: set[str] = {"drug_1", "drug_2", "similarity"}

# ---------------------------------------------------------------------------
# CBP measurement defaults
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CBP_INPUT_DIR = PROJECT_ROOT / "data_for_cbp"
DEFAULT_CBP_OUTPUT_DIR = PROJECT_ROOT / "results" / "cbp_measurement"
DEFAULT_CBP_EDGE_LIST = DEFAULT_CBP_INPUT_DIR / "edge_list.csv"

REQUIRED_CBP_EDGE_COLUMNS: set[str] = {"drug_1", "drug_2"}


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _validate_columns(df: pd.DataFrame, required: set[str], file_path: Path) -> None:
    missing = required - set(df.columns)
    if missing:
        missing_str = ", ".join(sorted(missing))
        raise ValueError(f"{file_path} is missing required columns: {missing_str}")


def _build_graph(edge_df: pd.DataFrame) -> nx.Graph:
    """Build a graph from an edge list DataFrame.

    The *similarity* column is used as edge weight when present.
    """
    graph = nx.Graph()
    has_weight = "similarity" in edge_df.columns
    for row in edge_df.itertuples(index=False):
        if has_weight:
            graph.add_edge(str(row.drug_1), str(row.drug_2), weight=float(row.similarity))
        else:
            graph.add_edge(str(row.drug_1), str(row.drug_2))
    return graph


def _extract_giant_component(graph: nx.Graph) -> nx.Graph:
    """Return the largest connected component of *graph* as a subgraph view."""
    largest_cc = max(nx.connected_components(graph), key=len)
    return graph.subgraph(largest_cc).copy()


# ---------------------------------------------------------------------------
# Community measurements
# ---------------------------------------------------------------------------

def compute_community_measurements(
    node_metrics_path: Path,
    edge_list_path: Path,
) -> pd.DataFrame:
    """Compute per-community metrics (density, size, clustering coefficient,
    weighted degree) from Louvain community assignments."""
    node_df = pd.read_csv(node_metrics_path)
    edge_df = pd.read_csv(edge_list_path)

    _validate_columns(node_df, REQUIRED_NODE_COLUMNS, node_metrics_path)
    _validate_columns(edge_df, REQUIRED_COMMUNITY_EDGE_COLUMNS, edge_list_path)

    node_df = node_df.copy()
    node_df["drug_id"] = node_df["drug_id"].astype(str)

    graph = _build_graph(edge_df)

    # Ensure isolated nodes are represented in the graph too.
    graph.add_nodes_from(node_df["drug_id"].tolist())

    rows: list[dict[str, float | int]] = []

    # Iterate over each community, extracting the subgraph induced by its members.
    grouped = node_df.groupby("community_id", sort=True)
    for community_id, group in grouped:
        members = group["drug_id"].tolist()
        # Induced subgraph retains only edges between nodes in this community.
        subgraph = graph.subgraph(members)
        size = int(subgraph.number_of_nodes())

        # Single-node communities have no edges, so density and clustering are 0.
        density = float(nx.density(subgraph)) if size > 1 else 0.0
        clustering = (
            float(nx.average_clustering(subgraph, weight="weight")) if size > 1 else 0.0
        )

        # Sum edge weights incident to each node, then average over the community size.
        weighted_degree_total = float(
            sum(weight for _, weight in subgraph.degree(weight="weight"))
        )
        weighted_degree = weighted_degree_total / size if size > 0 else 0.0

        rows.append(
            {
                "community_id": int(community_id),
                "size": size,
                "density": density,
                "clustering_coefficient": clustering,
                "weighted_degree": weighted_degree,
                "weighted_degree_total": weighted_degree_total,
            }
        )

    return pd.DataFrame(rows).sort_values("community_id").reset_index(drop=True)


# ---------------------------------------------------------------------------
# CBP measurements (betweenness centrality, closeness centrality, PageRank)
# ---------------------------------------------------------------------------

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


def compute_degree(
    graph: nx.Graph,
    output_path: Path,
    weight: str | None = "weight",
) -> pd.DataFrame:
    """Compute degree (and weighted degree when weights are present) for all nodes and save results to CSV."""
    degree = dict(graph.degree())
    weighted_degree = dict(graph.degree(weight=weight)) if weight else degree
    df = pd.DataFrame(
        {
            "node": list(degree.keys()),
            "degree": list(degree.values()),
            "weighted_degree": [weighted_degree[n] for n in degree],
        }
    ).sort_values("degree", ascending=False).reset_index(drop=True)

    stats = pd.DataFrame(
        [
            {"node": "mean", "degree": df["degree"].mean(), "weighted_degree": df["weighted_degree"].mean()},
            {"node": "std",  "degree": df["degree"].std(),  "weighted_degree": df["weighted_degree"].std()},
        ]
    )
    output_df = pd.concat([df, stats], ignore_index=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(output_path, index=False)
    return df


def compute_cbp_measurements(
    edge_list_path: Path = DEFAULT_CBP_EDGE_LIST,
    output_dir: Path = DEFAULT_CBP_OUTPUT_DIR,
    giant_component_only: bool = True,
) -> dict[str, pd.DataFrame]:
    """Compute Betweenness Centrality, Closeness Centrality, PageRank and Degree.

    Reads the graph from *edge_list_path* and writes four CSV files into
    *output_dir*:
      - betweenness_centrality.csv
      - closeness_centrality.csv
      - pagerank.csv
      - degree.csv

    If *giant_component_only* is True (default), metrics are computed only on
    the largest connected component of the graph.

    Returns a dict mapping metric name to the corresponding DataFrame.
    """
    if not edge_list_path.exists():
        raise FileNotFoundError(f"Edge list not found: {edge_list_path}")

    edge_df = pd.read_csv(edge_list_path)
    _validate_columns(edge_df, REQUIRED_CBP_EDGE_COLUMNS, edge_list_path)

    graph = _build_graph(edge_df)

    if giant_component_only:
        # Restrict metrics to the largest connected component to avoid
        # inflating centrality scores with disconnected fragments.
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
    degree_df = compute_degree(graph, output_dir / "degree.csv")

    return {
        "betweenness_centrality": betweenness_df,
        "closeness_centrality": closeness_df,
        "pagerank": pagerank_df,
        "degree": degree_df,
    }


# ---------------------------------------------------------------------------
# CLI (community measurements)
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-community metrics (density, size, clustering coefficient, "
            "weighted degree) from Louvain assignments."
        )
    )
    parser.add_argument(
        "--node-metrics",
        type=Path,
        default=DEFAULT_COMMUNITY_NODE_METRICS,
        help=f"Path to node_metrics.csv (default: {DEFAULT_COMMUNITY_NODE_METRICS})",
    )
    parser.add_argument(
        "--edge-list",
        type=Path,
        default=DEFAULT_COMMUNITY_EDGE_LIST,
        help=f"Path to edge_list.csv (default: {DEFAULT_COMMUNITY_EDGE_LIST})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_COMMUNITY_OUTPUT,
        help=f"Output CSV path (default: {DEFAULT_COMMUNITY_OUTPUT})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.node_metrics.exists():
        raise FileNotFoundError(f"Node metrics file not found: {args.node_metrics}")
    if not args.edge_list.exists():
        raise FileNotFoundError(f"Edge list file not found: {args.edge_list}")

    output_df = compute_community_measurements(
        node_metrics_path=args.node_metrics,
        edge_list_path=args.edge_list,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(args.output, index=False)

    print(f"Saved community measurements to {args.output}")
    print(f"Communities processed: {len(output_df)}")


if __name__ == "__main__":
    main()