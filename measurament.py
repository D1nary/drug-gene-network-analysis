from __future__ import annotations

import argparse
from pathlib import Path

import networkx as nx
import pandas as pd


DEFAULT_BASE_DIR = Path("results") / "similarity"
DEFAULT_NODE_METRICS = DEFAULT_BASE_DIR / "node_metrics.csv"
DEFAULT_EDGE_LIST = DEFAULT_BASE_DIR / "edge_list.csv"
DEFAULT_OUTPUT = DEFAULT_BASE_DIR / "community_measurament.csv"


REQUIRED_NODE_COLUMNS = {"drug_id", "community_id"}
REQUIRED_EDGE_COLUMNS = {"drug_1", "drug_2", "similarity"}


def _validate_columns(df: pd.DataFrame, required: set[str], file_path: Path) -> None:
    missing = required - set(df.columns)
    if missing:
        missing_str = ", ".join(sorted(missing))
        raise ValueError(f"{file_path} is missing required columns: {missing_str}")


def _build_graph(edge_df: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()
    for row in edge_df.itertuples(index=False):
        graph.add_edge(str(row.drug_1), str(row.drug_2), weight=float(row.similarity))
    return graph


def compute_community_measuraments(
    node_metrics_path: Path,
    edge_list_path: Path,
) -> pd.DataFrame:
    node_df = pd.read_csv(node_metrics_path)
    edge_df = pd.read_csv(edge_list_path)

    _validate_columns(node_df, REQUIRED_NODE_COLUMNS, node_metrics_path)
    _validate_columns(edge_df, REQUIRED_EDGE_COLUMNS, edge_list_path)

    node_df = node_df.copy()
    node_df["drug_id"] = node_df["drug_id"].astype(str)

    graph = _build_graph(edge_df)

    # Ensure isolated nodes are represented in the graph too.
    graph.add_nodes_from(node_df["drug_id"].tolist())

    rows: list[dict[str, float | int]] = []

    grouped = node_df.groupby("community_id", sort=True)
    for community_id, group in grouped:
        members = group["drug_id"].tolist()
        subgraph = graph.subgraph(members)
        size = int(subgraph.number_of_nodes())

        density = float(nx.density(subgraph)) if size > 1 else 0.0
        clustering = (
            float(nx.average_clustering(subgraph, weight="weight")) if size > 1 else 0.0
        )

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
        default=DEFAULT_NODE_METRICS,
        help=f"Path to node_metrics.csv (default: {DEFAULT_NODE_METRICS})",
    )
    parser.add_argument(
        "--edge-list",
        type=Path,
        default=DEFAULT_EDGE_LIST,
        help=f"Path to edge_list.csv (default: {DEFAULT_EDGE_LIST})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output CSV path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.node_metrics.exists():
        raise FileNotFoundError(f"Node metrics file not found: {args.node_metrics}")
    if not args.edge_list.exists():
        raise FileNotFoundError(f"Edge list file not found: {args.edge_list}")

    output_df = compute_community_measuraments(
        node_metrics_path=args.node_metrics,
        edge_list_path=args.edge_list,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(args.output, index=False)

    print(f"Saved community measuraments to {args.output}")
    print(f"Communities processed: {len(output_df)}")


if __name__ == "__main__":
    main()
