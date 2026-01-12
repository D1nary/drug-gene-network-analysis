#!/usr/bin/env python3
"""Utility script to inspect and clean the ChG-InterDecagon targets dataset."""

from __future__ import annotations

import argparse
from pathlib import Path
import networkx as nx
import pandas as pd


from network import (
    build_drug_target_network,
    build_drug_similarity_network,
    build_community_network,
)
from saves import save_network_parameters, save_community_data
from visualizzation import (
    visualize_random_drug_target_subgraph,
    visualize_similarity_subgraph,
)


PROJECT_ROOT = Path(__file__).resolve().parent

# Default location of the compressed targets dataset
DEFAULT_DATA_PATH = PROJECT_ROOT / "data" / "ChG-InterDecagon_targets.csv.gz"
RESULTS_DIR = PROJECT_ROOT / "results"
COMMUNITY_DIR = RESULTS_DIR / "community"
COMMUNITY_METRICS_DIR = COMMUNITY_DIR / "community_network_metrics"


def read_targets(path: Path) -> pd.DataFrame:
    """Load the targets CSV/GZ file and perform initial sanity checks."""

    # Load the raw targets with explicit column names
    df = pd.read_csv(path, sep=",", header=None, names=["Drug", "Gene"])

    first_drug = str(df.iloc[0]["Drug"])
    if first_drug.startswith("#"):
        print("Detected header artifact in the first row. Dropping it.")

        # Drop the artifact row and reset indices to keep downstream logic simple
        df = df.iloc[1:].reset_index(drop=True)

    # Ensure the Gene column is numeric; invalid values become NaN
    df["Gene"] = pd.to_numeric(df["Gene"], errors="coerce")
    return df


def preprocess(df: pd.DataFrame) -> pd.DataFrame:
    """Clean, normalize, and deduplicate the dataset."""

    # Remove any rows that still look like commented headers
    df = df[~df["Drug"].astype(str).str.startswith("#")]

    # Keep only rows with both Drug and Gene values
    df = df.dropna(subset=["Drug", "Gene"])

    # Normalize text spacing and types up front
    df["Drug"] = df["Drug"].astype(str).str.strip()
    df["Gene"] = pd.to_numeric(df["Gene"], errors="coerce").astype(float)

    # Avoid repeated Drug-Gene pairs
    df = df.drop_duplicates(subset=["Drug", "Gene"])

    # Strip the CID prefix and turn Drug into an integer identifier for future analysis
    df["Drug"] = df["Drug"].str.replace("CID", "", regex=False).astype(int)

    return df


def describe(df: pd.DataFrame) -> None:
    """Print a short report about the dataset."""

    # Provide an immediate preview for manual inspection
    print("Preview of the dataset:")
    print(df.head())

    print("\nDataset shape (rows, columns):", df.shape)

    print(f"Number of unique drugs: {df['Drug'].nunique()}")
    print(f"Number of unique genes: {df['Gene'].nunique()}")

    print("\nBasic info:")
    print(df.info())
    print(df.head())
    
    print(f"Total unique drugs: {df['Drug'].nunique()}")
    print(f"Total unique genes: {df['Gene'].nunique()}")
    print(f"Total interactions: {len(df)}")


def print_graph_sample(graph: nx.Graph) -> None:
    """Display one sample drug node, one gene node, and one connecting edge."""
    first_node = next(iter(graph.nodes(data=True)), None)
    sample_drug = next(
        (node for node, data in graph.nodes(data=True) if data.get("bipartite") == "drug"),
        None,
    )
    sample_gene = next(
        (node for node, data in graph.nodes(data=True) if data.get("bipartite") == "gene"),
        None,
    )
    sample_edge = next(iter(graph.edges()), None)

    if not all([sample_drug, sample_gene, sample_edge]):
        print("Graph sample unavailable: missing nodes or edges.")
        return

    print("\nGraph sample:")
    if first_node:
        node_name, attrs = first_node
        print(f"  First node: {node_name}")
        print(f"  First node attributes: {attrs}")
    print(f"  Drug node example: {sample_drug}")
    print(f"  Gene node example: {sample_gene}")
    print(f"  Edge example: {sample_edge}")
    print("\n")


def parse_args() -> argparse.Namespace:
    # Build the CLI interface for selecting alternative dataset paths
    parser = argparse.ArgumentParser(
        description="Inspect and preprocess the ChG-InterDecagon targets dataset."
    )
    parser.add_argument(
        "--data-path",
        type=Path,
        default=DEFAULT_DATA_PATH,
        help="Path to the ChG-InterDecagon_targets.csv(.gz) file "
        "(defaults to ./data/ChG-InterDecagon_targets.csv.gz)",
    )
    return parser.parse_args()


def main() -> None:
    # Parse CLI options and resolve the dataset path
    args = parse_args()
    data_path = args.data_path.expanduser().resolve()

    # Make sure the results directories exist before further processing
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "drug_gene").mkdir(parents=True, exist_ok=True)
    COMMUNITY_DIR.mkdir(parents=True, exist_ok=True)

    if not data_path.exists():
        raise FileNotFoundError(f"Dataset not found at {data_path}")

    # Sequentially process and summarize the dataset
    df = read_targets(data_path)
    df = preprocess(df)
    describe(df)

    graph = build_drug_target_network(df)
    print_graph_sample(graph)
    drug_spotlight = visualize_random_drug_target_subgraph(
        graph,
        title="Mid-degree drug spotlight",
        focus="mid",
    )
    drug_param_paths = save_network_parameters(
        drug_spotlight,
        label="mid_degree_drug_spotlight",
    )
    print(
        "Saved mid-degree spotlight parameters to",
        drug_param_paths["global"].parent,
    )

    similarity_threshold = 0.40
    similarity_graph = build_drug_similarity_network(
        df,
        similarity_threshold=similarity_threshold,
    )
    filtering_details = {
        "similarity_threshold": similarity_threshold,
        "nodes_removed": similarity_graph.graph.get("removed_drugs"),
        "edges_filtered": similarity_graph.graph.get("filtered_edges"),
        "original_node_count": similarity_graph.graph.get("original_drug_count"),
        "retained_node_count": similarity_graph.graph.get("retained_drug_count"),
        "potential_edges": similarity_graph.graph.get("potential_edges"),
    }
    if similarity_graph.number_of_nodes() == 0:
        print("Similarity graph is empty; skipping visualization.")
    else:
        snapshot_seed = 2
        similarity_snapshot = visualize_similarity_subgraph(
            similarity_graph,
            max_nodes=500,
            title="Random drug similarity snapshot",
        )
        similarity_param_paths = save_network_parameters(
            similarity_snapshot,
            label="random_similarity_snapshot",
            filtering_details=filtering_details,
        )
        print(
            "Saved random similarity snapshot parameters to",
            similarity_param_paths["global"].parent,
        )

        community_graph, membership, communities = build_community_network(
            similarity_graph,
            weight="weight",
            resolution=1.0,
            seed=42,
        )
        community_count = len(communities)
        largest = max((len(c) for c in communities), default=0)
        print(
            f"Louvain communities detected: {community_count} "
            f"(largest size {largest})"
        )

        save_community_data(community_graph, COMMUNITY_METRICS_DIR)

        visualize_similarity_subgraph(
            similarity_snapshot,
            max_nodes=similarity_snapshot.number_of_nodes(),
            title="Drug similarity snapshot by community",
            seed=snapshot_seed,
            community_membership=membership,
            output_dir=COMMUNITY_DIR,
            max_legend_items=20,
            legend_columns=1,
        )



if __name__ == "__main__":
    main()
