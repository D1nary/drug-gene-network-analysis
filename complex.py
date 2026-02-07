#!/usr/bin/env python3
"""Utility script to inspect and clean the ChG-InterDecagon targets dataset."""

from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path
import networkx as nx
import pandas as pd


from network import (
    build_drug_target_network,
    build_gene_cooccurrence_network,
    build_drug_similarity_network,
    build_community_network,
    compute_community_normalized_laplacian_spectra,
    build_target_inclusion_dag,
    reduce_transitive_edges,
)
from saves import (
    compute_cooccurrence_parameters,
    save_community_normalized_laplacian_spectra,
    save_community_data,
    save_community_members_by_density,
    save_community_profile_summary,
    save_community_node_info,
    save_cooccurrence_parameters_by_community,
    save_network_parameters,
    save_dag_parameters,
    compute_dag_global_parameters,
)
from visualizzation import (
    visualize_random_drug_target_subgraph,
    visualize_similarity_subgraph,
    visualize_community_dag,
)


PROJECT_ROOT = Path(__file__).resolve().parent

# Default location of the compressed targets dataset
DEFAULT_DATA_PATH = PROJECT_ROOT / "data" / "ChG-InterDecagon_targets.csv.gz"
RESULTS_DIR = PROJECT_ROOT / "results"
COMMUNITY_DIR = RESULTS_DIR / "community"
COMMUNITY_METRICS_DIR = COMMUNITY_DIR / "community_network_metrics"
DAG_DIR = RESULTS_DIR / "dag"
DAG_COMMUNITY_INFO_DIR = DAG_DIR / "communities"


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


def _parse_target_list(value: object) -> list[int] | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    if isinstance(value, list):
        return value
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        try:
            parsed = ast.literal_eval(text)
        except (ValueError, SyntaxError):
            return None
    if not isinstance(parsed, list):
        return None
    return parsed


def build_and_save_community_dag(
    community_id: int,
    community_dir: Path,
    min_set_difference: int = 3,
    min_depth: int = 1,
    max_depth: int = 3,
    remove_transitive_edges: bool = True,
) -> None:
    community_nodes_path = community_dir / "community_nodes.csv"
    if not community_nodes_path.exists():
        print(f"Community {community_id} has no node info at {community_nodes_path}")
        return

    df = pd.read_csv(community_nodes_path)
    targets_by_node: dict[int, set[int]] = {}
    for row in df.itertuples(index=False):
        drug_id = getattr(row, "drug_id", None)
        if pd.isna(drug_id):
            continue
        targets = _parse_target_list(getattr(row, "target_list", None))
        if not targets:
            continue
        try:
            node_id = int(drug_id)
        except (TypeError, ValueError):
            continue
        targets_by_node[node_id] = {int(t) for t in targets}

    if not targets_by_node:
        print(f"Community {community_id} has no targets to build a DAG.")
        return

    dag = build_target_inclusion_dag(
        targets_by_node,
        min_set_difference=min_set_difference,
    )
    if remove_transitive_edges:
        dag = reduce_transitive_edges(dag)
    visualize_community_dag(
        dag,
        community_id=community_id,
        output_dir=community_dir / "graph",
    )

    global_params = compute_dag_global_parameters(dag)
    if global_params["max_depth"] < min_depth:
        print(
            f"Community {community_id} skipped: max_depth "
            f"{global_params['max_depth']} < {min_depth}."
        )
        return
    if global_params["max_depth"] > max_depth:
        print(
            f"Community {community_id} skipped: max_depth "
            f"{global_params['max_depth']} > {max_depth}."
        )
        return

    dag_paths = save_dag_parameters(dag, community_dir)
    print(
        f"Saved DAG parameters for community {community_id} to {dag_paths['global'].parent}"
    )


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
    parser.add_argument(
        "--community-density-threshold",
        type=float,
        default=0.99,
        help=(
            "Save community member lists for communities whose density exceeds "
            "this threshold (e.g. 0.9)."
        ),
    )
    parser.add_argument(
        "--community-min-size",
        type=int,
        default=20,
        help=(
            "Minimum community size for saving member lists when density filtering "
            "is enabled. Used as the minimum size to compute per-community "
            "clustering coefficients."
        ),
    )
    parser.add_argument(
        "--networks",
        nargs="+",
        choices=["similarity", "community", "cooccurence"],
        default=["similarity", "community"],
        help=(
            "Networks to build and save. "
            "Choose from: similarity, community, cooccurence."
        ),
    )
    parser.add_argument(
        "--cooccurrence-min-drugs-per-gene",
        type=int,
        default=1,
        help="Minimum number of drugs per gene to keep genes in co-occurrence.",
    )
    parser.add_argument(
        "--cooccurrence-max-drugs-percentile",
        type=float,
        default=95.0,
        help=(
            "Percentile cutoff for max drugs per gene when filtering the "
            "co-occurrence network (e.g. 95)."
        ),
    )
    parser.add_argument(
        "--cooccurrence-community-min-size",
        type=int,
        default=15,
        help=(
            "Minimum community size required to build per-community "
            "co-occurrence networks."
        ),
    )
    parser.add_argument(
        "--similarity-min-degree",
        type=int,
        default=0,
        help="Minimum node degree required to visualize similarity network nodes.",
    )
    parser.add_argument(
        "--laplacian-community-min-size",
        type=int,
        default=4,
        help="Minimum community size required to compute Laplacian spectra.",
    )
    parser.add_argument(
        "--community-ids",
        nargs="+",
        type=int,
        default=[21],
        help=(
            "Community IDs to save under results/dag/communities "
            "(default: 21)."
        ),
    )
    return parser.parse_args()


def main() -> None:
    # Parse CLI options and resolve the dataset path
    args = parse_args()
    data_path = args.data_path.expanduser().resolve()
    requested = set(args.networks)

    # Make sure the results directories exist before further processing
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "drug_gene").mkdir(parents=True, exist_ok=True)
    COMMUNITY_DIR.mkdir(parents=True, exist_ok=True)
    DAG_COMMUNITY_INFO_DIR.mkdir(parents=True, exist_ok=True)

    for community_id in sorted(set(args.community_ids)):
        community_dir = DAG_COMMUNITY_INFO_DIR / str(community_id)
        community_dir.mkdir(parents=True, exist_ok=True)
        (community_dir / "graph").mkdir(parents=True, exist_ok=True)

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
    print(
        "Mid-degree drug spotlight graph:",
        f"{drug_spotlight.number_of_nodes()} nodes,",
        f"{drug_spotlight.number_of_edges()} edges",
    )
    drug_param_paths = save_network_parameters(
        drug_spotlight,
        label="mid_degree_drug_spotlight",
        output_root=RESULTS_DIR / "drug_gene",
    )
    print(
        "Saved mid-degree spotlight parameters to",
        drug_param_paths["global"].parent,
    )

    similarity_graph = None
    similarity_snapshot = None
    community_graph = None
    membership = {}
    communities = []
    if {"similarity", "community", "cooccurence"} & requested:
        print("SIMILARITY NETWORK CREATION")
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
                seed=snapshot_seed,
                min_degree=args.similarity_min_degree,
            )
            similarity_param_paths = save_network_parameters(
                similarity_graph,
                label="similarity_network",
                output_root=RESULTS_DIR / "similarity",
                filtering_details=filtering_details,
            )
            print(
                "Saved similarity network parameters to",
                similarity_param_paths["global"].parent,
            )

    if {"community", "cooccurence"} & requested:
        if similarity_graph is None or similarity_graph.number_of_nodes() == 0:
            print("Community/co-occurrence requested but similarity graph is empty.")
        else:
            print("COMMUNITY NETWORK CREATION")
            community_graph, membership, communities = build_community_network(
                similarity_graph,
                weight="weight",
                resolution=1.0,
                seed=42,
                min_clustering_size=args.community_min_size,
            )

    if "community" in requested:
        if community_graph is None or community_graph.number_of_nodes() == 0:
            print("Community requested but no communities were generated.")
            return
        community_count = len(communities)
        largest = max((len(c) for c in communities), default=0)
        print(
            f"Louvain communities detected: {community_count} "
            f"(largest size {largest})"
        )

        save_community_data(community_graph, COMMUNITY_METRICS_DIR)
        total_targets = int(df["Gene"].nunique())
        for community_id in sorted(set(args.community_ids)):
            if community_id < 0 or community_id >= community_count:
                print(
                    f"Community id {community_id} out of range "
                    f"(available: 0..{community_count - 1})."
                )
                continue
            members = communities[community_id]
            output_path = save_community_node_info(
                community_id,
                members,
                similarity_graph,
                total_targets,
                DAG_COMMUNITY_INFO_DIR / str(community_id),
            )
            print(f"Saved community {community_id} node info to {output_path}")
        laplacian_spectra = compute_community_normalized_laplacian_spectra(
            similarity_graph,
            communities,
            min_size=args.laplacian_community_min_size,
        )
        laplacian_path = save_community_normalized_laplacian_spectra(
            laplacian_spectra,
            RESULTS_DIR / "similarity" / "similarity_network",
        )
        print(f"Saved normalized Laplacian spectra to {laplacian_path}")
        if args.community_density_threshold is not None:
            drug_targets = {
                node: data.get("targets", [])
                for node, data in similarity_graph.nodes(data=True)
                if data.get("bipartite") == "drug"
            }
            members_path = save_community_members_by_density(
                community_graph,
                args.community_density_threshold,
                args.community_min_size,
                COMMUNITY_METRICS_DIR,
            )
            print(
                "Saved community member lists (density >= "
                f"{args.community_density_threshold}, size >= "
                f"{args.community_min_size}) to {members_path}"
            )
            summary_path = save_community_profile_summary(
                community_graph,
                drug_targets,
                args.community_density_threshold,
                args.community_min_size,
                COMMUNITY_METRICS_DIR,
            )
            print(
                "Saved community profile summary (density >= "
                f"{args.community_density_threshold}, size >= "
                f"{args.community_min_size}) to {summary_path}"
            )

        if similarity_snapshot is not None:
            visualize_similarity_subgraph(
                similarity_snapshot,
                max_nodes=similarity_snapshot.number_of_nodes(),
                title="Drug similarity snapshot by community",
                seed=snapshot_seed,
                min_degree=-1,
                community_membership=membership,
                output_dir=COMMUNITY_DIR,
                max_legend_items=20,
                legend_columns=1,
            )

    if "cooccurence" in requested:
        if not communities:
            print("Co-occurrence requested but no communities were generated.")
        else:
            print("CO-OCCURRENCE NETWORK CREATION")
            community_parameters = {}
            for idx, community in enumerate(communities):
                if len(community) < args.cooccurrence_community_min_size:
                    continue
                label = f"Community_{idx}"
                drug_ids = [
                    int(node.replace("Drug_", ""))
                    for node in community
                    if str(node).startswith("Drug_")
                ]
                if not drug_ids:
                    continue
                community_df = df[df["Drug"].isin(drug_ids)]
                cooccurrence_graph = build_gene_cooccurrence_network(
                    community_df,
                    min_drugs_per_gene=args.cooccurrence_min_drugs_per_gene,
                    max_drugs_per_gene_percentile=args.cooccurrence_max_drugs_percentile,
                )
                params = compute_cooccurrence_parameters(
                    cooccurrence_graph,
                    weight_ge_threshold=1,
                )
                params["community_size"] = len(community)
                community_parameters[label] = params

            cooccurrence_path = save_cooccurrence_parameters_by_community(
                community_parameters
            )
            print(f"Saved co-occurrence parameters to {cooccurrence_path}")

    community_count = len(communities)
    for community_id in sorted(set(args.community_ids)):
        if community_count and (community_id < 0 or community_id >= community_count):
            print(
                f"Community id {community_id} out of range "
                f"(available: 0..{community_count - 1})."
            )
            continue
        build_and_save_community_dag(
            community_id,
            DAG_COMMUNITY_INFO_DIR / str(community_id),
            min_set_difference=3,
            min_depth=1,
            max_depth=10,
            remove_transitive_edges=True,
        )



if __name__ == "__main__":
    main()
