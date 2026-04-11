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
    compute_bipartite_graph_stats,
    build_node2vec_embeddings,
    build_metapath2vec_pq,
    build_metapath2vec_standard,
    build_metapath2vec_pp,
    build_metapath2vec_pp_v0,
    build_cosine_similarity_network_from_node2vec_embeddings,
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
from measurament import compute_community_measuraments
from measurement import compute_cbp_measurements


PROJECT_ROOT = Path(__file__).resolve().parent

# Default location of the compressed targets dataset
DEFAULT_DATA_PATH = PROJECT_ROOT / "data" / "ChG-InterDecagon_targets.csv.gz"
RESULTS_DIR = PROJECT_ROOT / "results"
PARAMETERS_DIR = PROJECT_ROOT / "parameters"
EMBEDDING_DIR = RESULTS_DIR / "embedding"
COMMUNITY_DIR = RESULTS_DIR / "community"
COMMUNITY_METRICS_DIR = COMMUNITY_DIR / "community_network_metrics"
DAG_DIR = RESULTS_DIR / "dag"
DAG_COMMUNITY_INFO_DIR = DAG_DIR / "communities"
NODE2VEC_PARAMETERS_DIR = PARAMETERS_DIR / "node2vec"
METAPATH2VEC_PARAMETERS_DIR = PARAMETERS_DIR / "metapath2vec"
NODE2VEC_HYPERPARAMETERS_PATH = NODE2VEC_PARAMETERS_DIR / "hyperparameters.json"
METAPATH2VEC_HYPERPARAMETERS_PATH = METAPATH2VEC_PARAMETERS_DIR / "hyperparameters.json"


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


def preprocess_dataset(df: pd.DataFrame) -> pd.DataFrame:
    """Syntactically clean and normalize Drug/Gene identifiers.

    Steps:
    - remove commented/header-like rows,
    - strip spaces and empty values,
    - remove identifier prefixes (e.g. CID),
    - convert identifiers to numeric values,
    - drop incomplete and duplicated rows.
    """

    required_cols = {"Drug", "Gene"}
    if not required_cols.issubset(df.columns):
        raise ValueError(
            f"DataFrame must contain columns {required_cols}, "
            f"but has {set(df.columns)}."
        )

    tmp = df[["Drug", "Gene"]].copy()

    # Initial syntactic cleaning: remove comments/artifacts and trim spaces.
    tmp["Drug"] = tmp["Drug"].astype(str).str.strip()
    tmp["Gene"] = tmp["Gene"].astype(str).str.strip()
    tmp = tmp[~tmp["Drug"].str.startswith("#")]
    tmp = tmp.replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})
    tmp = tmp.dropna(subset=["Drug", "Gene"])

    # Identifier normalization: strip CID-like prefixes and convert to numbers.
    tmp["Drug"] = tmp["Drug"].str.replace("CID", "", regex=False).str.strip()
    tmp["Drug"] = pd.to_numeric(tmp["Drug"], errors="coerce")
    tmp["Gene"] = pd.to_numeric(tmp["Gene"], errors="coerce")
    tmp = tmp.dropna(subset=["Drug", "Gene"])

    # Canonical numeric formats for downstream graph creation.
    tmp["Drug"] = tmp["Drug"].astype(int)
    tmp["Gene"] = tmp["Gene"].astype(float)

    # Remove duplicates to avoid artificial node/edge inflation.
    tmp = tmp.drop_duplicates(subset=["Drug", "Gene"]).reset_index(drop=True)
    return tmp


def preprocess(df: pd.DataFrame) -> pd.DataFrame:
    """Backward-compatible alias for dataset preprocessing."""
    return preprocess_dataset(df)


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
    # Normalize heterogeneous CSV representations (JSON-like strings, Python lists, NaN)
    # into a clean Python list so target sets can be safely built for DAG generation.
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
    # Load per-community drug metadata generated in the community analysis step.
    community_nodes_path = community_dir / "community_nodes.csv"
    if not community_nodes_path.exists():
        print(f"Community {community_id} has no node info at {community_nodes_path}")
        return

    df = pd.read_csv(community_nodes_path)
    targets_by_node: dict[int, set[int]] = {}
    # Build a mapping: drug_id -> set(target_gene_ids), skipping malformed rows.
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

    # Create the inclusion DAG where edges represent target-set containment relations.
    dag = build_target_inclusion_dag(
        targets_by_node,
        min_set_difference=min_set_difference,
    )
    if remove_transitive_edges:
        # Optional transitive reduction keeps only informative direct dependencies.
        dag = reduce_transitive_edges(dag)
    # Persist a visual rendering before filtering by depth constraints.
    visualize_community_dag(
        dag,
        community_id=community_id,
        output_dir=community_dir / "graph",
    )

    # Filter communities by DAG depth to keep only structures within the chosen range.
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


def load_node2vec_hyperparameters(path: Path) -> dict[str, int | float | None]:
    """Load and validate the Node2Vec hyperparameter configuration from JSON."""

    if not path.exists():
        raise FileNotFoundError(f"Node2Vec hyperparameter file not found at {path}")

    with path.open("r", encoding="utf-8") as handle:
        hyperparameters = json.load(handle)

    required_keys = {
        "dimensions": int,
        "walk_length": int,
        "num_walks": int,
        "p": (int, float),
        "q": (int, float),
        "window_size": int,
        "epochs": int,
        "seed": int,
        "max_start_nodes": (int, type(None)),
        "workers": int,
        "min_count": int,
        "sg": int,
    }
    missing_keys = [key for key in required_keys if key not in hyperparameters]
    if missing_keys:
        raise ValueError(
            "Node2Vec hyperparameter file is missing required keys: "
            + ", ".join(sorted(missing_keys))
        )

    for key, expected_type in required_keys.items():
        if not isinstance(hyperparameters[key], expected_type):
            raise TypeError(
                f"Invalid type for Node2Vec hyperparameter '{key}': "
                f"expected {expected_type}, got {type(hyperparameters[key])}."
            )

    if hyperparameters["dimensions"] <= 0:
        raise ValueError("Node2Vec hyperparameter 'dimensions' must be > 0.")
    if hyperparameters["walk_length"] <= 1:
        raise ValueError("Node2Vec hyperparameter 'walk_length' must be > 1.")
    if hyperparameters["num_walks"] <= 0:
        raise ValueError("Node2Vec hyperparameter 'num_walks' must be > 0.")
    if hyperparameters["p"] <= 0 or hyperparameters["q"] <= 0:
        raise ValueError("Node2Vec hyperparameters 'p' and 'q' must be > 0.")
    if hyperparameters["window_size"] <= 0:
        raise ValueError("Node2Vec hyperparameter 'window_size' must be > 0.")
    if hyperparameters["epochs"] <= 0:
        raise ValueError("Node2Vec hyperparameter 'epochs' must be > 0.")
    if hyperparameters["workers"] <= 0:
        raise ValueError("Node2Vec hyperparameter 'workers' must be > 0.")
    if hyperparameters["min_count"] < 0:
        raise ValueError("Node2Vec hyperparameter 'min_count' must be >= 0.")
    if hyperparameters["sg"] not in {0, 1}:
        raise ValueError("Node2Vec hyperparameter 'sg' must be 0 or 1.")
    if (
        hyperparameters["max_start_nodes"] is not None
        and hyperparameters["max_start_nodes"] <= 0
    ):
        raise ValueError(
            "Node2Vec hyperparameter 'max_start_nodes' must be > 0 when provided."
        )

    return hyperparameters


def load_metapath2vec_hyperparameters(
    path: Path,
) -> dict[str, int | float | list[str]]:
    """Load and validate the metapath2vec hyperparameter configuration from JSON."""

    if not path.exists():
        raise FileNotFoundError(f"metapath2vec hyperparameter file not found at {path}")

    with path.open("r", encoding="utf-8") as handle:
        hyperparameters = json.load(handle)

    required_keys = {
        "embedding_dim": int,
        "window_size": int,
        "walk_length": int,
        "num_walks": int,
        "metapath": list,
        "negative_samples": int,
        "epochs": int,
        "learning_rate": (int, float),
        "min_count": int,
        "p": (int, float),
        "q": (int, float),
    }
    missing_keys = [key for key in required_keys if key not in hyperparameters]
    if missing_keys:
        raise ValueError(
            "metapath2vec hyperparameter file is missing required keys: "
            + ", ".join(sorted(missing_keys))
        )

    for key, expected_type in required_keys.items():
        if not isinstance(hyperparameters[key], expected_type):
            raise TypeError(
                f"Invalid type for metapath2vec hyperparameter '{key}': "
                f"expected {expected_type}, got {type(hyperparameters[key])}."
            )

    if hyperparameters["embedding_dim"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'embedding_dim' must be > 0.")
    if hyperparameters["window_size"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'window_size' must be > 0.")
    if hyperparameters["walk_length"] <= 1:
        raise ValueError("metapath2vec hyperparameter 'walk_length' must be > 1.")
    if hyperparameters["num_walks"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'num_walks' must be > 0.")
    if hyperparameters["negative_samples"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'negative_samples' must be > 0.")
    if hyperparameters["epochs"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'epochs' must be > 0.")
    if hyperparameters["learning_rate"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'learning_rate' must be > 0.")
    if hyperparameters["min_count"] < 0:
        raise ValueError("metapath2vec hyperparameter 'min_count' must be >= 0.")
    if hyperparameters["p"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'p' must be > 0.")
    if hyperparameters["q"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'q' must be > 0.")
    if len(hyperparameters["metapath"]) < 2:
        raise ValueError("metapath2vec hyperparameter 'metapath' must have length >= 2.")
    if not all(isinstance(item, str) and item.strip() for item in hyperparameters["metapath"]):
        raise TypeError(
            "metapath2vec hyperparameter 'metapath' must contain non-empty strings."
        )

    return hyperparameters


def load_metapath2vec_standard_hyperparameters(
    path: Path,
) -> dict[str, int | float | list[str]]:
    """Load and validate the standard Metapath2Vec hyperparameters from JSON.

    Reads the same file as the pq variant but only requires the parameters
    relevant to the standard algorithm (p and q are ignored).
    """

    if not path.exists():
        raise FileNotFoundError(f"metapath2vec hyperparameter file not found at {path}")

    with path.open("r", encoding="utf-8") as handle:
        hyperparameters = json.load(handle)

    required_keys = {
        "embedding_dim": int,
        "window_size": int,
        "walk_length": int,
        "num_walks": int,
        "metapath": list,
        "negative_samples": int,
        "epochs": int,
        "learning_rate": (int, float),
        "min_count": int,
    }
    missing_keys = [key for key in required_keys if key not in hyperparameters]
    if missing_keys:
        raise ValueError(
            "metapath2vec hyperparameter file is missing required keys: "
            + ", ".join(sorted(missing_keys))
        )

    for key, expected_type in required_keys.items():
        if not isinstance(hyperparameters[key], expected_type):
            raise TypeError(
                f"Invalid type for metapath2vec hyperparameter '{key}': "
                f"expected {expected_type}, got {type(hyperparameters[key])}."
            )

    if hyperparameters["embedding_dim"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'embedding_dim' must be > 0.")
    if hyperparameters["window_size"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'window_size' must be > 0.")
    if hyperparameters["walk_length"] <= 1:
        raise ValueError("metapath2vec hyperparameter 'walk_length' must be > 1.")
    if hyperparameters["num_walks"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'num_walks' must be > 0.")
    if hyperparameters["negative_samples"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'negative_samples' must be > 0.")
    if hyperparameters["epochs"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'epochs' must be > 0.")
    if hyperparameters["learning_rate"] <= 0:
        raise ValueError("metapath2vec hyperparameter 'learning_rate' must be > 0.")
    if hyperparameters["min_count"] < 0:
        raise ValueError("metapath2vec hyperparameter 'min_count' must be >= 0.")
    if len(hyperparameters["metapath"]) < 2:
        raise ValueError("metapath2vec hyperparameter 'metapath' must have length >= 2.")
    if not all(isinstance(item, str) and item.strip() for item in hyperparameters["metapath"]):
        raise TypeError(
            "metapath2vec hyperparameter 'metapath' must contain non-empty strings."
        )

    return {k: hyperparameters[k] for k in required_keys}


def parse_args() -> argparse.Namespace:
    # Build the CLI interface for selecting alternative dataset paths
    parser = argparse.ArgumentParser(
        description="Inspect and preprocess the ChG-InterDecagon targets dataset."
    )
    parser.add_argument(
        "--embedding-algorithm",
        choices=["node2vec", "metapath2vec", "metapath2vec-standard", "metapath2vec-pp", "metapath2vec-pp-v0"],
        default="metapath2vec-pp-v0",
        help=(
            "Embedding algorithm to execute when --run-embedding is enabled. "
            "Choose from: node2vec, metapath2vec (pq-biased), metapath2vec-standard "
            "(uniform transitions), metapath2vec-pp (heterogeneous skip-gram). "
            "Default: metapath2vec-standard."
        ),
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
        default=10,
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
        default=["similarity", "community", "cooccurence"],
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
        default=1,
        help="Minimum node degree required to visualize similarity network nodes.",
    )
    parser.add_argument(
        "--cosine-threshold",
        type=float,
        default=0.60,
        help="Cosine threshold in [0,1] for drug-drug similarity on node2vec embeddings.",
    )
    parser.add_argument(
        "--similarity-path-centralities",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Whether to compute asnd save similarity-node betweenness/closeness "
            "centralities. Disable to skip their calculation."
        ),
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
    parser.add_argument(
        "--run-graph-stats",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Compute and save bipartite graph statistics to results/drug_gene/graph_stats.json. "
            "Default: enabled (use --no-run-graph-stats to skip)."
        ),
    )
    parser.add_argument(
        "--run-embedding",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Run embedding generation before similarity network creation. "
            "Default: disabled."
        ),
    )
    parser.add_argument(
        "--run-similarity",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Run drug-drug similarity network creation from embeddings and save outputs "
            "under results/similarity. Default: disabled."
        ),
    )
    parser.add_argument(
        "--run-similarity-visualization",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Create similarity network visualization from results/similarity outputs. "
            "Default: enabled (use --no-run-similarity-visualization to skip)."
        ),
    )
    parser.add_argument(
        "--run-measurament",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Compute and save per-community measurament metrics from "
            "results/similarity/node_metrics.csv and edge_list.csv. "
            "Default: enabled (use --no-run-measurament to skip)."
        ),
    )
    parser.add_argument(
        "--cbp",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Compute Betweenness Centrality, Closeness Centrality and PageRank "
            "from data_for_cbp/edge_list.csv and save results to cbp_measurement/. "
            "Default: enabled (use --no-cbp to skip)."
        ),
    )
    return parser.parse_args()


def main() -> None:
    # Parse CLI options and resolve the dataset path
    args = parse_args()
    data_path = args.data_path.expanduser().resolve()  

    # Make sure the results directories exist before further processing
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    EMBEDDING_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "drug_gene").mkdir(parents=True, exist_ok=True)

    if not data_path.exists():
        raise FileNotFoundError(f"Dataset not found at {data_path}")

    # Load, clean, and summarize the source interaction table before graph construction.
    df = read_targets(data_path)
    print("Perform pre-processing\n")
    df = preprocess_dataset(df)
    describe(df)

    # Build the bipartite drug-target graph and save a representative spotlight snapshot.
    graph = build_drug_target_network(df)
    print_graph_sample(graph)
    bipartite_snapshot = visualize_random_drug_target_subgraph(
        graph,
        title="bipartite_graph",
        focus="mid",
    )
    print(
        "Bipartite graph visualization:",
        f"{bipartite_snapshot.number_of_nodes()} nodes,",
        f"{bipartite_snapshot.number_of_edges()} edges",
    )
    embedding_algorithm = args.embedding_algorithm
    if embedding_algorithm == "node2vec":
        embedding_output_path = EMBEDDING_DIR / "node2vec_embeddings.csv"
        node2vec_hyperparameters = load_node2vec_hyperparameters(
            NODE2VEC_HYPERPARAMETERS_PATH
        )
        metapath2vec_hyperparameters = None
        metapath2vec_standard_hyperparameters = None
    elif embedding_algorithm == "metapath2vec":
        embedding_output_path = EMBEDDING_DIR / "metapath2vec_pq_embeddings.csv"
        node2vec_hyperparameters = None
        metapath2vec_hyperparameters = load_metapath2vec_hyperparameters(
            METAPATH2VEC_HYPERPARAMETERS_PATH
        )
        metapath2vec_standard_hyperparameters = None
    elif embedding_algorithm == "metapath2vec-standard":
        embedding_output_path = EMBEDDING_DIR / "metapath2vec_standard_embeddings.csv"
        node2vec_hyperparameters = None
        metapath2vec_hyperparameters = None
        metapath2vec_standard_hyperparameters = load_metapath2vec_standard_hyperparameters(
            METAPATH2VEC_HYPERPARAMETERS_PATH
        )
    elif embedding_algorithm == "metapath2vec-pp":
        embedding_output_path = EMBEDDING_DIR / "metapath2vec_pp_embeddings.csv"
        node2vec_hyperparameters = None
        metapath2vec_hyperparameters = None
        metapath2vec_standard_hyperparameters = load_metapath2vec_standard_hyperparameters(
            METAPATH2VEC_HYPERPARAMETERS_PATH
        )
    elif embedding_algorithm == "metapath2vec-pp-v0":
        embedding_output_path = EMBEDDING_DIR / "metapath2vec_pp_v0_embeddings.csv"
        node2vec_hyperparameters = None
        metapath2vec_hyperparameters = None
        metapath2vec_standard_hyperparameters = load_metapath2vec_standard_hyperparameters(
            METAPATH2VEC_HYPERPARAMETERS_PATH
        )
    else:
        raise ValueError(f"Unsupported embedding algorithm: {embedding_algorithm}")

    generated_walks: list[list[str]] | None = None

    if args.run_embedding:
        print(f"\nPerform embedding algorithm: {embedding_algorithm}")
        if embedding_algorithm == "node2vec":
            embeddings_df, generated_walks = build_node2vec_embeddings(
                graph,
                dimensions=node2vec_hyperparameters["dimensions"],
                walk_length=node2vec_hyperparameters["walk_length"],
                num_walks=node2vec_hyperparameters["num_walks"],
                p=node2vec_hyperparameters["p"],
                q=node2vec_hyperparameters["q"],
                window_size=node2vec_hyperparameters["window_size"],
                epochs=node2vec_hyperparameters["epochs"],
                seed=node2vec_hyperparameters["seed"],
                max_start_nodes=node2vec_hyperparameters["max_start_nodes"],
                workers=node2vec_hyperparameters["workers"],
                min_count=node2vec_hyperparameters["min_count"],
                sg=node2vec_hyperparameters["sg"],
                output_path=embedding_output_path,
            )
            print(
                "Node2Vec embeddings saved to",
                embedding_output_path,
                (
                    f"({len(embeddings_df)} nodes, "
                    f"{node2vec_hyperparameters['dimensions']} dims)"
                ),
            )
        elif embedding_algorithm == "metapath2vec":
            embeddings_df, generated_walks = build_metapath2vec_pq(
                graph,
                embedding_dim=metapath2vec_hyperparameters["embedding_dim"],
                window_size=metapath2vec_hyperparameters["window_size"],
                walk_length=metapath2vec_hyperparameters["walk_length"],
                num_walks=metapath2vec_hyperparameters["num_walks"],
                metapath=metapath2vec_hyperparameters["metapath"],
                negative_samples=metapath2vec_hyperparameters["negative_samples"],
                epochs=metapath2vec_hyperparameters["epochs"],
                learning_rate=metapath2vec_hyperparameters["learning_rate"],
                min_count=metapath2vec_hyperparameters["min_count"],
                p=metapath2vec_hyperparameters["p"],
                q=metapath2vec_hyperparameters["q"],
                output_path=embedding_output_path,
            )
            print(
                "metapath2vec embeddings saved to",
                embedding_output_path,
                (
                    f"({len(embeddings_df)} nodes, "
                    f"{metapath2vec_hyperparameters['embedding_dim']} dims)"
                ),
            )
        elif embedding_algorithm == "metapath2vec-standard":
            embeddings_df, generated_walks = build_metapath2vec_standard(
                graph,
                embedding_dim=metapath2vec_standard_hyperparameters["embedding_dim"],
                window_size=metapath2vec_standard_hyperparameters["window_size"],
                walk_length=metapath2vec_standard_hyperparameters["walk_length"],
                num_walks=metapath2vec_standard_hyperparameters["num_walks"],
                metapath=metapath2vec_standard_hyperparameters["metapath"],
                negative_samples=metapath2vec_standard_hyperparameters["negative_samples"],
                epochs=metapath2vec_standard_hyperparameters["epochs"],
                learning_rate=metapath2vec_standard_hyperparameters["learning_rate"],
                min_count=metapath2vec_standard_hyperparameters["min_count"],
                output_path=embedding_output_path,
            )
            print(
                "metapath2vec-standard embeddings saved to",
                embedding_output_path,
                (
                    f"({len(embeddings_df)} nodes, "
                    f"{metapath2vec_standard_hyperparameters['embedding_dim']} dims)"
                ),
            )
        elif embedding_algorithm == "metapath2vec-pp":
            embeddings_df, generated_walks = build_metapath2vec_pp(
                graph,
                embedding_dim=metapath2vec_standard_hyperparameters["embedding_dim"],
                window_size=metapath2vec_standard_hyperparameters["window_size"],
                walk_length=metapath2vec_standard_hyperparameters["walk_length"],
                num_walks=metapath2vec_standard_hyperparameters["num_walks"],
                metapath=metapath2vec_standard_hyperparameters["metapath"],
                negative_samples=metapath2vec_standard_hyperparameters["negative_samples"],
                epochs=metapath2vec_standard_hyperparameters["epochs"],
                learning_rate=metapath2vec_standard_hyperparameters["learning_rate"],
                min_count=metapath2vec_standard_hyperparameters["min_count"],
                output_path=embedding_output_path,
            )
            print(
                "metapath2vec-pp embeddings saved to",
                embedding_output_path,
                (
                    f"({len(embeddings_df)} nodes, "
                    f"{metapath2vec_standard_hyperparameters['embedding_dim']} dims)"
                ),
            )
        elif embedding_algorithm == "metapath2vec-pp-v0":
            embeddings_df, generated_walks = build_metapath2vec_pp_v0(
                graph,
                embedding_dim=metapath2vec_standard_hyperparameters["embedding_dim"],
                window_size=metapath2vec_standard_hyperparameters["window_size"],
                walk_length=metapath2vec_standard_hyperparameters["walk_length"],
                num_walks=metapath2vec_standard_hyperparameters["num_walks"],
                metapath=metapath2vec_standard_hyperparameters["metapath"],
                negative_samples=metapath2vec_standard_hyperparameters["negative_samples"],
                epochs=metapath2vec_standard_hyperparameters["epochs"],
                learning_rate=metapath2vec_standard_hyperparameters["learning_rate"],
                min_count=metapath2vec_standard_hyperparameters["min_count"],
                output_path=embedding_output_path,
            )
            print(
                "metapath2vec-pp-v0 embeddings saved to",
                embedding_output_path,
                (
                    f"({len(embeddings_df)} nodes, "
                    f"{metapath2vec_standard_hyperparameters['embedding_dim']} dims)"
                ),
            )
        else:
            raise ValueError(f"Unsupported embedding algorithm: {embedding_algorithm}")
    else:
        print(
            f"Skipping {embedding_algorithm} embedding generation (--no-run-embedding set).",
            "Using existing embeddings at",
            embedding_output_path,
        )
        if not embedding_output_path.exists():
            raise FileNotFoundError(
                "Embedding generation is disabled and no embedding file was found at "
                f"{embedding_output_path}. Enable --run-embedding or provide the file."
            )

    if args.run_graph_stats:
        stats = compute_bipartite_graph_stats(
            graph,
            output_path=RESULTS_DIR / "drug_gene" / "graph_stats.json",
            all_walks=generated_walks,
        )
        print(
            "Bipartite graph stats saved to results/drug_gene/graph_stats.json",
            f"({stats['global']['n_drug']} drugs, {stats['global']['n_gene']} genes,"
            f" {stats['global']['n_edges']} edges)",
        )
    else:
        print("Skipping bipartite graph stats (--no-run-graph-stats set).")

    if args.run_similarity:
        print("Similarity netwok creation:")
        similarity_graph, similarity_global_params, similarity_global_path = (
            build_cosine_similarity_network_from_node2vec_embeddings(
                embedding_path=embedding_output_path,
                cosine_threshold=args.cosine_threshold,
                output_path=RESULTS_DIR / "similarity" / "global_parameters.json",
                include_path_centralities=args.similarity_path_centralities,
            )
        )
        print(
            "Drug-drug similarity network built:",
            f"{similarity_graph.number_of_nodes()} nodes,",
            f"{similarity_graph.number_of_edges()} edges",
        )
        print("Similarity global parameters saved to", similarity_global_path)
        print("Similarity summary:", similarity_global_params)
    else:
        print(
            "Skipping similarity network creation (--no-run-similarity set).",
            "No similarity outputs were generated in this run.",
        )

    if args.run_similarity_visualization:
        similarity_snapshot = visualize_similarity_subgraph(
            min_degree=args.similarity_min_degree,
            title="Similarity-network (embedding-based algorithm)",
            filename="similarity_network",
        )
        print(
            "Similarity network visualization:",
            f"{similarity_snapshot.number_of_nodes()} nodes,",
            f"{similarity_snapshot.number_of_edges()} edges",
        )
    else:
        print(
            "Skipping similarity network visualization "
            "(--no-run-similarity-visualization set)."
        )

    if args.run_measurament:
        node_metrics_path = RESULTS_DIR / "similarity" / "node_metrics.csv"
        edge_list_path = RESULTS_DIR / "similarity" / "edge_list.csv"
        measurament_output_path = (
            RESULTS_DIR / "similarity" / "community_measurament.csv"
        )

        if not node_metrics_path.exists():
            raise FileNotFoundError(
                "Measurament step enabled but node metrics file was not found at "
                f"{node_metrics_path}."
            )
        if not edge_list_path.exists():
            raise FileNotFoundError(
                "Measurament step enabled but edge list file was not found at "
                f"{edge_list_path}."
            )

        measurament_df = compute_community_measuraments(
            node_metrics_path=node_metrics_path,
            edge_list_path=edge_list_path,
        )
        measurament_output_path.parent.mkdir(parents=True, exist_ok=True)
        measurament_df.to_csv(measurament_output_path, index=False)
        print(
            "Community measurament metrics saved to",
            measurament_output_path,
            f"({len(measurament_df)} communities)",
        )
    else:
        print("Skipping measurament step (--no-run-measurament set).")

    if args.cbp:
        cbp_input = PROJECT_ROOT / "data_for_cbp" / "edge_list.csv"
        cbp_output = RESULTS_DIR / "cbp_measurement"

        if not cbp_input.exists():
            raise FileNotFoundError(
                "CBP step enabled but edge list was not found at "
                f"{cbp_input}. Provide the file or disable with --no-cbp."
            )

        cbp_results = compute_cbp_measurements(
            edge_list_path=cbp_input,
            output_dir=cbp_output,
        )
        print(
            "CBP measurements saved to",
            cbp_output,
            f"(betweenness: {len(cbp_results['betweenness_centrality'])} nodes,"
            f" closeness: {len(cbp_results['closeness_centrality'])} nodes,"
            f" pagerank: {len(cbp_results['pagerank'])} nodes)",
        )
    else:
        print("Skipping CBP measurements (--no-cbp set).")


if __name__ == "__main__":
    main()
