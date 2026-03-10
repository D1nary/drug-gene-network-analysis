#!/usr/bin/env python3
"""Utilities to build a Drug–Target bipartite network from 
the ChG-InterDecagon targets dataset."""

from __future__ import annotations

import heapq
from itertools import combinations
from pathlib import Path
import pandas as pd
import networkx as nx
import numpy as np
import scipy.sparse as sp
from saves import (
    compute_similarity_communities,
    save_similarity_edge_list,
    save_similarity_global_parameters,
    save_similarity_node_metrics,
)


def build_drug_target_network(df: pd.DataFrame) -> nx.Graph:
    """Build a bipartite Drug–Target Network from a preprocessed dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame already preprocessed (returned by `preprocess` in complex.py),
        which must contain at least the columns:
        - 'Drug' : numeric identifier for the drug (int)
        - 'Gene' : numeric identifier for the gene (float/int)

    Returns
    -------
    nx.Graph
        A bipartite NetworkX graph where:
        - nodes with attribute ``bipartite='drug'`` represent drugs,
        - nodes with attribute ``bipartite='gene'`` represent genes,
        - an edge connects a drug node and a gene node whenever the dataset
          contains at least one drug–gene interaction.
    """

    # Basic validation that the required columns are present
    required_cols = {"Drug", "Gene"}
    if not required_cols.issubset(df.columns):
        raise ValueError(
            f"DataFrame must contain columns {required_cols}, "
            f"but has {set(df.columns)}."
        )

    # Work on a copy that only keeps the columns of interest
    tmp = df[["Drug", "Gene"]].copy()

    # Ensure IDs are integers (Gene is stored as float in preprocess)
    tmp["Drug"] = tmp["Drug"].astype(int)
    tmp["Gene"] = tmp["Gene"].astype(int)

    # Drop duplicates again as a safeguard, even if preprocess already did it
    tmp = tmp.drop_duplicates(subset=["Drug", "Gene"])

    # Build the bipartite graph
    G = nx.Graph()

    # Sets of nodes
    unique_drugs = tmp["Drug"].unique()
    unique_genes = tmp["Gene"].unique()

    # Add drug nodes with attributes
    G.add_nodes_from(
        (
            f"Drug_{drug_id}",
            {"bipartite": "drug", "original_id": int(drug_id)},
        )
        for drug_id in unique_drugs
    )

    # Add gene nodes with attributes
    G.add_nodes_from(
        (
            f"Gene_{gene_id}",
            {"bipartite": "gene", "original_id": int( gene_id)},
        )
        for gene_id in unique_genes
    )

    # Connect each Drug–Gene interaction with an edge
    for row in tmp.itertuples(index=False):
        drug_id = int(row.Drug)
        gene_id = int(row.Gene)
        G.add_edge(f"Drug_{drug_id}", f"Gene_{gene_id}")

    return G


def build_node2vec_embeddings(
    graph: nx.Graph,
    dimensions: int = 128,
    walk_length: int = 10,
    num_walks: int = 3,
    p: float = 1.0,
    q: float = 0.5,
    window_size: int = 3,
    seed: int = 42,
    max_start_nodes: int | None = None,
    max_cooccurrence_pairs: int = 2_000_000,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Compute Node2Vec-style embeddings and save them to CSV.

    The output schema is:
    ``node_id,node_type,emb_1,...,emb_d`` where ``d = dimensions``.
    """

    if dimensions <= 0:
        raise ValueError("dimensions must be > 0.")
    if walk_length <= 1:
        raise ValueError("walk_length must be > 1.")
    if num_walks <= 0:
        raise ValueError("num_walks must be > 0.")
    if window_size <= 0:
        raise ValueError("window_size must be > 0.")
    if p <= 0 or q <= 0:
        raise ValueError("p and q must be > 0.")
    if max_start_nodes is not None and max_start_nodes <= 0:
        raise ValueError("max_start_nodes must be > 0 when provided.")
    if max_cooccurrence_pairs <= 0:
        raise ValueError("max_cooccurrence_pairs must be > 0.")

    nodes = sorted(graph.nodes())
    if not nodes:
        raise ValueError("Cannot compute embeddings on an empty graph.")

    rng = np.random.default_rng(seed)
    node_to_idx = {node: idx for idx, node in enumerate(nodes)}
    n_nodes = len(nodes)

    # Precompute neighborhood lists once to reduce overhead during walks.
    neighbors = {node: list(graph.neighbors(node)) for node in nodes}

    def _next_step(prev_node: str | None, curr_node: str) -> str | None:
        curr_neighbors = neighbors.get(curr_node, [])
        if not curr_neighbors:
            return None
        if prev_node is None:
            probs = np.full(len(curr_neighbors), 1.0 / len(curr_neighbors))
            return curr_neighbors[int(rng.choice(len(curr_neighbors), p=probs))]

        weights = []
        for candidate in curr_neighbors:
            edge_weight = float(graph[curr_node][candidate].get("weight", 1.0))
            if candidate == prev_node:
                bias = 1.0 / p
            elif graph.has_edge(candidate, prev_node):
                bias = 1.0
            else:
                bias = 1.0 / q
            weights.append(edge_weight * bias)

        probs = np.asarray(weights, dtype=float)
        probs_sum = probs.sum()
        if probs_sum <= 0:
            probs = np.full(len(curr_neighbors), 1.0 / len(curr_neighbors))
        else:
            probs = probs / probs_sum
        return curr_neighbors[int(rng.choice(len(curr_neighbors), p=probs))]

    start_node_pool = nodes
    if max_start_nodes is not None and len(nodes) > max_start_nodes:
        # Limit starting points to keep runtime bounded on large graphs.
        start_node_pool = sorted(
            rng.choice(nodes, size=max_start_nodes, replace=False).tolist()
        )

    walks: list[list[str]] = []
    for _ in range(num_walks):
        start_nodes = start_node_pool.copy()
        rng.shuffle(start_nodes)
        for start in start_nodes:
            walk = [start]
            prev = None
            curr = start
            for _ in range(walk_length - 1):
                nxt = _next_step(prev, curr)
                if nxt is None:
                    break
                walk.append(nxt)
                prev, curr = curr, nxt
            walks.append(walk)

    # Build a sparse co-occurrence matrix from random walks (skip-gram context window).
    cooc_counts: dict[tuple[int, int], float] = {}
    for walk in walks:
        walk_len = len(walk)
        for i, center in enumerate(walk):
            center_idx = node_to_idx[center]
            left = max(0, i - window_size)
            right = min(walk_len, i + window_size + 1)
            for j in range(left, right):
                if i == j:
                    continue
                context_idx = node_to_idx[walk[j]]
                key = (center_idx, context_idx)
                cooc_counts[key] = cooc_counts.get(key, 0.0) + 1.0

    if not cooc_counts:
        cooc = sp.csr_matrix((n_nodes, n_nodes), dtype=float)
    else:
        if len(cooc_counts) > max_cooccurrence_pairs:
            # Keep only the strongest co-occurrences to bound memory/time.
            top_pairs = heapq.nlargest(
                max_cooccurrence_pairs,
                cooc_counts.items(),
                key=lambda item: item[1],
            )
            cooc_counts = dict(top_pairs)
        rows, cols, data = zip(
            *((i, j, v) for (i, j), v in cooc_counts.items())
        )
        cooc = sp.csr_matrix((data, (rows, cols)), shape=(n_nodes, n_nodes), dtype=float)

    total = float(cooc.sum())
    if total <= 0:
        # Fallback for fully disconnected graphs.
        embedding_matrix = np.zeros((n_nodes, dimensions), dtype=float)
    else:
        # Normalize co-occurrences row-wise and project to a dense latent space.
        row_sums = np.asarray(cooc.sum(axis=1)).ravel()
        inv_row_sums = np.zeros_like(row_sums, dtype=float)
        nonzero_rows = row_sums > 0
        inv_row_sums[nonzero_rows] = 1.0 / row_sums[nonzero_rows]
        normalized_cooc = sp.diags(inv_row_sums) @ cooc

        random_projection = rng.normal(
            loc=0.0,
            scale=1.0 / max(1, dimensions),
            size=(n_nodes, dimensions),
        )
        embedding_matrix = normalized_cooc @ random_projection
        embedding_matrix = np.asarray(embedding_matrix, dtype=float)

        # L2 normalize each embedding vector for numerical stability.
        norms = np.linalg.norm(embedding_matrix, axis=1, keepdims=True)
        nonzero_norms = norms.squeeze() > 0
        embedding_matrix[nonzero_norms] = (
            embedding_matrix[nonzero_norms] / norms[nonzero_norms]
        )

    rows = []
    for node in nodes:
        node_idx = node_to_idx[node]
        attrs = graph.nodes[node]
        node_type = attrs.get("bipartite", "unknown")
        row = {"node_id": node, "node_type": node_type}
        for dim_idx in range(dimensions):
            row[f"emb_{dim_idx + 1}"] = float(embedding_matrix[node_idx, dim_idx])
        rows.append(row)

    embeddings_df = pd.DataFrame(rows)

    if output_path is None:
        output_path = Path(__file__).resolve().parent / "results" / "embedding" / "node_embeddings.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    embeddings_df.to_csv(output_path, index=False)
    return embeddings_df


def build_drug_similarity_network(
    df: pd.DataFrame,
    similarity_threshold: float = 0.30,
) -> nx.Graph:
    """Create a weighted drug–drug similarity network using Jaccard similarity.

    Each drug is represented by the set of its targets. Two drugs are linked if
    the Jaccard similarity between their target sets is greater than or equal to
    ``similarity_threshold``. Drugs connected to only one target are removed
    before computing similarities so that low-information nodes do not skew the
    network.

    Parameters
    ----------
    df : pandas.DataFrame
        Preprocessed dataframe containing at least ``Drug`` and ``Gene`` columns.
    similarity_threshold : float, optional
        Absolute similarity threshold applied to Jaccard similarities (default 0.30).

    Returns
    -------
    nx.Graph
        Weighted drug–drug network where edges carry a ``weight`` attribute
        equal to the Jaccard similarity of the incident nodes.
    """

    required_cols = {"Drug", "Gene"}
    if not required_cols.issubset(df.columns):
        raise ValueError(
            f"DataFrame must contain columns {required_cols}, "
            f"but has {set(df.columns)}."
        )

    # Normalize input schema before set-based similarity computations.
    tmp = df[["Drug", "Gene"]].copy()
    tmp["Drug"] = tmp["Drug"].astype(int)
    tmp["Gene"] = tmp["Gene"].astype(int)
    tmp = tmp.drop_duplicates(subset=["Drug", "Gene"])
    original_drug_count = tmp["Drug"].nunique()

    # Map each drug to its unique targets and filter out mono-target drugs.
    drug_targets = tmp.groupby("Drug")["Gene"].apply(lambda genes: set(genes)).to_dict()
    drug_targets = {
        drug: targets for drug, targets in drug_targets.items() if len(targets) > 1
    }
    retained_drug_count = len(drug_targets)
    removed_drugs = max(original_drug_count - retained_drug_count, 0)

    G = nx.Graph()
    if not drug_targets:
        # Keep empty-graph metadata so downstream reporting still has context.
        G.graph.update(
            {
                "similarity_threshold": similarity_threshold,
                "original_drug_count": original_drug_count,
                "retained_drug_count": 0,
                "removed_drugs": removed_drugs,
                "potential_edges": 0,
                "filtered_edges": 0,
            }
        )
        return G

    for drug_id, targets in drug_targets.items():
        # Persist each drug's target list as node metadata for later community profiling.
        G.add_node(
            f"Drug_{drug_id}",
            bipartite="drug",
            original_id=int(drug_id),
            targets=sorted(targets),
        )

    # Compare each drug pair once and keep only edges above the threshold.
    for drug_a, drug_b in combinations(drug_targets.keys(), 2):
        set_a = drug_targets[drug_a]
        set_b = drug_targets[drug_b]
        intersection = len(set_a & set_b)
        union = len(set_a | set_b)
        if union == 0:
            continue

        similarity = float(intersection / union)
        if similarity >= similarity_threshold:
            G.add_edge(
                f"Drug_{drug_a}",
                f"Drug_{drug_b}",
                weight=similarity,
            )

    potential_edges = retained_drug_count * (retained_drug_count - 1) // 2
    filtered_edges = max(potential_edges - G.number_of_edges(), 0)
    G.graph.update(
        {
            "similarity_threshold": similarity_threshold,
            "original_drug_count": original_drug_count,
            "retained_drug_count": retained_drug_count,
            "removed_drugs": removed_drugs,
            "potential_edges": potential_edges,
            "filtered_edges": filtered_edges,
            "similarity_metric": "jaccard",
        }
    )

    return G


def build_cosine_similarity_network_from_node2vec_embeddings(
    embedding_path: Path | str | None = None,
    cosine_threshold: float = 0.75,
    output_path: Path | str | None = None,
    block_size: int = 256,
) -> tuple[nx.Graph, dict[str, float | int], Path]:
    """Build a cosine-similarity network from node2vec embeddings CSV.

    The input CSV is expected to contain at least:
    - ``node_id`` column
    - one or more embedding columns named ``emb_*``
    """

    if cosine_threshold < -1.0 or cosine_threshold > 1.0:
        raise ValueError("cosine_threshold must be in [-1, 1].")
    if block_size <= 0:
        raise ValueError("block_size must be > 0.")

    if embedding_path is None:
        embedding_path = (
            Path(__file__).resolve().parent
            / "results"
            / "embedding"
            / "node_embeddings.csv"
        )
    embedding_path = Path(embedding_path)
    if not embedding_path.exists():
        raise FileNotFoundError(f"Embedding file not found at {embedding_path}")

    embeddings_df = pd.read_csv(embedding_path)
    if "node_id" not in embeddings_df.columns:
        raise ValueError("Embedding file must contain a 'node_id' column.")

    emb_cols = [col for col in embeddings_df.columns if col.startswith("emb_")]
    if not emb_cols:
        raise ValueError("Embedding file must contain embedding columns named 'emb_*'.")

    node_ids = embeddings_df["node_id"].astype(str).tolist()
    node_attrs = embeddings_df.drop(columns=emb_cols).set_index("node_id")

    emb_matrix = embeddings_df[emb_cols].to_numpy(dtype=float)
    if emb_matrix.size == 0:
        raise ValueError("Embedding matrix is empty.")

    # Robust L2 normalization ensures cosine == dot product.
    norms = np.linalg.norm(emb_matrix, axis=1, keepdims=True)
    nonzero = norms.squeeze() > 0
    emb_matrix[nonzero] = emb_matrix[nonzero] / norms[nonzero]

    graph = nx.Graph()
    for node in node_ids:
        attrs = node_attrs.loc[node].to_dict() if node in node_attrs.index else {}
        cleaned_attrs = {
            str(key): (value.item() if isinstance(value, np.generic) else value)
            for key, value in attrs.items()
            if not pd.isna(value)
        }
        graph.add_node(node, **cleaned_attrs)

    n_nodes = len(node_ids)
    for start in range(0, n_nodes, block_size):
        end = min(start + block_size, n_nodes)
        similarities = emb_matrix[start:end] @ emb_matrix.T
        for local_idx, source_idx in enumerate(range(start, end)):
            row = similarities[local_idx]
            row[: source_idx + 1] = -np.inf
            target_indices = np.flatnonzero(row >= cosine_threshold)
            for target_idx in target_indices:
                graph.add_edge(
                    node_ids[source_idx],
                    node_ids[int(target_idx)],
                    weight=float(row[target_idx]),
                )

    if output_path is None:
        output_path = (
            Path(__file__).resolve().parent
            / "results"
            / "similarity"
            / "global_parameters.json"
        )
    communities = compute_similarity_communities(
        graph,
        weight_attr="weight",
        seed=42,
    )
    global_params, saved_output_path = save_similarity_global_parameters(
        graph,
        output_path=output_path,
        weight_attr="weight",
        communities=communities,
        community_seed=42,
    )
    _, node_metrics_output_path = save_similarity_node_metrics(
        graph,
        output_path=Path(__file__).resolve().parent
        / "results"
        / "similarity"
        / "node_metrics.csv",
        weight_attr="weight",
        community_seed=42,
        communities=communities,
    )
    _, edge_list_output_path = save_similarity_edge_list(
        graph,
        output_path=Path(__file__).resolve().parent
        / "results"
        / "similarity"
        / "edge_list.csv",
        weight_attr="weight",
    )

    graph.graph.update(
        {
            "similarity_metric": "cosine",
            "cosine_threshold": cosine_threshold,
            "embedding_path": str(embedding_path),
            "global_parameters_path": str(saved_output_path),
            "node_metrics_path": str(node_metrics_output_path),
            "edge_list_path": str(edge_list_output_path),
        }
    )

    return graph, global_params, saved_output_path


def build_gene_cooccurrence_network(
    df: pd.DataFrame,
    min_drugs_per_gene: int = 3,
    max_drugs_per_gene_percentile: float = 95.0,
) -> nx.Graph:
    """Create a weighted gene co-occurrence network via shared drugs."""

    required_cols = {"Drug", "Gene"}
    if not required_cols.issubset(df.columns):
        raise ValueError(
            f"DataFrame must contain columns {required_cols}, "
            f"but has {set(df.columns)}."
        )

    # Clean identifiers and remove duplicated drug-gene interactions.
    tmp = df[["Drug", "Gene"]].copy()
    tmp["Drug"] = tmp["Drug"].astype(int)
    tmp["Gene"] = tmp["Gene"].astype(int)
    tmp = tmp.drop_duplicates(subset=["Drug", "Gene"])

    # Keep informative genes by frequency: remove both rare and overly generic genes.
    gene_drug_counts = tmp.groupby("Gene")["Drug"].nunique()
    original_gene_count = int(gene_drug_counts.shape[0])
    if gene_drug_counts.empty:
        empty_graph = nx.Graph()
        empty_graph.graph.update(
            {
                "metric": "shared_drugs_count",
                "min_drugs_per_gene": min_drugs_per_gene,
                "max_drugs_per_gene_percentile": max_drugs_per_gene_percentile,
                "max_drugs_per_gene_value": 0.0,
                "cooccurrence_weight_threshold": weight_threshold,
                "original_gene_count": 0,
                "retained_gene_count": 0,
                "removed_genes": 0,
                "filtered_edges": 0,
            }
        )
        return empty_graph

    percentile_value = float(
        gene_drug_counts.quantile(max_drugs_per_gene_percentile / 100.0)
    )
    informative_genes = gene_drug_counts[
        (gene_drug_counts >= min_drugs_per_gene)
        & (gene_drug_counts <= percentile_value)
    ].index
    retained_gene_count = int(len(informative_genes))
    removed_genes = max(original_gene_count - retained_gene_count, 0)

    filtered_tmp = tmp[tmp["Gene"].isin(informative_genes)]
    gene_ids = sorted(set(filtered_tmp["Gene"].astype(int).tolist()))
    cooccurrence_graph = nx.Graph()
    cooccurrence_graph.add_nodes_from(f"Gene_{gene_id}" for gene_id in gene_ids)

    # Count how many drugs induce each gene pair co-occurrence.
    shared_drug_counts: dict[tuple[int, int], int] = {}
    for _, group in filtered_tmp.groupby("Drug"):
        genes = sorted(set(group["Gene"].astype(int).tolist()))
        if len(genes) < 2:
            continue
        for gene_a, gene_b in combinations(genes, 2):
            key = (gene_a, gene_b)
            shared_drug_counts[key] = shared_drug_counts.get(key, 0) + 1

    for (gene_a, gene_b), weight in shared_drug_counts.items():
        cooccurrence_graph.add_edge(
            f"Gene_{gene_a}",
            f"Gene_{gene_b}",
            weight=float(weight),
        )

    # Save filtering metadata in graph attributes for reproducibility.
    cooccurrence_graph.graph.update(
        {
            "metric": "shared_drugs_count",
            "min_drugs_per_gene": min_drugs_per_gene,
            "max_drugs_per_gene_percentile": max_drugs_per_gene_percentile,
            "max_drugs_per_gene_value": percentile_value,
            "original_gene_count": original_gene_count,
            "retained_gene_count": retained_gene_count,
            "removed_genes": removed_genes,
            "filtered_edges": 0,
        }
    )

    return cooccurrence_graph


def build_community_network(
    graph: nx.Graph,
    weight: str = "weight",
    resolution: float = 1.0,
    seed: int | None = 42,
    min_clustering_size: int | None = None,
) -> tuple[nx.Graph, dict[str, str], list[set[str]]]:
    """Aggregate a graph into a community-level network via the Louvain method.

    Parameters
    ----------
    graph : nx.Graph
        Input graph whose nodes will be clustered in communities.
    weight : str, optional
        Edge attribute containing weights used for modularity optimization.
        Defaults to ``"weight"`` (unweighted if attribute is missing).
    resolution : float, optional
        Resolution parameter passed to Louvain; larger values yield more
        communities. Defaults to 1.0.
    seed : int | None, optional
        Random seed forwarded to NetworkX for reproducibility.
    min_clustering_size : int | None, optional
        Minimum community size required to compute internal clustering.
        Communities smaller than this get a 0.0 value. Defaults to None
        (compute for all sizes).

    Returns
    -------
    tuple
        ``(community_graph, membership, communities)`` where:
        - community_graph is a graph whose nodes represent detected communities
          and edges carry the summed weight between communities;
        - membership maps each original node to the community node label;
        - communities is the list of sets returned by Louvain.
    """

    if graph.number_of_nodes() == 0:
        return nx.Graph(), {}, []

    # Detect communities in the similarity graph using Louvain modularity optimization.
    communities = nx.algorithms.community.louvain_communities(
        graph,
        weight=weight,
        resolution=resolution,
        seed=seed,
    )

    membership: dict[str, str] = {}
    internal_weights: dict[str, float] = {}
    internal_edges: dict[str, int] = {}
    community_labels: list[str] = []

    for idx, community in enumerate(communities):
        label = f"Community_{idx}"
        community_labels.append(label)
        internal_weights[label] = 0.0
        internal_edges[label] = 0
        # Track node -> community assignment for later coloring/export steps.
        for node in community:
            membership[node] = label

    inter_weights: dict[tuple[str, str], float] = {}
    # Aggregate original edges into intra-community and inter-community totals.
    for u, v, data in graph.edges(data=True):
        edge_weight = float(data.get(weight, 1.0))
        comm_u = membership[u]
        comm_v = membership[v]

        if comm_u == comm_v:
            internal_weights[comm_u] += edge_weight
            internal_edges[comm_u] += 1
            continue

        key = tuple(sorted((comm_u, comm_v)))
        inter_weights[key] = inter_weights.get(key, 0.0) + edge_weight

    community_graph = nx.Graph()
    for label, community in zip(community_labels, communities):
        subgraph = graph.subgraph(community)
        # Optionally suppress clustering values for very small communities.
        if min_clustering_size is not None and len(community) < min_clustering_size:
            internal_clustering = 0.0
        else:
            internal_clustering = nx.average_clustering(subgraph, weight=weight)
        community_graph.add_node(
            label,
            size=len(community),
            members=sorted(community),
            internal_weight=internal_weights.get(label, 0.0),
            internal_edge_count=internal_edges.get(label, 0),
            internal_clustering_coefficient=internal_clustering,
        )

    for (comm_u, comm_v), edge_weight in inter_weights.items():
        community_graph.add_edge(
            comm_u,
            comm_v,
            weight=edge_weight,
        )

    modularity = nx.algorithms.community.modularity(
        graph,
        communities,
        weight=weight,
    )
    community_graph.graph["method"] = "louvain"
    community_graph.graph["resolution"] = resolution
    community_graph.graph["modularity"] = modularity

    return community_graph, membership, communities


def build_target_inclusion_dag(
    targets_by_node: dict[int | str, set[int]],
    min_set_difference: int = 3,
) -> nx.DiGraph:
    """Build a DAG using target-set inclusion (A -> B if T(B) ⊂ T(A))."""

    dag = nx.DiGraph()
    if not targets_by_node:
        return dag

    # Store target sets on nodes for interpretability of exported DAG artifacts.
    for node, targets in targets_by_node.items():
        dag.add_node(node, targets=sorted(targets))

    # Sort by target-set size so candidate supersets can be scanned efficiently.
    items = sorted(
        ((node, targets) for node, targets in targets_by_node.items()),
        key=lambda item: len(item[1]),
    )

    for i, (node_small, targets_small) in enumerate(items):
        size_small = len(targets_small)
        for node_large, targets_large in items[i + 1 :]:
            size_large = len(targets_large)
            size_diff = size_large - size_small
            # Early stop: list is size-ordered, so larger differences only increase.
            if size_diff > min_set_difference:
                break
            if size_diff <= 0:
                continue
            if targets_small.issubset(targets_large):
                dag.add_edge(node_large, node_small)

    return dag


def reduce_transitive_edges(dag: nx.DiGraph) -> nx.DiGraph:
    """Return the transitive reduction of a DAG, preserving node attributes."""

    if dag.number_of_edges() == 0:
        return dag.copy()

    # Remove redundant paths while preserving original node-level annotations.
    reduced = nx.algorithms.dag.transitive_reduction(dag)
    for node, data in dag.nodes(data=True):
        if node in reduced:
            reduced.nodes[node].update(data)
    return reduced


def compute_community_normalized_laplacian_spectra(
    graph: nx.Graph,
    communities: list[set[str]],
    weight_attr: str = "weight",
    min_size: int = 1,
) -> dict[str, dict[str, object]]:
    """Compute normalized Laplacian matrices and eigen spectra per community."""

    spectra: dict[str, dict[str, object]] = {}
    for idx, community in enumerate(communities):
        label = f"Community_{idx}"
        if not community or len(community) < min_size:
            continue
        subgraph = graph.subgraph(community)
        nodes = sorted(subgraph.nodes())
        if not nodes:
            continue
        # Export both matrix and eigenvalues to support downstream spectral analysis.
        laplacian = nx.normalized_laplacian_matrix(
            subgraph,
            nodelist=nodes,
            weight=weight_attr,
        ).toarray()
        eigenvalues = np.linalg.eigvalsh(laplacian)
        spectra[label] = {
            "node_count": len(nodes),
            "nodes": nodes,
            "normalized_laplacian": laplacian.astype(float).tolist(),
            "eigenvalues": eigenvalues.astype(float).tolist(),
        }

    return spectra
