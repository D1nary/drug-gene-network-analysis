#!/usr/bin/env python3
"""Utilities to build a Drug-Target bipartite network from
the ChG-InterDecagon targets dataset."""

from __future__ import annotations

from itertools import combinations
import json
from pathlib import Path
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite as nx_bipartite
import numpy as np
try:
    from gensim.models import Word2Vec
except ImportError:
    Word2Vec = None
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


def compute_bipartite_graph_stats(
    graph: nx.Graph,
    output_path: Path,
    all_walks: list[list[str]] | None = None,
    rare_walk_threshold: int = 5,
) -> dict:
    """Compute and save statistics of the bipartite drug-target graph.

    Saves a JSON file at *output_path* and returns the same dict.
    """
    drug_nodes = [n for n, d in graph.nodes(data=True) if d.get("bipartite") == "drug"]
    gene_nodes = [n for n, d in graph.nodes(data=True) if d.get("bipartite") == "gene"]

    n_drug = len(drug_nodes)
    n_gene = len(gene_nodes)
    n_edges = graph.number_of_edges()
    density = n_edges / (n_drug * n_gene) if n_drug and n_gene else 0.0

    drug_degrees = np.array([graph.degree(n) for n in drug_nodes], dtype=float)

    drug_proj = nx_bipartite.weighted_projected_graph(graph, drug_nodes)
    edge_weights = [d["weight"] for _, _, d in drug_proj.edges(data=True)]

    k = float(np.mean(drug_degrees)) if len(drug_degrees) > 0 else 0.0

    stats: dict = {
        "global": {
            "n_drug": n_drug,
            "n_gene": n_gene,
            "n_edges": n_edges,
            "density": round(density, 6),
        },
        "drug_degree_distribution": {
            "min": float(drug_degrees.min()) if len(drug_degrees) else 0.0,
            "max": float(drug_degrees.max()) if len(drug_degrees) else 0.0,
            "mean": round(float(np.mean(drug_degrees)), 4) if len(drug_degrees) else 0.0,
            "median": round(float(np.median(drug_degrees)), 4) if len(drug_degrees) else 0.0,
            "std": round(float(np.std(drug_degrees)), 4) if len(drug_degrees) else 0.0,
            "p10": round(float(np.percentile(drug_degrees, 10)), 4) if len(drug_degrees) else 0.0,
            "p25": round(float(np.percentile(drug_degrees, 25)), 4) if len(drug_degrees) else 0.0,
            "p75": round(float(np.percentile(drug_degrees, 75)), 4) if len(drug_degrees) else 0.0,
            "p90": round(float(np.percentile(drug_degrees, 90)), 4) if len(drug_degrees) else 0.0,
            "p95": round(float(np.percentile(drug_degrees, 95)), 4) if len(drug_degrees) else 0.0,
        },
        "gene_overlap_between_drugs": {
            "min": float(min(edge_weights)) if edge_weights else 0.0,
            "max": float(max(edge_weights)) if edge_weights else 0.0,
            "mean": round(float(np.mean(edge_weights)), 4) if edge_weights else 0.0,
            "median": round(float(np.median(edge_weights)), 4) if edge_weights else 0.0,
        },
        "walk_reachability_estimate": {
            "mean_drug_degree_k": round(k, 4),
            "drugs_reachable_1_hop": round(k, 4),
            "drugs_reachable_2_hop": int(min(k ** 2, n_drug)),
        },
    }

    if all_walks is not None:
        from collections import Counter
        drug_node_set = {str(n) for n in drug_nodes}
        walk_counts = Counter(
            node for walk in all_walks for node in walk if node in drug_node_set
        )
        n_undersampled = sum(
            1 for d in drug_nodes if walk_counts[str(d)] < rare_walk_threshold
        )
        coverage_ratio = round(1.0 - n_undersampled / n_drug, 6) if n_drug else 0.0
        stats["rare_drug_walk_coverage"] = {
            "threshold": rare_walk_threshold,
            "n_drug_total": n_drug,
            "n_drug_undersampled": n_undersampled,
            "coverage_ratio": coverage_ratio,
        }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(stats, fh, indent=2)

    return stats


def build_node2vec_embeddings(
    graph: nx.Graph,
    dimensions: int,
    walk_length: int,
    num_walks: int,
    p: float,
    q: float,
    window_size: int,
    epochs: int,
    seed: int,
    max_start_nodes: int | None,
    workers: int,
    min_count: int,
    sg: int,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Compute classic Node2Vec embeddings and save them to CSV.

    The output schema is:
    ``node_id,node_type,emb_1,...,emb_d`` where ``d = dimensions``.
    """

    if Word2Vec is None:
        raise ImportError(
            "gensim is required to train classic Node2Vec embeddings. "
            "Install it with `pip install gensim`."
        )
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
    if epochs <= 0:
        raise ValueError("epochs must be > 0.")
    if workers <= 0:
        raise ValueError("workers must be > 0.")
    if min_count < 0:
        raise ValueError("min_count must be >= 0.")
    if sg not in {0, 1}:
        raise ValueError("sg must be either 0 (CBOW) or 1 (Skip-gram).")
    if max_start_nodes is not None and max_start_nodes <= 0:
        raise ValueError("max_start_nodes must be > 0 when provided.")

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

    model = Word2Vec(
        sentences=walks,
        vector_size=dimensions,
        window=window_size,
        min_count=min_count,
        sg=sg,
        workers=workers,
        epochs=epochs,
        seed=seed,
    )

    embedding_matrix = np.zeros((n_nodes, dimensions), dtype=float)
    for node in nodes:
        if node in model.wv:
            embedding_matrix[node_to_idx[node]] = model.wv[node]

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
    return embeddings_df, walks


def build_metapath2vec_pq(
    graph: nx.Graph,
    embedding_dim: int,
    window_size: int,
    walk_length: int,
    num_walks: int,
    metapath: list[str],
    negative_samples: int,
    epochs: int,
    learning_rate: float,
    min_count: int,
    p: float = 1.0,
    q: float = 1.0,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Compute metapath2vec embeddings on a bipartite graph and save them to CSV.

    The output schema is:
    ``node_id,node_type,emb_1,...,emb_d`` where ``d = embedding_dim``.
    """

    if Word2Vec is None:
        raise ImportError(
            "gensim is required to train metapath2vec embeddings. "
            "Install it with `pip install gensim`."
        )
    if embedding_dim <= 0:
        raise ValueError("embedding_dim must be > 0.")
    if window_size <= 0:
        raise ValueError("window_size must be > 0.")
    if walk_length <= 1:
        raise ValueError("walk_length must be > 1.")
    if num_walks <= 0:
        raise ValueError("num_walks must be > 0.")
    if negative_samples <= 0:
        raise ValueError("negative_samples must be > 0.")
    if epochs <= 0:
        raise ValueError("epochs must be > 0.")
    if learning_rate <= 0:
        raise ValueError("learning_rate must be > 0.")
    if min_count < 0:
        raise ValueError("min_count must be >= 0.")
    if p <= 0:
        raise ValueError("p must be > 0.")
    if q <= 0:
        raise ValueError("q must be > 0.")
    if len(metapath) < 2:
        raise ValueError("metapath must contain at least 2 node types.")

    nodes = sorted(graph.nodes())
    if not nodes:
        raise ValueError("Cannot compute embeddings on an empty graph.")

    node_types = {
        node: str(graph.nodes[node].get("bipartite", "unknown")).strip().lower()
        for node in nodes
    }
    normalized_metapath = [str(node_type).strip().lower() for node_type in metapath]
    if any(not node_type for node_type in normalized_metapath):
        raise ValueError("metapath cannot contain empty node types.")

    unique_node_types = set(node_types.values())
    unknown_types = sorted(set(normalized_metapath) - unique_node_types)
    if unknown_types:
        raise ValueError(
            "metapath contains node types that are not present in the graph: "
            + ", ".join(unknown_types)
        )

    metapath_cycle = normalized_metapath[:-1]
    if not metapath_cycle:
        metapath_cycle = normalized_metapath
    cycle_length = len(metapath_cycle)

    type_to_offsets: dict[str, list[int]] = {}
    for idx, node_type in enumerate(metapath_cycle):
        type_to_offsets.setdefault(node_type, []).append(idx)

    missing_types = sorted(unique_node_types - set(type_to_offsets))
    if missing_types:
        raise ValueError(
            "metapath does not define traversal positions for graph node types: "
            + ", ".join(missing_types)
        )

    neighbors_by_type: dict[str, dict[str, list[str]]] = {}
    for node in nodes:
        typed_neighbors: dict[str, list[str]] = {}
        for neighbor in graph.neighbors(node):
            neighbor_type = node_types[neighbor]
            typed_neighbors.setdefault(neighbor_type, []).append(neighbor)
        neighbors_by_type[node] = typed_neighbors

    rng = np.random.default_rng(42)
    node_to_idx = {node: idx for idx, node in enumerate(nodes)}

    def _generate_walk(start_node: str, start_offset: int) -> list[str]:
        walk = [start_node]
        current_node = start_node
        current_offset = start_offset
        prev_node: str | None = None

        for _ in range(walk_length - 1):
            next_offset = (current_offset + 1) % cycle_length
            expected_type = metapath_cycle[next_offset]
            candidate_neighbors = neighbors_by_type[current_node].get(expected_type, [])
            if not candidate_neighbors:
                break

            if prev_node is None:
                next_node = str(rng.choice(candidate_neighbors))
            else:
                prev_typed_neighbors = set(
                    neighbors_by_type[prev_node].get(expected_type, [])
                )
                weights = []
                for candidate in candidate_neighbors:
                    if candidate == prev_node:
                        weights.append(1.0 / p)
                    elif candidate in prev_typed_neighbors:
                        weights.append(1.0)
                    else:
                        weights.append(1.0 / q)
                probs = np.asarray(weights, dtype=float)
                probs = probs / probs.sum()
                next_node = str(candidate_neighbors[int(rng.choice(len(candidate_neighbors), p=probs))])

            walk.append(next_node)
            prev_node = current_node
            current_node = next_node
            current_offset = next_offset

        return walk

    walks: list[list[str]] = []
    for _ in range(num_walks):
        shuffled_nodes = nodes.copy()
        rng.shuffle(shuffled_nodes)
        for node in shuffled_nodes:
            node_type = node_types[node]
            start_offset = type_to_offsets[node_type][0]
            walks.append(_generate_walk(node, start_offset))

    model = Word2Vec(
        sentences=walks,
        vector_size=embedding_dim,
        window=window_size,
        min_count=min_count,
        sg=1,
        negative=negative_samples,
        workers=1,
        epochs=epochs,
        alpha=learning_rate,
        min_alpha=learning_rate,
        seed=42,
    )

    embedding_matrix = np.zeros((len(nodes), embedding_dim), dtype=float)
    for node in nodes:
        if node in model.wv:
            embedding_matrix[node_to_idx[node]] = model.wv[node]

    norms = np.linalg.norm(embedding_matrix, axis=1, keepdims=True)
    nonzero_norms = norms.squeeze() > 0
    embedding_matrix[nonzero_norms] = (
        embedding_matrix[nonzero_norms] / norms[nonzero_norms]
    )

    rows = []
    for node in nodes:
        node_idx = node_to_idx[node]
        row = {"node_id": node, "node_type": node_types[node]}
        for dim_idx in range(embedding_dim):
            row[f"emb_{dim_idx + 1}"] = float(embedding_matrix[node_idx, dim_idx])
        rows.append(row)

    embeddings_df = pd.DataFrame(rows)

    if output_path is None:
        output_path = (
            Path(__file__).resolve().parent
            / "results"
            / "embedding"
            / "metapath2vec_pq_embeddings.csv"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    embeddings_df.to_csv(output_path, index=False)
    return embeddings_df, walks


def build_metapath2vec_standard(
    graph: nx.Graph,
    embedding_dim: int,
    window_size: int,
    walk_length: int,
    num_walks: int,
    metapath: list[str],
    negative_samples: int,
    epochs: int,
    learning_rate: float,
    min_count: int,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Compute standard Metapath2Vec embeddings on a bipartite graph and save them to CSV.

    Unlike the MP2Vec-pq variant, walks use uniform transition probabilities: at
    each step the next node is chosen uniformly at random among the neighbors that
    match the expected metapath type.  There is no p/q bias and no prev_node
    tracking.

    The output schema is:
    ``node_id,node_type,emb_1,...,emb_d`` where ``d = embedding_dim``.
    """

    if Word2Vec is None:
        raise ImportError(
            "gensim is required to train metapath2vec embeddings. "
            "Install it with `pip install gensim`."
        )
    if embedding_dim <= 0:
        raise ValueError("embedding_dim must be > 0.")
    if window_size <= 0:
        raise ValueError("window_size must be > 0.")
    if walk_length <= 1:
        raise ValueError("walk_length must be > 1.")
    if num_walks <= 0:
        raise ValueError("num_walks must be > 0.")
    if negative_samples <= 0:
        raise ValueError("negative_samples must be > 0.")
    if epochs <= 0:
        raise ValueError("epochs must be > 0.")
    if learning_rate <= 0:
        raise ValueError("learning_rate must be > 0.")
    if min_count < 0:
        raise ValueError("min_count must be >= 0.")
    if len(metapath) < 2:
        raise ValueError("metapath must contain at least 2 node types.")

    nodes = sorted(graph.nodes())
    if not nodes:
        raise ValueError("Cannot compute embeddings on an empty graph.")

    node_types = {
        node: str(graph.nodes[node].get("bipartite", "unknown")).strip().lower()
        for node in nodes
    }
    normalized_metapath = [str(node_type).strip().lower() for node_type in metapath]
    if any(not node_type for node_type in normalized_metapath):
        raise ValueError("metapath cannot contain empty node types.")

    unique_node_types = set(node_types.values())
    unknown_types = sorted(set(normalized_metapath) - unique_node_types)
    if unknown_types:
        raise ValueError(
            "metapath contains node types that are not present in the graph: "
            + ", ".join(unknown_types)
        )

    metapath_cycle = normalized_metapath[:-1]
    if not metapath_cycle:
        metapath_cycle = normalized_metapath
    cycle_length = len(metapath_cycle)

    type_to_offsets: dict[str, list[int]] = {}
    for idx, node_type in enumerate(metapath_cycle):
        type_to_offsets.setdefault(node_type, []).append(idx)

    missing_types = sorted(unique_node_types - set(type_to_offsets))
    if missing_types:
        raise ValueError(
            "metapath does not define traversal positions for graph node types: "
            + ", ".join(missing_types)
        )

    neighbors_by_type: dict[str, dict[str, list[str]]] = {}
    for node in nodes:
        typed_neighbors: dict[str, list[str]] = {}
        for neighbor in graph.neighbors(node):
            neighbor_type = node_types[neighbor]
            typed_neighbors.setdefault(neighbor_type, []).append(neighbor)
        neighbors_by_type[node] = typed_neighbors

    rng = np.random.default_rng(42)
    node_to_idx = {node: idx for idx, node in enumerate(nodes)}

    def _generate_walk(start_node: str, start_offset: int) -> list[str]:
        walk = [start_node]
        current_node = start_node
        current_offset = start_offset

        for _ in range(walk_length - 1):
            next_offset = (current_offset + 1) % cycle_length
            expected_type = metapath_cycle[next_offset]
            candidate_neighbors = neighbors_by_type[current_node].get(expected_type, [])
            if not candidate_neighbors:
                break

            next_node = str(rng.choice(candidate_neighbors))
            walk.append(next_node)
            current_node = next_node
            current_offset = next_offset

        return walk

    walks: list[list[str]] = []
    for _ in range(num_walks):
        shuffled_nodes = nodes.copy()
        rng.shuffle(shuffled_nodes)
        for node in shuffled_nodes:
            node_type = node_types[node]
            start_offset = type_to_offsets[node_type][0]
            walks.append(_generate_walk(node, start_offset))

    model = Word2Vec(
        sentences=walks,
        vector_size=embedding_dim,
        window=window_size,
        min_count=min_count,
        sg=1,
        negative=negative_samples,
        workers=1,
        epochs=epochs,
        alpha=learning_rate,
        min_alpha=learning_rate,
        seed=42,
    )

    embedding_matrix = np.zeros((len(nodes), embedding_dim), dtype=float)
    for node in nodes:
        if node in model.wv:
            embedding_matrix[node_to_idx[node]] = model.wv[node]

    norms = np.linalg.norm(embedding_matrix, axis=1, keepdims=True)
    nonzero_norms = norms.squeeze() > 0
    embedding_matrix[nonzero_norms] = (
        embedding_matrix[nonzero_norms] / norms[nonzero_norms]
    )

    rows = []
    for node in nodes:
        node_idx = node_to_idx[node]
        row = {"node_id": node, "node_type": node_types[node]}
        for dim_idx in range(embedding_dim):
            row[f"emb_{dim_idx + 1}"] = float(embedding_matrix[node_idx, dim_idx])
        rows.append(row)

    embeddings_df = pd.DataFrame(rows)

    if output_path is None:
        output_path = (
            Path(__file__).resolve().parent
            / "results"
            / "embedding"
            / "metapath2vec_standard_embeddings.csv"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    embeddings_df.to_csv(output_path, index=False)
    return embeddings_df, walks


def build_metapath2vec_pp(
    graph: nx.Graph,
    embedding_dim: int,
    window_size: int,
    walk_length: int,
    num_walks: int,
    metapath: list[str],
    negative_samples: int,
    epochs: int,
    learning_rate: float,
    min_count: int,
    output_path: Path | None = None,
) -> tuple[pd.DataFrame, list[list[str]]]:
    """Compute Metapath2Vec++ embeddings (Dong, Chawla & Swami, KDD '17).

    Walk generation is identical to :func:`build_metapath2vec_standard` — uniform
    metapath-guided random walks.  The difference is in the training objective:
    metapath2vec++ uses *heterogeneous negative sampling*, where for every
    (center, context) pair the negative nodes are drawn exclusively from V_t —
    the subset of nodes that share the same type t as the context node —
    rather than from the full node set V.

    Formally the conditional probability becomes::

        p(c_t | v; θ) = exp(X_c_t · X_v) / Σ_{u_t ∈ V_t} exp(X_u_t · X_v)

    Because gensim's Word2Vec only supports global negative sampling, the
    skip-gram update is implemented directly with NumPy.

    Output schema: ``node_id,node_type,emb_1,...,emb_d``  (d = embedding_dim).
    """
    from collections import Counter

    # ------------------------------------------------------------------ #
    # Validation                                                           #
    # ------------------------------------------------------------------ #
    if embedding_dim <= 0:
        raise ValueError("embedding_dim must be > 0.")
    if window_size <= 0:
        raise ValueError("window_size must be > 0.")
    if walk_length <= 1:
        raise ValueError("walk_length must be > 1.")
    if num_walks <= 0:
        raise ValueError("num_walks must be > 0.")
    if negative_samples <= 0:
        raise ValueError("negative_samples must be > 0.")
    if epochs <= 0:
        raise ValueError("epochs must be > 0.")
    if learning_rate <= 0:
        raise ValueError("learning_rate must be > 0.")
    if min_count < 0:
        raise ValueError("min_count must be >= 0.")
    if len(metapath) < 2:
        raise ValueError("metapath must contain at least 2 node types.")

    nodes = sorted(graph.nodes())
    if not nodes:
        raise ValueError("Cannot compute embeddings on an empty graph.")

    # ------------------------------------------------------------------ #
    # Graph structures                                                     #
    # ------------------------------------------------------------------ #
    node_types: dict[str, str] = {
        node: str(graph.nodes[node].get("bipartite", "unknown")).strip().lower()
        for node in nodes
    }
    normalized_metapath = [str(t).strip().lower() for t in metapath]
    if any(not t for t in normalized_metapath):
        raise ValueError("metapath cannot contain empty node types.")

    unique_node_types = set(node_types.values())
    unknown_types = sorted(set(normalized_metapath) - unique_node_types)
    if unknown_types:
        raise ValueError(
            "metapath contains node types not present in the graph: "
            + ", ".join(unknown_types)
        )

    # The metapath must be in cyclic notation: the first and last type must be
    # the same (e.g. ['drug','target','drug']), as per the paper's definition.
    if normalized_metapath[0] != normalized_metapath[-1]:
        raise ValueError(
            "metapath must be in cyclic notation: the last type must equal the "
            f"first (got '{normalized_metapath[0]}' != '{normalized_metapath[-1]}')."
        )

    # Drop the repeated last element so the cycle length equals the period.
    metapath_cycle = normalized_metapath[:-1] or normalized_metapath
    cycle_length = len(metapath_cycle)

    type_to_offsets: dict[str, list[int]] = {}
    for idx, t in enumerate(metapath_cycle):
        type_to_offsets.setdefault(t, []).append(idx)

    # Pre-bucket each node's neighbours by their type for O(1) look-up.
    neighbors_by_type: dict[str, dict[str, list[str]]] = {}
    for node in nodes:
        buckets: dict[str, list[str]] = {}
        for nb in graph.neighbors(node):
            buckets.setdefault(node_types[nb], []).append(nb)
        neighbors_by_type[node] = buckets

    rng = np.random.default_rng(42)
    node_to_idx: dict[str, int] = {node: i for i, node in enumerate(nodes)}
    n_nodes = len(nodes)

    # ------------------------------------------------------------------ #
    # Per-type node index (key addition for metapath2vec++)               #
    # Maps each type label → sorted list of node ids belonging to it.    #
    # Built once here; used both for validation and for negative sampling.#
    # ------------------------------------------------------------------ #
    type_to_nodes: dict[str, list[str]] = {}
    for node in nodes:
        type_to_nodes.setdefault(node_types[node], []).append(node)

    # ------------------------------------------------------------------ #
    # Meta-path random walks (uniform, identical to standard variant)     #
    # ------------------------------------------------------------------ #
    def _walk(start_node: str, start_offset: int) -> list[str]:
        path = [start_node]
        cur_node, cur_offset = start_node, start_offset
        for _ in range(walk_length - 1):
            nxt_offset = (cur_offset + 1) % cycle_length
            candidates = neighbors_by_type[cur_node].get(metapath_cycle[nxt_offset], [])
            if not candidates:
                break
            cur_node = str(rng.choice(candidates))
            path.append(cur_node)
            cur_offset = nxt_offset
        return path

    walks: list[list[str]] = []
    for _ in range(num_walks):
        shuffled = nodes.copy()
        rng.shuffle(shuffled)
        for node in shuffled:
            if node_types[node] not in type_to_offsets:
                continue  # node type not covered by metapath — not a valid walk start
            walks.append(_walk(node, type_to_offsets[node_types[node]][0]))

    # ------------------------------------------------------------------ #
    # min_count filtering                                                  #
    # ------------------------------------------------------------------ #
    node_freq: Counter = Counter(n for walk in walks for n in walk)
    effective_min = max(min_count, 1)
    valid_nodes: set[str] = {n for n, c in node_freq.items() if c >= effective_min}
    if not valid_nodes:
        valid_nodes = set(nodes)

    # ------------------------------------------------------------------ #
    # Per-type negative-sampling tables  (distribution ∝ freq^0.75)      #
    # Each type gets its own pool and probability vector so that negatives #
    # for a context node of type t are drawn solely from V_t.             #
    # ------------------------------------------------------------------ #
    type_neg_pool: dict[str, np.ndarray] = {}   # type → array of global node indices
    type_neg_prob: dict[str, np.ndarray] = {}   # type → sampling probabilities

    for t, t_nodes in type_to_nodes.items():
        valid_t = [n for n in t_nodes if n in valid_nodes]
        if not valid_t:
            continue
        freqs = np.array([node_freq.get(n, 1) for n in valid_t], dtype=float)
        probs = freqs ** 0.75
        probs /= probs.sum()
        type_neg_pool[t] = np.array([node_to_idx[n] for n in valid_t], dtype=np.intp)
        type_neg_prob[t] = probs

    # ------------------------------------------------------------------ #
    # Embedding matrices                                                   #
    # W     — input  (center) embeddings                                  #
    # W_ctx — output (context) embeddings                                 #
    # ------------------------------------------------------------------ #
    scale = 1.0 / embedding_dim
    W = rng.uniform(-scale, scale, (n_nodes, embedding_dim))
    W_ctx = np.zeros((n_nodes, embedding_dim), dtype=float)

    def _sigmoid(x: np.ndarray) -> np.ndarray:
        return 1.0 / (1.0 + np.exp(-np.clip(x, -10.0, 10.0)))

    # ------------------------------------------------------------------ #
    # Heterogeneous skip-gram with type-constrained negative sampling     #
    #                                                                     #
    # The per-pair Python loop is kept (online SGD semantics preserved)   #
    # but the three per-pair allocation bottlenecks are eliminated:       #
    #   • labels buffer pre-allocated once (not np.zeros per pair).       #
    #   • all_ctx_buf pre-allocated once (not np.concatenate per pair).   #
    #   • W_ctx update uses broadcasting (not np.outer per pair).         #
    # rng.choice is batched per type per center: all m_t contexts of the  #
    # same type get their negatives in a single (m_t, k) call.            #
    # ------------------------------------------------------------------ #

    # Pre-allocate fixed-size buffers reused across every (center, context) pair.
    _labels = np.zeros(1 + negative_samples)
    _labels[0] = 1.0                                   # label for the positive context
    _all_ctx = np.empty(1 + negative_samples, dtype=np.intp)

    for epoch in range(epochs):
        epoch_rng = np.random.default_rng(42 + epoch)
        for walk in walks:
            walk_len = len(walk)
            for center_pos in range(walk_len):
                center_node = walk[center_pos]
                if center_node not in valid_nodes:
                    continue
                center_idx = node_to_idx[center_node]

                # --- Collect context nodes and pre-sample all negatives ---
                # Group context positions by type so that all m_t contexts
                # sharing the same type get their (m_t × k) negatives in a
                # single rng.choice call instead of m_t separate calls.
                ctx_start = max(0, center_pos - window_size)
                ctx_end   = min(walk_len, center_pos + window_size + 1)

                type_ctx: dict[str, list[tuple[int, int]]] = {}  # type → [(ctx_idx, pos)]
                for p in range(ctx_start, ctx_end):
                    if p == center_pos:
                        continue
                    cn = walk[p]
                    if cn not in valid_nodes:
                        continue
                    t = node_types[cn]
                    type_ctx.setdefault(t, []).append((node_to_idx[cn], p))

                # Build per-context negative index arrays in one pass per type
                ctx_neg: list[tuple[int, np.ndarray]] = []  # (ctx_idx, neg_indices[k])
                for t, ctx_list in type_ctx.items():
                    pool  = type_neg_pool.get(t)
                    probs = type_neg_prob.get(t)
                    if pool is None or len(pool) < 2:
                        continue
                    n_t = len(ctx_list)
                    # Single call: (n_t, k) negatives for all contexts of type t
                    neg_local = epoch_rng.choice(
                        len(pool), size=(n_t, negative_samples), p=probs, replace=True
                    )
                    neg_pool_idx = pool[neg_local]          # (n_t, k)
                    for row, (ctx_idx, _) in enumerate(ctx_list):
                        ctx_neg.append((ctx_idx, neg_pool_idx[row]))

                # --- Online SGD: one update per (center, context) pair ---
                # Exact same update order as the original loop; W is read
                # fresh each time so online semantics are fully preserved.
                for ctx_idx, neg_indices in ctx_neg:
                    _all_ctx[0]  = ctx_idx
                    _all_ctx[1:] = neg_indices

                    center_vec = W[center_idx]              # (d,)  — re-read each pair
                    ctx_vecs   = W_ctx[_all_ctx]            # (1+k, d)
                    sig        = _sigmoid(ctx_vecs @ center_vec)       # (1+k,)

                    err        = sig - _labels              # (1+k,)
                    grad_ctr   = err @ ctx_vecs             # (d,)

                    # Broadcasting: (-lr * err)[:, None] * center_vec → (1+k, d)
                    # avoids np.outer and reuses the pre-allocated _all_ctx index.
                    np.add.at(
                        W_ctx, _all_ctx,
                        (-learning_rate * err[:, np.newaxis]) * center_vec,
                    )
                    W[center_idx] -= learning_rate * grad_ctr

    # ------------------------------------------------------------------ #
    # L2-normalise input embeddings                                       #
    # ------------------------------------------------------------------ #
    norms = np.linalg.norm(W, axis=1, keepdims=True)
    nonzero = norms.squeeze() > 0
    W[nonzero] = W[nonzero] / norms[nonzero]

    # ------------------------------------------------------------------ #
    # Build output DataFrame                                               #
    # ------------------------------------------------------------------ #
    rows = []
    for node in nodes:
        idx = node_to_idx[node]
        row: dict = {"node_id": node, "node_type": node_types[node]}
        for dim_idx in range(embedding_dim):
            row[f"emb_{dim_idx + 1}"] = float(W[idx, dim_idx])
        rows.append(row)

    embeddings_df = pd.DataFrame(rows)

    if output_path is None:
        output_path = (
            Path(__file__).resolve().parent
            / "results"
            / "embedding"
            / "metapath2vec_pp_embeddings.csv"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    embeddings_df.to_csv(output_path, index=False)
    return embeddings_df, walks


def build_metapath2vec_pp_v0(
    graph: nx.Graph,
    embedding_dim: int,
    window_size: int,
    walk_length: int,
    num_walks: int,
    metapath: list[str],
    negative_samples: int,
    epochs: int,
    learning_rate: float,
    min_count: int,
    output_path: Path | None = None,
) -> tuple[pd.DataFrame, list[list[str]]]:
    """Unoptimised reference implementation of Metapath2Vec++ (Dong et al., KDD '17).

    Identical to :func:`build_metapath2vec_pp` in every algorithmic respect.
    The only difference is the training loop: this version uses the naive
    per-pair implementation without any of the three bottleneck optimisations:

    * **B1** — no pre-allocation of ``labels`` / ``all_ctx_idx`` buffers.
    * **B2** — ``rng.choice`` is called once per (center, context) pair
      instead of once per type per center.
    * **B3** — ``np.concatenate``, ``np.zeros``, and ``np.outer`` are
      called inside the inner loop, allocating new arrays on every pair.

    Intended as a correctness baseline and benchmark reference.
    Output saved to ``metapath2vec_pp_v0_embeddings.csv``.
    """
    from collections import Counter

    # ------------------------------------------------------------------ #
    # Validation                                                           #
    # ------------------------------------------------------------------ #
    if embedding_dim <= 0:
        raise ValueError("embedding_dim must be > 0.")
    if window_size <= 0:
        raise ValueError("window_size must be > 0.")
    if walk_length <= 1:
        raise ValueError("walk_length must be > 1.")
    if num_walks <= 0:
        raise ValueError("num_walks must be > 0.")
    if negative_samples <= 0:
        raise ValueError("negative_samples must be > 0.")
    if epochs <= 0:
        raise ValueError("epochs must be > 0.")
    if learning_rate <= 0:
        raise ValueError("learning_rate must be > 0.")
    if min_count < 0:
        raise ValueError("min_count must be >= 0.")
    if len(metapath) < 2:
        raise ValueError("metapath must contain at least 2 node types.")

    nodes = sorted(graph.nodes())
    if not nodes:
        raise ValueError("Cannot compute embeddings on an empty graph.")

    # ------------------------------------------------------------------ #
    # Graph structures                                                     #
    # ------------------------------------------------------------------ #
    node_types: dict[str, str] = {
        node: str(graph.nodes[node].get("bipartite", "unknown")).strip().lower()
        for node in nodes
    }
    normalized_metapath = [str(t).strip().lower() for t in metapath]
    if any(not t for t in normalized_metapath):
        raise ValueError("metapath cannot contain empty node types.")

    unique_node_types = set(node_types.values())
    unknown_types = sorted(set(normalized_metapath) - unique_node_types)
    if unknown_types:
        raise ValueError(
            "metapath contains node types not present in the graph: "
            + ", ".join(unknown_types)
        )

    if normalized_metapath[0] != normalized_metapath[-1]:
        raise ValueError(
            "metapath must be in cyclic notation: the last type must equal the "
            f"first (got '{normalized_metapath[0]}' != '{normalized_metapath[-1]}')."
        )

    metapath_cycle = normalized_metapath[:-1] or normalized_metapath
    cycle_length = len(metapath_cycle)

    type_to_offsets: dict[str, list[int]] = {}
    for idx, t in enumerate(metapath_cycle):
        type_to_offsets.setdefault(t, []).append(idx)

    neighbors_by_type: dict[str, dict[str, list[str]]] = {}
    for node in nodes:
        buckets: dict[str, list[str]] = {}
        for nb in graph.neighbors(node):
            buckets.setdefault(node_types[nb], []).append(nb)
        neighbors_by_type[node] = buckets

    rng = np.random.default_rng(42)
    node_to_idx: dict[str, int] = {node: i for i, node in enumerate(nodes)}
    n_nodes = len(nodes)

    # ------------------------------------------------------------------ #
    # Per-type node index                                                  #
    # ------------------------------------------------------------------ #
    type_to_nodes: dict[str, list[str]] = {}
    for node in nodes:
        type_to_nodes.setdefault(node_types[node], []).append(node)

    # ------------------------------------------------------------------ #
    # Meta-path random walks                                               #
    # ------------------------------------------------------------------ #
    def _walk(start_node: str, start_offset: int) -> list[str]:
        path = [start_node]
        cur_node, cur_offset = start_node, start_offset
        for _ in range(walk_length - 1):
            nxt_offset = (cur_offset + 1) % cycle_length
            candidates = neighbors_by_type[cur_node].get(metapath_cycle[nxt_offset], [])
            if not candidates:
                break
            cur_node = str(rng.choice(candidates))
            path.append(cur_node)
            cur_offset = nxt_offset
        return path

    walks: list[list[str]] = []
    for _ in range(num_walks):
        shuffled = nodes.copy()
        rng.shuffle(shuffled)
        for node in shuffled:
            if node_types[node] not in type_to_offsets:
                continue
            walks.append(_walk(node, type_to_offsets[node_types[node]][0]))

    # ------------------------------------------------------------------ #
    # min_count filtering                                                  #
    # ------------------------------------------------------------------ #
    node_freq: Counter = Counter(n for walk in walks for n in walk)
    effective_min = max(min_count, 1)
    valid_nodes: set[str] = {n for n, c in node_freq.items() if c >= effective_min}
    if not valid_nodes:
        valid_nodes = set(nodes)

    # ------------------------------------------------------------------ #
    # Per-type negative-sampling tables  (distribution ∝ freq^0.75)      #
    # ------------------------------------------------------------------ #
    type_neg_pool: dict[str, np.ndarray] = {}
    type_neg_prob: dict[str, np.ndarray] = {}

    for t, t_nodes in type_to_nodes.items():
        valid_t = [n for n in t_nodes if n in valid_nodes]
        if not valid_t:
            continue
        freqs = np.array([node_freq.get(n, 1) for n in valid_t], dtype=float)
        probs = freqs ** 0.75
        probs /= probs.sum()
        type_neg_pool[t] = np.array([node_to_idx[n] for n in valid_t], dtype=np.intp)
        type_neg_prob[t] = probs

    # ------------------------------------------------------------------ #
    # Embedding matrices                                                   #
    # ------------------------------------------------------------------ #
    scale = 1.0 / embedding_dim
    W = rng.uniform(-scale, scale, (n_nodes, embedding_dim))
    W_ctx = np.zeros((n_nodes, embedding_dim), dtype=float)

    def _sigmoid(x: np.ndarray) -> np.ndarray:
        return 1.0 / (1.0 + np.exp(-np.clip(x, -10.0, 10.0)))

    # ------------------------------------------------------------------ #
    # Heterogeneous skip-gram — naive per-pair loop (no optimisations)    #
    #                                                                      #
    # B1: labels and all_ctx_idx allocated fresh on every pair.           #
    # B2: rng.choice called once per pair (not batched by type).          #
    # B3: np.concatenate, np.zeros, np.outer called inside the loop.      #
    # ------------------------------------------------------------------ #
    for epoch in range(epochs):
        epoch_rng = np.random.default_rng(42 + epoch)
        for walk in walks:
            for center_pos, center_node in enumerate(walk):
                if center_node not in valid_nodes:
                    continue
                center_idx = node_to_idx[center_node]
                ctx_start = max(0, center_pos - window_size)
                ctx_end = min(len(walk), center_pos + window_size + 1)

                for ctx_pos in range(ctx_start, ctx_end):
                    if ctx_pos == center_pos:
                        continue
                    ctx_node = walk[ctx_pos]
                    if ctx_node not in valid_nodes:
                        continue

                    ctx_idx = node_to_idx[ctx_node]
                    ctx_type = node_types[ctx_node]

                    pool = type_neg_pool.get(ctx_type)
                    probs = type_neg_prob.get(ctx_type)
                    if pool is None or len(pool) < 2:
                        continue

                    # B2: one rng.choice per pair
                    neg_local = epoch_rng.choice(
                        len(pool), size=negative_samples, p=probs, replace=True
                    )
                    neg_indices = pool[neg_local]

                    # B3: per-pair allocations
                    all_ctx_idx = np.concatenate([[ctx_idx], neg_indices])
                    labels = np.zeros(1 + negative_samples)
                    labels[0] = 1.0

                    center_vec = W[center_idx]
                    ctx_vecs = W_ctx[all_ctx_idx]
                    sig = _sigmoid(ctx_vecs @ center_vec)

                    err = sig - labels
                    grad_center = err @ ctx_vecs

                    # B3: np.outer per pair
                    np.add.at(W_ctx, all_ctx_idx, -learning_rate * np.outer(err, center_vec))
                    W[center_idx] -= learning_rate * grad_center

    # ------------------------------------------------------------------ #
    # L2-normalise input embeddings                                       #
    # ------------------------------------------------------------------ #
    norms = np.linalg.norm(W, axis=1, keepdims=True)
    nonzero = norms.squeeze() > 0
    W[nonzero] = W[nonzero] / norms[nonzero]

    # ------------------------------------------------------------------ #
    # Build output DataFrame                                               #
    # ------------------------------------------------------------------ #
    rows = []
    for node in nodes:
        idx = node_to_idx[node]
        row: dict = {"node_id": node, "node_type": node_types[node]}
        for dim_idx in range(embedding_dim):
            row[f"emb_{dim_idx + 1}"] = float(W[idx, dim_idx])
        rows.append(row)

    embeddings_df = pd.DataFrame(rows)

    if output_path is None:
        output_path = (
            Path(__file__).resolve().parent
            / "results"
            / "embedding"
            / "metapath2vec_pp_v0_embeddings.csv"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    embeddings_df.to_csv(output_path, index=False)
    return embeddings_df, walks


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
    include_path_centralities: bool = True,
) -> tuple[nx.Graph, dict[str, float | int], Path]:
    """Build a drug-drug network from node2vec embeddings.

    Workflow:
    - keep only drug nodes from the embeddings file;
    - connect two drugs when cosine(drug_i, drug_j) >= threshold.

    The input CSV is expected to contain at least:
    - ``node_id`` column
    - one or more embedding columns named ``emb_*``
    """

    if cosine_threshold < 0.0 or cosine_threshold > 1.0:
        raise ValueError("cosine_threshold must be in [0, 1].")
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
    node_id_series = embeddings_df["node_id"].astype(str)
    if "node_type" in embeddings_df.columns:
        node_type_series = embeddings_df["node_type"].astype(str).str.lower()
    else:
        node_type_series = pd.Series("unknown", index=embeddings_df.index)
        node_type_series[node_id_series.str.startswith("Drug_")] = "drug"
        node_type_series[node_id_series.str.startswith("Gene_")] = "gene"
    node_attrs = embeddings_df.drop(columns=emb_cols).set_index("node_id")

    emb_matrix = embeddings_df[emb_cols].to_numpy(dtype=float)
    if emb_matrix.size == 0:
        raise ValueError("Embedding matrix is empty.")

    # Robust L2 normalization ensures cosine == dot product.
    norms = np.linalg.norm(emb_matrix, axis=1, keepdims=True)
    nonzero = norms.squeeze() > 0
    emb_matrix[nonzero] = emb_matrix[nonzero] / norms[nonzero]

    drug_indices = np.flatnonzero(node_type_series.to_numpy() == "drug")
    if drug_indices.size == 0:
        raise ValueError("Embedding file does not contain drug nodes.")

    drug_node_ids = [node_ids[int(idx)] for idx in drug_indices]
    drug_matrix = emb_matrix[drug_indices]

    graph = nx.Graph()
    for drug_id in drug_node_ids:
        attrs = node_attrs.loc[drug_id].to_dict() if drug_id in node_attrs.index else {}
        cleaned_attrs = {
            str(key): (value.item() if isinstance(value, np.generic) else value)
            for key, value in attrs.items()
            if not pd.isna(value)
        }
        graph.add_node(drug_id, **cleaned_attrs)

    n_drugs = len(drug_node_ids)
    for start in range(0, n_drugs, block_size):
        end = min(start + block_size, n_drugs)
        similarities = drug_matrix[start:end] @ drug_matrix.T
        for local_idx, source_idx in enumerate(range(start, end)):
            row = similarities[local_idx]
            row[: source_idx + 1] = -np.inf
            target_indices = np.flatnonzero(row >= cosine_threshold)
            for target_idx in target_indices:
                graph.add_edge(
                    drug_node_ids[source_idx],
                    drug_node_ids[int(target_idx)],
                    weight=float(row[target_idx]),
                )

    original_node_count = n_drugs
    retained_node_count = sum(1 for n in graph.nodes() if graph.degree(n) > 0)
    nodes_removed = original_node_count - retained_node_count
    potential_edges = retained_node_count * (retained_node_count - 1) // 2
    edges_filtered = potential_edges - graph.number_of_edges()

    filtering_stats = {
        "similarity_threshold": cosine_threshold,
        "nodes_removed": nodes_removed,
        "edges_filtered": edges_filtered,
        "original_node_count": original_node_count,
        "retained_node_count": retained_node_count,
        "potential_edges": potential_edges,
    }
    filtering_output_path = (
        Path(__file__).resolve().parent / "results" / "similarity" / "filtering.json"
    )
    filtering_output_path.parent.mkdir(parents=True, exist_ok=True)
    with filtering_output_path.open("w", encoding="utf-8") as fh:
        json.dump(filtering_stats, fh, indent=2)

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
        include_path_centralities=include_path_centralities,
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
            "similarity_metric": "cosine_on_drug_embeddings",
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
