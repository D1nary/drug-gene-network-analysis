#!/usr/bin/env python3
"""Utilities to build a Drug–Target bipartite network from 
the ChG-InterDecagon targets dataset."""

from __future__ import annotations

from itertools import combinations
import pandas as pd
import networkx as nx


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
        G.add_node(
            f"Drug_{drug_id}",
            bipartite="drug",
            original_id=int(drug_id),
            targets=sorted(targets),
        )

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


def build_gene_cooccurrence_network(
    df: pd.DataFrame,
    min_drugs_per_gene: int = 3,
    max_drugs_per_gene_percentile: float = 95.0,
    weight_threshold: int = 1,
) -> nx.Graph:
    """Create a weighted gene co-occurrence network via shared drugs."""

    required_cols = {"Drug", "Gene"}
    if not required_cols.issubset(df.columns):
        raise ValueError(
            f"DataFrame must contain columns {required_cols}, "
            f"but has {set(df.columns)}."
        )

    tmp = df[["Drug", "Gene"]].copy()
    tmp["Drug"] = tmp["Drug"].astype(int)
    tmp["Gene"] = tmp["Gene"].astype(int)
    tmp = tmp.drop_duplicates(subset=["Drug", "Gene"])

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

    shared_drug_counts: dict[tuple[int, int], int] = {}
    for _, group in filtered_tmp.groupby("Drug"):
        genes = sorted(set(group["Gene"].astype(int).tolist()))
        if len(genes) < 2:
            continue
        for gene_a, gene_b in combinations(genes, 2):
            key = (gene_a, gene_b)
            shared_drug_counts[key] = shared_drug_counts.get(key, 0) + 1

    filtered_edges = 0
    for (gene_a, gene_b), weight in shared_drug_counts.items():
        if weight_threshold is not None and weight < weight_threshold:
            filtered_edges += 1
            continue
        cooccurrence_graph.add_edge(
            f"Gene_{gene_a}",
            f"Gene_{gene_b}",
            weight=float(weight),
        )

    if weight_threshold is not None and weight_threshold > 1:
        isolated = [node for node, degree in cooccurrence_graph.degree() if degree == 0]
        if isolated:
            cooccurrence_graph.remove_nodes_from(isolated)

    cooccurrence_graph.graph.update(
        {
            "metric": "shared_drugs_count",
            "min_drugs_per_gene": min_drugs_per_gene,
            "max_drugs_per_gene_percentile": max_drugs_per_gene_percentile,
            "max_drugs_per_gene_value": percentile_value,
            "cooccurrence_weight_threshold": weight_threshold,
            "original_gene_count": original_gene_count,
            "retained_gene_count": retained_gene_count,
            "removed_genes": removed_genes,
            "filtered_edges": filtered_edges,
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
        for node in community:
            membership[node] = label

    inter_weights: dict[tuple[str, str], float] = {}
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
#prova
