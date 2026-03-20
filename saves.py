#!/usr/bin/env python3
"""Utilities for persisting network parameter summaries."""

from __future__ import annotations

from pathlib import Path
import csv
import json
from collections import Counter
from statistics import mean, median, pstdev

import networkx as nx
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
RESULTS_DIR = PROJECT_ROOT / "results"
PARAMETERS_DIR = RESULTS_DIR
COMMUNITY_DIR = RESULTS_DIR / "community"
COMMUNITY_METRICS_DIR = COMMUNITY_DIR / "community_network_metrics"
CO_OCCURENCE_DIR = RESULTS_DIR / "co_occurence"
CO_OCCURENCE_PARAMETERS_DIR = CO_OCCURENCE_DIR / "parameters"


def _sanitize_label(label: str) -> str:
    # Keep output folder names filesystem-safe and deterministic.
    cleaned = label.strip().lower() or "graph"
    return "".join(
        ch if ch.isalnum() or ch in {"-", "_"} else "_"
        for ch in cleaned
    )


def compute_global_parameters(
    graph: nx.Graph,
    weight_attr: str = "weight",
) -> dict[str, float | int]:
    """Compute network-wide statistics."""

    # Core topology descriptors.
    node_count = graph.number_of_nodes()
    edge_count = graph.number_of_edges()
    density = nx.density(graph) if node_count > 1 else 0.0

    # Weight distribution summary (defaults to 1.0 for unweighted edges).
    weights = [
        float(data.get(weight_attr, 1.0))
        for _, _, data in graph.edges(data=True)
    ]
    if weights:
        mean_weight = float(np.mean(weights))
        median_weight = float(np.median(weights))
        min_weight = float(np.min(weights))
        max_weight = float(np.max(weights))
        std_weight = float(np.std(weights))
    else:
        mean_weight = median_weight = min_weight = max_weight = std_weight = 0.0

    components = list(nx.connected_components(graph))
    component_count = len(components)
    largest_component_size = max((len(c) for c in components), default=0)

    return {
        "node_count": node_count,
        "edge_count": edge_count,
        "density": density,
        "mean_weight": mean_weight,
        "median_weight": median_weight,
        "min_weight": min_weight,
        "max_weight": max_weight,
        "std_weight": std_weight,
        "component_count": component_count,
        "largest_component_size": largest_component_size,
    }


def save_similarity_global_parameters(
    graph: nx.Graph,
    output_path: Path | str | None = None,
    weight_attr: str = "weight",
    communities: list[set[object]] | None = None,
    community_seed: int = 42,
) -> tuple[dict[str, float | int], Path]:
    """Compute and save global parameters for a cosine similarity network."""

    components = list(nx.connected_components(graph))
    n_components = len(components)
    giant_component_nodes = max(components, key=len) if components else set()
    giant_component_size = len(giant_component_nodes)
    giant_component = (
        graph.subgraph(giant_component_nodes).copy()
        if giant_component_nodes
        else nx.Graph()
    )

    if giant_component_size > 1:
        diameter = int(nx.diameter(giant_component))
        avg_path_length = float(nx.average_shortest_path_length(giant_component))
    else:
        diameter = 0
        avg_path_length = 0.0

    if graph.number_of_nodes() > 0:
        communities = communities or compute_similarity_communities(
            graph,
            weight_attr=weight_attr,
            seed=community_seed,
        )
        modularity = float(
            nx.algorithms.community.modularity(
                graph,
                communities,
                weight=weight_attr,
            )
        )
        n_communities = len(communities)
    else:
        modularity = 0.0
        n_communities = 0

    global_params: dict[str, float | int] = {
        "n_nodes": int(graph.number_of_nodes()),
        "n_edges": int(graph.number_of_edges()),
        "density": float(nx.density(graph)) if graph.number_of_nodes() > 1 else 0.0,
        "n_components": int(n_components),
        "giant_component_size": int(giant_component_size),
        "avg_clustering": float(nx.average_clustering(graph, weight=weight_attr))
        if graph.number_of_nodes() > 0
        else 0.0,
        "diameter": int(diameter),
        "avg_path_length": float(avg_path_length),
        "modularity": float(modularity),
        "n_communities": int(n_communities),
    }

    if output_path is None:
        output_path = RESULTS_DIR / "similarity" / "global_parameters.json"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(global_params, fh, indent=2, ensure_ascii=False)

    return global_params, output_path


def compute_similarity_communities(
    graph: nx.Graph,
    weight_attr: str = "weight",
    seed: int = 42,
) -> list[set[object]]:
    """Return Louvain communities for the similarity graph."""

    if graph.number_of_nodes() == 0:
        return []
    if graph.number_of_edges() == 0:
        return [{node} for node in sorted(graph.nodes())]
    return nx.algorithms.community.louvain_communities(
        graph,
        weight=weight_attr,
        seed=seed,
    )


def compute_similarity_node_metrics(
    graph: nx.Graph,
    weight_attr: str = "weight",
    community_seed: int = 42,
    communities: list[set[object]] | None = None,
    include_path_centralities: bool = True,
) -> pd.DataFrame:
    """Compute per-node metrics required for similarity-network reporting."""

    columns = [
        "drug_id",
        "degree",
        "weighted_degree",
        "clustering",
        "community_id",
    ]
    if include_path_centralities:
        columns.insert(3, "betweenness")
        columns.insert(4, "closeness")
    nodes = list(graph.nodes())
    if not nodes:
        return pd.DataFrame(columns=columns)

    degree = dict(graph.degree())
    weighted_degree = dict(graph.degree(weight=weight_attr))
    clustering = nx.clustering(graph, weight=weight_attr)
    if include_path_centralities:
        # Keep shortest-path semantics unweighted; similarity weights are affinities.
        betweenness = nx.betweenness_centrality(graph, weight=None)
        closeness = nx.closeness_centrality(graph, distance=None)
    else:
        betweenness = {}
        closeness = {}

    if graph.number_of_nodes() == 0:
        selected_communities: list[set[object]] = []
    else:
        selected_communities = communities or compute_similarity_communities(
            graph,
            weight_attr=weight_attr,
            seed=community_seed,
        )
    community_membership: dict[object, int] = {}
    for idx, community_nodes in enumerate(selected_communities):
        for node in community_nodes:
            community_membership[node] = idx

    records = []
    for node in nodes:
        record = {
            "drug_id": str(node),
            "degree": int(degree.get(node, 0)),
            "weighted_degree": float(weighted_degree.get(node, 0.0)),
            "clustering": float(clustering.get(node, 0.0)),
            "community_id": int(community_membership.get(node, -1)),
        }
        if include_path_centralities:
            record["betweenness"] = float(betweenness.get(node, 0.0))
            record["closeness"] = float(closeness.get(node, 0.0))
        records.append(record)

    return pd.DataFrame(records, columns=columns)


def save_similarity_node_metrics(
    graph: nx.Graph,
    output_path: Path | str | None = None,
    weight_attr: str = "weight",
    community_seed: int = 42,
    communities: list[set[object]] | None = None,
    include_path_centralities: bool = True,
) -> tuple[pd.DataFrame, Path]:
    """Compute and save per-node similarity metrics to CSV."""

    node_df = compute_similarity_node_metrics(
        graph,
        weight_attr=weight_attr,
        community_seed=community_seed,
        communities=communities,
        include_path_centralities=include_path_centralities,
    )

    if output_path is None:
        output_path = RESULTS_DIR / "similarity" / "node_metrics.csv"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    node_df.to_csv(output_path, index=False)

    return node_df, output_path


def compute_similarity_edge_list(
    graph: nx.Graph,
    weight_attr: str = "weight",
) -> pd.DataFrame:
    """Return one row per similarity edge with derived distance."""

    columns = ["drug_1", "drug_2", "similarity", "distance"]
    if graph.number_of_edges() == 0:
        return pd.DataFrame(columns=columns)

    records = []
    for source, target, data in graph.edges(data=True):
        similarity = float(data.get(weight_attr, 1.0))
        records.append(
            {
                "drug_1": str(source),
                "drug_2": str(target),
                "similarity": similarity,
                "distance": float(1.0 - similarity),
            }
        )

    edge_df = pd.DataFrame(records, columns=columns)
    return edge_df.sort_values(by=["drug_1", "drug_2"]).reset_index(drop=True)


def save_similarity_edge_list(
    graph: nx.Graph,
    output_path: Path | str | None = None,
    weight_attr: str = "weight",
) -> tuple[pd.DataFrame, Path]:
    """Compute and save the similarity edge list to CSV."""

    edge_df = compute_similarity_edge_list(
        graph,
        weight_attr=weight_attr,
    )

    if output_path is None:
        output_path = RESULTS_DIR / "similarity" / "edge_list.csv"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    edge_df.to_csv(output_path, index=False)

    return edge_df, output_path


def compute_node_parameters(
    graph: nx.Graph,
    weight_attr: str = "weight",
) -> pd.DataFrame:
    """Return node-level metrics as a DataFrame."""

    # Return a stable schema even for empty graphs.
    nodes = list(graph.nodes())
    if not nodes:
        return pd.DataFrame(
            columns=[
                "node",
                "degree",
                "weighted_degree",
                "clustering_coefficient",
                "betweenness_centrality",
                "closeness_centrality",
            ]
        )

    degree = dict(graph.degree())
    weighted_degree = dict(graph.degree(weight=weight_attr))
    clustering = nx.clustering(graph, weight=weight_attr)
    # Use unweighted betweenness/closeness to avoid interpreting similarity
    # weights as distances (which would invert their semantics).
    betweenness = nx.betweenness_centrality(graph, weight=None)
    closeness = nx.closeness_centrality(graph, distance=None)

    records = []
    # Aggregate per-node centrality and local cohesion metrics.
    for node in nodes:
        records.append(
            {
                "node": node,
                "degree": degree.get(node, 0),
                "weighted_degree": float(weighted_degree.get(node, 0.0)),
                "clustering_coefficient": float(clustering.get(node, 0.0)),
                "betweenness_centrality": float(betweenness.get(node, 0.0)),
                "closeness_centrality": float(closeness.get(node, 0.0)),
            }
        )

    return pd.DataFrame(records)


def compute_edge_parameters(
    graph: nx.Graph,
    weight_attr: str = "weight",
) -> pd.DataFrame:
    """Return edge-level metrics (weight only, no target-overlap)."""

    if graph.number_of_edges() == 0:
        return pd.DataFrame(columns=["source", "target", weight_attr])

    records = []
    # Export one row per edge for downstream tabular analysis.
    for u, v, data in graph.edges(data=True):
        records.append(
            {
                "source": u,
                "target": v,
                weight_attr: float(data.get(weight_attr, 1.0)),
            }
        )

    return pd.DataFrame(records)


def compute_filtering_parameters(
    graph: nx.Graph,
    similarity_threshold: float | None = None,
    nodes_removed: int | None = None,
    edges_filtered: int | None = None,
) -> dict[str, float | int | None]:
    """Describe preprocessing/filtering applied while building the graph."""

    # Prefer explicit overrides; otherwise read the metadata stored on the graph.
    similarity_threshold = (
        similarity_threshold
        if similarity_threshold is not None
        else graph.graph.get("similarity_threshold")
    )
    nodes_removed = (
        nodes_removed
        if nodes_removed is not None
        else graph.graph.get("removed_drugs")
    )
    edges_filtered = (
        edges_filtered
        if edges_filtered is not None
        else graph.graph.get("filtered_edges")
    )

    original_nodes = graph.graph.get("original_drug_count")
    retained_nodes = graph.graph.get("retained_drug_count")
    potential_edges = graph.graph.get("potential_edges")

    return {
        "similarity_threshold": similarity_threshold,
        "nodes_removed": nodes_removed,
        "edges_filtered": edges_filtered,
        "original_node_count": original_nodes,
        "retained_node_count": retained_nodes,
        "potential_edges": potential_edges,
    }


def save_network_parameters(
    graph: nx.Graph,
    label: str,
    output_root: Path | str | None = None,
    weight_attr: str = "weight",
    filtering_details: dict | None = None,
) -> dict[str, Path]:
    """Persist parameter groups (global/node/edge/filtering) in dedicated files."""

    filtering_details = filtering_details or {}
    output_root = Path(output_root) if output_root else PARAMETERS_DIR
    target_dir = output_root / _sanitize_label(label)
    target_dir.mkdir(parents=True, exist_ok=True)

    # Save each parameter family in a dedicated artifact for modular consumption.
    global_params = compute_global_parameters(graph, weight_attr=weight_attr)
    global_path = target_dir / "global_parameters.json"
    with global_path.open("w", encoding="utf-8") as fh:
        json.dump(global_params, fh, indent=2, ensure_ascii=False)

    node_df = compute_node_parameters(graph, weight_attr=weight_attr)
    node_path = target_dir / "node_parameters.csv"
    node_df.to_csv(node_path, index=False)

    edge_df = compute_edge_parameters(graph, weight_attr=weight_attr)
    edge_path = target_dir / "edge_parameters.csv"
    edge_df.to_csv(edge_path, index=False)

    filtering_params = compute_filtering_parameters(
        graph,
        similarity_threshold=filtering_details.get("similarity_threshold"),
        nodes_removed=filtering_details.get("nodes_removed"),
        edges_filtered=filtering_details.get("edges_filtered"),
    )
    filtering_path = target_dir / "filtering.json"
    with filtering_path.open("w", encoding="utf-8") as fh:
        json.dump(filtering_params, fh, indent=2, ensure_ascii=False)

    return {
        "global": global_path,
        "node": node_path,
        "edge": edge_path,
        "filtering": filtering_path,
    }


def compute_cooccurrence_parameters(
    graph: nx.Graph,
    weight_attr: str = "weight",
    weight_ge_threshold: int = 3,
) -> dict[str, object]:
    """Compute co-occurrence network metrics (excluding assortativity)."""

    # Global structure and connectivity statistics.
    node_count = graph.number_of_nodes()
    edge_count = graph.number_of_edges()
    density = nx.density(graph) if node_count > 1 else 0.0
    component_count = (
        nx.number_connected_components(graph)
        if node_count > 0
        else 0
    )
    largest_component_size = (
        max((len(c) for c in nx.connected_components(graph)), default=0)
        if node_count > 0
        else 0
    )
    global_clustering = nx.transitivity(graph) if node_count > 2 else 0.0

    # Weight-level summaries capture sparsity and interaction strength.
    weights = [
        float(data.get(weight_attr, 1.0))
        for _, _, data in graph.edges(data=True)
    ]
    if weights:
        min_weight = float(np.min(weights))
        median_weight = float(np.median(weights))
        mean_weight = float(np.mean(weights))
        max_weight = float(np.max(weights))
    else:
        min_weight = median_weight = mean_weight = max_weight = 0.0

    weight_eq_1 = 0
    weight_eq_2 = 0
    weight_ge_k = 0
    # Bucket edge weights to quantify how many links are weak vs recurrent.
    for weight in weights:
        weight_int = int(round(weight))
        if weight_int == 1:
            weight_eq_1 += 1
        elif weight_int == 2:
            weight_eq_2 += 1
        if weight_int >= weight_ge_threshold:
            weight_ge_k += 1

    if edge_count > 0:
        pct_eq_1 = (weight_eq_1 / edge_count) * 100.0
        pct_eq_2 = (weight_eq_2 / edge_count) * 100.0
        pct_ge_k = (weight_ge_k / edge_count) * 100.0
    else:
        pct_eq_1 = pct_eq_2 = pct_ge_k = 0.0

    return {
        "n_nodes": node_count,
        "n_edges": edge_count,
        "density": density,
        "component_count": component_count,
        "giant_component_size": largest_component_size,
        "global_clustering_coefficient": global_clustering,
        "weight_distribution": {
            "min": min_weight,
            "median": median_weight,
            "mean": mean_weight,
            "max": max_weight,
        },
        "weight_sparsity": {
            "weight_eq_1_pct": pct_eq_1,
            "weight_eq_2_pct": pct_eq_2,
            "weight_ge_threshold": weight_ge_threshold,
            "weight_ge_threshold_pct": pct_ge_k,
        },
    }


def save_cooccurrence_parameters(
    graph: nx.Graph,
    output_dir: Path | str | None = None,
    weight_attr: str = "weight",
    weight_ge_threshold: int = 3,
) -> Path:
    """Save co-occurrence network parameters to a single file."""

    output_dir = Path(output_dir) if output_dir else CO_OCCURENCE_PARAMETERS_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    params = compute_cooccurrence_parameters(
        graph,
        weight_attr=weight_attr,
        weight_ge_threshold=weight_ge_threshold,
    )
    output_path = output_dir / "co_occurence_parameters.json"
    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(params, fh, indent=2, ensure_ascii=False)

    return output_path


def save_cooccurrence_parameters_by_community(
    community_parameters: dict[str, dict[str, object]],
    output_dir: Path | str | None = None,
) -> Path:
    """Save per-community co-occurrence parameters to a single JSON file."""

    output_dir = Path(output_dir) if output_dir else CO_OCCURENCE_PARAMETERS_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / "co_occurence_parameters.json"
    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(community_parameters, fh, indent=2, ensure_ascii=False)

    return output_path


def _community_weight_distribution(weights: list[float]) -> dict[str, float]:
    # Provide robust defaults so summaries stay valid when no inter-community edges exist.
    if not weights:
        return {
            "min_weight": 0.0,
            "max_weight": 0.0,
            "mean_weight": 0.0,
            "median_weight": 0.0,
            "std_weight": 0.0,
        }

    return {
        "min_weight": float(min(weights)),
        "max_weight": float(max(weights)),
        "mean_weight": float(mean(weights)),
        "median_weight": float(median(weights)),
        "std_weight": float(pstdev(weights)) if len(weights) > 1 else 0.0,
    }


def _community_global_parameters(community_graph: nx.Graph) -> dict[str, float | int]:
    # Summarize the community graph as a network-of-communities.
    community_count = community_graph.number_of_nodes()
    edge_count = community_graph.number_of_edges()
    density = nx.density(community_graph) if community_count > 1 else 0.0
    weights = [
        float(data.get("weight", 0.0))
        for _, _, data in community_graph.edges(data=True)
    ]
    weight_stats = _community_weight_distribution(weights)
    component_count = (
        nx.number_connected_components(community_graph)
        if community_count > 0
        else 0
    )

    return {
        "community_count": community_count,
        "inter_community_edge_count": edge_count,
        "density": density,
        "weight_distribution": weight_stats,
        "component_count": component_count,
    }


def _community_parameters_df(community_graph: nx.Graph) -> pd.DataFrame:
    # Build per-community rows with size, connectivity, and internal density.
    if community_graph.number_of_nodes() == 0:
        return pd.DataFrame(
            columns=[
                "community_id",
                "size",
                "degree",
                "weighted_degree",
                "clustering_coefficient",
                "density",
            ]
        )

    degree = dict(community_graph.degree())
    weighted_degree = dict(community_graph.degree(weight="weight"))
    rows = []
    for node, attrs in community_graph.nodes(data=True):
        size = int(attrs.get("size", 0))
        internal_edge_count = int(attrs.get("internal_edge_count", 0))
        internal_clustering = float(
            attrs.get("internal_clustering_coefficient", 0.0)
        )
        max_internal_edges = size * (size - 1) / 2 if size > 1 else 0
        density = (
            float(internal_edge_count) / max_internal_edges
            if max_internal_edges > 0
            else 0.0
        )
        rows.append(
            {
                "community_id": node,
                "size": size,
                "degree": int(degree.get(node, 0)),
                "weighted_degree": float(weighted_degree.get(node, 0.0)),
                "clustering_coefficient": internal_clustering,
                "density": density,
            }
        )

    return pd.DataFrame(rows)


def _community_edge_parameters_df(community_graph: nx.Graph) -> pd.DataFrame:
    # Transform inter-community links to tabular form and add a normalized distance proxy.
    if community_graph.number_of_edges() == 0:
        return pd.DataFrame(
            columns=["source", "target", "weight", "biological_distance"]
        )

    weights = [
        float(data.get("weight", 0.0))
        for _, _, data in community_graph.edges(data=True)
    ]
    max_weight = max(weights) if weights else 0.0
    records = []
    for u, v, data in community_graph.edges(data=True):
        weight = float(data.get("weight", 0.0))
        if max_weight > 0:
            biological_distance = max(0.0, 1.0 - (weight / max_weight))
        else:
            biological_distance = 0.0
        records.append(
            {
                "source": u,
                "target": v,
                "weight": weight,
                "biological_distance": biological_distance,
            }
        )
    return pd.DataFrame(records)


def _louvain_parameters(community_graph: nx.Graph) -> dict[str, float | int]:
    # Extract detection metadata plus descriptive statistics of community sizes.
    sizes = [
        int(attrs.get("size", len(attrs.get("members", []))))
        for _, attrs in community_graph.nodes(data=True)
    ]
    if sizes:
        min_size = min(sizes)
        max_size = max(sizes)
        mean_size = float(mean(sizes))
        median_size = float(median(sizes))
    else:
        min_size = max_size = 0
        mean_size = median_size = 0.0

    return {
        "method": community_graph.graph.get("method"),
        "resolution": community_graph.graph.get("resolution"),
        "modularity": community_graph.graph.get("modularity"),
        "min_community_size": min_size,
        "max_community_size": max_size,
        "mean_community_size": mean_size,
        "median_community_size": median_size,
    }


def save_community_data(
    community_graph: nx.Graph,
    output_dir: Path | str | None = None,
) -> dict[str, Path]:
    """Persist multi-file summary of the community network."""

    output_dir = Path(output_dir) if output_dir else COMMUNITY_METRICS_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save global stats, per-community table, edge table, and Louvain summary.
    global_params = _community_global_parameters(community_graph)
    global_path = output_dir / "community_global_parameters.json"
    with global_path.open("w", encoding="utf-8") as fh:
        json.dump(global_params, fh, indent=2, ensure_ascii=False)

    node_df = _community_parameters_df(community_graph)
    community_path = output_dir / "community_parameters.csv"
    node_df.to_csv(community_path, index=False)

    edge_df = _community_edge_parameters_df(community_graph)
    edge_path = output_dir / "community_edge_parameters.csv"
    edge_df.to_csv(edge_path, index=False)

    louvain_summary = _louvain_parameters(community_graph)
    louvain_path = output_dir / "louvain_parameters.json"
    with louvain_path.open("w", encoding="utf-8") as fh:
        json.dump(louvain_summary, fh, indent=2, ensure_ascii=False)

    print(f"Community network metrics written to {output_dir}")

    return {
        "global": global_path,
        "community": community_path,
        "edge": edge_path,
        "louvain": louvain_path,
    }


def save_community_normalized_laplacian_spectra(
    spectra: dict[str, dict[str, object]],
    output_dir: Path | str,
) -> Path:
    """Save normalized Laplacian matrices and eigen spectra per community."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "normalized_laplacian_spectra.json"

    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(spectra, fh, indent=2, ensure_ascii=False)

    return output_path


def _community_density(size: int, internal_edge_count: int) -> float:
    # Density of an undirected simple graph induced by community members.
    max_internal_edges = size * (size - 1) / 2 if size > 1 else 0
    if max_internal_edges == 0:
        return 0.0
    return float(internal_edge_count) / max_internal_edges


def _format_count(value: int) -> str:
    return f"{value:,}".replace(",", ".")


def save_community_members_by_density(
    community_graph: nx.Graph,
    density_threshold: float,
    min_size: int,
    output_dir: Path | str | None = None,
) -> Path:
    """Save community member lists for communities over density and size limits."""

    output_dir = Path(output_dir) if output_dir else COMMUNITY_METRICS_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    density_tag = f"{density_threshold:g}".replace(".", "p")
    output_path = output_dir / (
        f"community_members_density_gt_{density_tag}_size_ge_{min_size}.csv"
    )

    community_members: dict[str, list[str]] = {}
    # Keep only dense and sufficiently large communities.
    for node, attrs in community_graph.nodes(data=True):
        size = int(attrs.get("size", 0))
        internal_edge_count = int(attrs.get("internal_edge_count", 0))
        density = _community_density(size, internal_edge_count)
        if density < density_threshold or size < min_size:
            continue
        members = attrs.get("members", [])
        if not isinstance(members, (list, tuple, set)):
            members = list(members) if members else []
        community_members[node] = [str(member) for member in sorted(members)]

    with output_path.open("w", encoding="utf-8", newline="") as fh:
        # Write one column per community and one row per member rank.
        fieldnames = list(community_members.keys())
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        max_len = max((len(members) for members in community_members.values()), default=0)
        for idx in range(max_len):
            row = {
                community_id: members[idx] if idx < len(members) else ""
                for community_id, members in community_members.items()
            }
            writer.writerow(row)

    return output_path


def save_community_profile_summary(
    community_graph: nx.Graph,
    drug_targets: dict[str, list[int] | set[int] | tuple[int, ...]],
    density_threshold: float,
    min_size: int,
    output_dir: Path | str | None = None,
) -> Path:
    """Save per-community unique target profile counts with filters applied."""

    output_dir = Path(output_dir) if output_dir else COMMUNITY_METRICS_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    density_tag = f"{density_threshold:g}".replace(".", "p")
    output_path = output_dir / (
        f"community_profile_summary_density_ge_{density_tag}_size_ge_{min_size}.json"
    )

    records: list[dict[str, object]] = []
    communities = sorted(community_graph.nodes())
    # For each eligible community, count distinct and shared target signatures.
    for node in communities:
        attrs = community_graph.nodes[node]
        size = int(attrs.get("size", 0))
        internal_edge_count = int(attrs.get("internal_edge_count", 0))
        density = _community_density(size, internal_edge_count)
        if density < density_threshold or size < min_size:
            continue

        members = attrs.get("members", [])
        if not isinstance(members, (list, tuple, set)):
            members = list(members) if members else []

        profiles = []
        for member in members:
            targets = drug_targets.get(str(member))
            if targets is None:
                continue
            profile = tuple(sorted(targets))
            profiles.append(profile)

        profile_counts = Counter(profiles)
        unique_profiles = len(profile_counts)
        shared_profiles = sum(1 for count in profile_counts.values() if count > 1)
        if profile_counts:
            most_profile, most_count = profile_counts.most_common(1)[0]
            most_genes = len(most_profile)
        else:
            most_count = 0
            most_genes = 0

        records.append(
            {
                "community_id": str(node),
                "density": density,
                "size": size,
                "unique_profiles": unique_profiles,
                "shared_profiles": shared_profiles,
                "most_frequent_profile": {
                    "drug_count": most_count,
                    "gene_count": most_genes,
                },
            }
        )

    payload = {
        "filters": {
            "density_gte": density_threshold,
            "size_gte": min_size,
        },
        "communities": records,
    }

    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)

    return output_path


def save_community_node_info(
    community_id: int,
    members: list[str] | set[str] | tuple[str, ...],
    similarity_graph: nx.Graph,
    total_targets: int,
    output_dir: Path | str,
) -> Path:
    """Save per-drug metadata for a single community into a CSV file."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = []
    # Normalize member metadata from similarity graph attributes into tabular rows.
    for member in sorted(members):
        attrs = similarity_graph.nodes.get(member, {})
        drug_id = attrs.get("original_id")
        if drug_id is None and isinstance(member, str) and member.startswith("Drug_"):
            try:
                drug_id = int(member.replace("Drug_", ""))
            except ValueError:
                drug_id = None

        drug_name = attrs.get("drug_name") or attrs.get("name")
        targets = attrs.get("targets")
        if targets is None:
            target_list = None
            n_target = None
            target_coverage = None
        else:
            target_values = list(targets)
            try:
                target_values = sorted(target_values)
            except TypeError:
                pass
            target_list = json.dumps(target_values)
            n_target = len(targets)
            target_coverage = (
                float(n_target / total_targets) if total_targets > 0 else None
            )

        records.append(
            {
                "community_id": community_id,
                "drug_id": drug_id,
                "drug_name": drug_name,
                "n_target": n_target,
                "target_list": target_list,
                "target_coverage": target_coverage,
            }
        )

    df = pd.DataFrame(
        records,
        columns=[
            "community_id",
            "drug_id",
            "drug_name",
            "n_target",
            "target_list",
            "target_coverage",
        ],
    )
    output_path = output_dir / "community_nodes.csv"
    df.to_csv(output_path, index=False)
    return output_path


def compute_dag_global_parameters(graph: nx.DiGraph) -> dict[str, float | int]:
    """Compute global metrics for a directed acyclic graph."""

    n_nodes = graph.number_of_nodes()
    n_edges = graph.number_of_edges()
    density_dag = nx.density(graph) if n_nodes > 1 else 0.0

    in_degrees = dict(graph.in_degree())
    out_degrees = dict(graph.out_degree())
    n_sources = sum(1 for degree in in_degrees.values() if degree == 0)
    n_sinks = sum(1 for degree in out_degrees.values() if degree == 0)

    # Estimate hierarchy depth via longest directed path in the acyclic graph.
    if n_nodes == 0 or n_edges == 0:
        max_depth = 0
    else:
        try:
            # NetworkX renamed/relocated this in newer versions.
            if hasattr(nx.algorithms.dag, "dag_longest_path_length"):
                max_depth = nx.algorithms.dag.dag_longest_path_length(graph)
            elif hasattr(nx.algorithms.dag, "longest_path_length"):
                max_depth = nx.algorithms.dag.longest_path_length(graph)
            else:
                max_depth = nx.dag_longest_path_length(graph)
        except nx.NetworkXUnfeasible:
            max_depth = 0

    return {
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "density_dag": float(density_dag),
        "n_sources": n_sources,
        "n_sinks": n_sinks,
        "max_depth": int(max_depth),
    }


def compute_dag_node_parameters(graph: nx.DiGraph) -> pd.DataFrame:
    """Return DAG node metrics as a DataFrame."""

    if graph.number_of_nodes() == 0:
        return pd.DataFrame(
            columns=[
                "drug_id",
                "in_degree",
                "out_degree",
                "degree_ratio",
                "topological_level",
            ]
        )

    levels: dict[object, int] = {}
    try:
        # Assign level as the longest predecessor chain ending at each node.
        for node in nx.topological_sort(graph):
            preds = list(graph.predecessors(node))
            if not preds:
                levels[node] = 0
            else:
                levels[node] = 1 + max(levels[pred] for pred in preds)
    except nx.NetworkXUnfeasible:
        levels = {node: 0 for node in graph.nodes()}

    records = []
    for node in graph.nodes():
        in_degree = int(graph.in_degree(node))
        out_degree = int(graph.out_degree(node))
        degree_ratio = float(out_degree / (in_degree + 1))
        records.append(
            {
                "drug_id": node,
                "in_degree": in_degree,
                "out_degree": out_degree,
                "degree_ratio": degree_ratio,
                "topological_level": int(levels.get(node, 0)),
            }
        )

    return pd.DataFrame(records)


def save_dag_parameters(
    graph: nx.DiGraph,
    output_dir: Path | str,
) -> dict[str, Path]:
    """Persist DAG global and node parameters under the given directory."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Persist both graph-level and node-level DAG diagnostics.
    global_params = compute_dag_global_parameters(graph)
    global_path = output_dir / "dag_global_parameters.json"
    with global_path.open("w", encoding="utf-8") as fh:
        json.dump(global_params, fh, indent=2, ensure_ascii=False)

    node_df = compute_dag_node_parameters(graph)
    node_path = output_dir / "dag_node_parameters.csv"
    node_df.to_csv(node_path, index=False)

    return {
        "global": global_path,
        "node": node_path,
    }
