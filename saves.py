#!/usr/bin/env python3
"""Utilities for persisting network parameter summaries."""

from __future__ import annotations

from pathlib import Path
import json
from statistics import mean, median, pstdev

import networkx as nx
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
RESULTS_DIR = PROJECT_ROOT / "results"
PARAMETERS_DIR = RESULTS_DIR / "network_parameters"
COMMUNITY_DIR = RESULTS_DIR / "community"
COMMUNITY_METRICS_DIR = COMMUNITY_DIR / "community_network_metrics"


def _sanitize_label(label: str) -> str:
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

    node_count = graph.number_of_nodes()
    edge_count = graph.number_of_edges()
    density = nx.density(graph) if node_count > 1 else 0.0

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


def compute_node_parameters(
    graph: nx.Graph,
    weight_attr: str = "weight",
) -> pd.DataFrame:
    """Return node-level metrics as a DataFrame."""

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


def _community_weight_distribution(weights: list[float]) -> dict[str, float]:
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
    if community_graph.number_of_nodes() == 0:
        return pd.DataFrame(
            columns=[
                "community_id",
                "size",
                "degree",
                "weighted_degree",
                "clustering_coefficient",
            ]
        )

    degree = dict(community_graph.degree())
    weighted_degree = dict(community_graph.degree(weight="weight"))
    clustering = nx.clustering(community_graph, weight="weight")

    rows = []
    for node, attrs in community_graph.nodes(data=True):
        rows.append(
            {
                "community_id": node,
                "size": int(attrs.get("size", 0)),
                "degree": int(degree.get(node, 0)),
                "weighted_degree": float(weighted_degree.get(node, 0.0)),
                "clustering_coefficient": float(clustering.get(node, 0.0)),
            }
        )

    return pd.DataFrame(rows)


def _community_edge_parameters_df(community_graph: nx.Graph) -> pd.DataFrame:
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
