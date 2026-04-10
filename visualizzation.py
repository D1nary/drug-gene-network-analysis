#!/usr/bin/env python3
"""Visualization utilities for the Drug–Target bipartite network."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Optional
import math
import random

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
RESULTS_DIR = PROJECT_ROOT / "results"
DRUG_GENE_DIR = RESULTS_DIR / "drug_gene"
SIMILARITY_DIR = RESULTS_DIR / "similarity"
COMMUNITY_DIR = RESULTS_DIR / "community"
MAX_DRUGS_PER_VIEW = 50


def _sanitize_title(title: str) -> str:
    # Convert titles to filesystem-safe filenames.
    return "".join(
        c if c.isalnum() or c in {"-", "_"} else "_"
        for c in title.lower()
    )


def _prepare_output_path(title: str) -> Path:
    DRUG_GENE_DIR.mkdir(parents=True, exist_ok=True)
    return DRUG_GENE_DIR / f"{_sanitize_title(title)}.png"


def _prepare_similarity_output_path(title: str) -> Path:
    SIMILARITY_DIR.mkdir(parents=True, exist_ok=True)
    return SIMILARITY_DIR / f"{_sanitize_title(title)}.png"


def _prepare_community_output_path(title: str) -> Path:
    COMMUNITY_DIR.mkdir(parents=True, exist_ok=True)
    return COMMUNITY_DIR / f"{_sanitize_title(title)}.png"


def _select_drug_gene_subgraph(
    graph: nx.Graph,
    max_nodes: int,
    focus: str,
) -> tuple[nx.Graph, list, list, str]:
    # Split bipartite partitions from node attributes.
    drug_nodes = [
        n for n, data in graph.nodes(data=True)
        if data.get("bipartite") == "drug"
    ]
    gene_nodes = [
        n for n, data in graph.nodes(data=True)
        if data.get("bipartite") == "gene"
    ]

    if not drug_nodes or not gene_nodes:
        raise ValueError(
            "Graph does not appear to be bipartite with 'drug' and 'gene' partitions."
        )

    focus = focus.lower()
    if focus not in {"high", "mid", "low"}:
        raise ValueError("focus must be one of: 'high', 'mid', 'low'.")

    if focus == "high":
        # Prioritize highly connected drugs.
        candidate_drugs = [n for n in drug_nodes if graph.degree[n] > 10]
        sort_reverse = True
        degree_label = "High-degree"
        empty_msg = "No drug nodes with more than 10 connections were found in the graph."
    elif focus == "mid":
        # Keep medium-connectivity drugs for a balanced spotlight.
        candidate_drugs = [
            n for n in drug_nodes
            if 5 <= graph.degree[n] <= 15
        ]
        sort_reverse = True
        degree_label = "Mid-degree"
        empty_msg = "No drug nodes with between 5 and 15 connections were found in the graph."
    else:
        # Focus on sparse-profile drugs.
        candidate_drugs = [n for n in drug_nodes if graph.degree[n] <= 10]
        sort_reverse = False
        degree_label = "Low-degree"
        empty_msg = "No drug nodes with 10 or fewer connections were found in the graph."

    if not candidate_drugs:
        raise ValueError(empty_msg)

    candidate_drugs.sort(key=lambda n: graph.degree[n], reverse=sort_reverse)

    # Cap the number of drugs and fill the remaining budget with their genes.
    max_drug_limit = max(1, max_nodes - 1)
    max_drug_quota = min(len(candidate_drugs), MAX_DRUGS_PER_VIEW, max_drug_limit)
    selected_drugs = candidate_drugs[:max_drug_quota]

    max_gene_quota = max(1, max_nodes - len(selected_drugs))
    selected_genes: set = set()

    for drug in selected_drugs:
        neighbors = [
            neighbor for neighbor in graph.neighbors(drug)
            if graph.nodes[neighbor].get("bipartite") == "gene"
        ]
        # Prefer high-degree genes to retain informative interactions.
        neighbors.sort(key=lambda n: graph.degree[n], reverse=True)

        for neighbor in neighbors:
            if neighbor not in selected_genes:
                selected_genes.add(neighbor)

            if len(selected_genes) >= max_gene_quota:
                break

        if len(selected_genes) >= max_gene_quota:
            break

    nodes_to_keep = set(selected_drugs) | selected_genes
    subgraph = graph.subgraph(nodes_to_keep).copy()

    sub_drugs = [
        n for n, data in subgraph.nodes(data=True)
        if data.get("bipartite") == "drug"
    ]
    sub_genes = [
        n for n, data in subgraph.nodes(data=True)
        if data.get("bipartite") == "gene"
    ]

    return subgraph, sub_drugs, sub_genes, degree_label


def visualize_random_drug_target_subgraph(
    graph: nx.Graph,
    max_nodes: int = 200,
    title: Optional[str] = None,
    focus: str = "high",
    drop_isolates: bool = True,
    keep_largest_component: bool = True,
) -> nx.Graph:
    
    """Visualize a subgraph of the Drug–Target network filtered by drug degree.

    Depending on ``focus`` the function selects either the drug hubs
    (>10 targets) or the low-degree drugs (<=10 targets) and includes
    their gene neighbors, capped by ``max_nodes``. This keeps the
    figure readable while surfacing specific interaction patterns.

    Parameters
    ----------
    graph : nx.Graph
        The full bipartite Drug–Target network, as returned by
        `build_drug_target_network` in network.py.
    max_nodes : int, optional
        Maximum total number of nodes to include in the visualization
        (both drugs and genes). Default is 200.
    title : str, optional
        Custom title for the plot. If None, a default title is used.
    focus : {"high", "mid", "low"}
        Select highly connected drugs (>10 edges), mid-degree drugs (5–15 edges),
        or low-degree drugs (<=10 edges).

    Returns
    -------
    nx.Graph
        The induced subgraph that has been visualized.
    """

    subgraph, sub_drugs, sub_genes, degree_label = _select_drug_gene_subgraph(
        graph,
        max_nodes,
        focus,
    )

    # Optionally simplify the rendering by removing isolates/components.
    if drop_isolates:
        isolates = list(nx.isolates(subgraph))
        if isolates:
            subgraph.remove_nodes_from(isolates)

    if keep_largest_component and subgraph.number_of_nodes() > 0:
        components = list(nx.connected_components(subgraph))
        if len(components) > 1:
            largest = max(components, key=len)
            subgraph = subgraph.subgraph(largest).copy()

    sub_drugs = [
        n for n, data in subgraph.nodes(data=True)
        if data.get("bipartite") == "drug"
    ]
    sub_genes = [
        n for n, data in subgraph.nodes(data=True)
        if data.get("bipartite") == "gene"
    ]

    # Force-directed layout for a more organic positioning.
    node_count = max(1, subgraph.number_of_nodes())
    desired_spacing = 5.0 / math.sqrt(node_count)
    pos = nx.spring_layout(subgraph, seed=42, k=desired_spacing, iterations=300)

    # Redistribute gene nodes around their connected drug to avoid overlap
    gene_radius = 0.09
    placed_genes: set = set()

    for drug in sub_drugs:
        neighbors = [
            neighbor
            for neighbor in subgraph.neighbors(drug)
            if neighbor in sub_genes
        ]
        count = len(neighbors)
        if count == 0:
            continue

        for idx, gene in enumerate(neighbors):
            if gene in placed_genes:
                continue
            angle = 2 * math.pi * idx / count
            offset_x = gene_radius * math.cos(angle)
            offset_y = gene_radius * math.sin(angle)
            base_x, base_y = pos[drug]
            pos[gene] = (base_x + offset_x, base_y + offset_y)
            placed_genes.add(gene)

    # Create the figure
    plt.figure(figsize=(10, 8))

    # Draw drug and gene nodes with different colors and sizes
    nx.draw_networkx_nodes(
        subgraph,
        pos,
        nodelist=sub_drugs,
        node_color="tab:blue",
        node_size=80,
        label="Drugs",
        alpha=0.8,
    )
    nx.draw_networkx_nodes(
        subgraph,
        pos,
        nodelist=sub_genes,
        node_color="tab:orange",
        node_size=60,
        label="Genes",
        alpha=0.8,
    )

    # Draw edges
    nx.draw_networkx_edges(
        subgraph,
        pos,
        width=0.5,
        alpha=0.4,
    )

    # Nothing (or very few) labels to keep the plot readable
    # If you want labels, you can uncomment:
    # nx.draw_networkx_labels(subgraph, pos, font_size=6)

    if title is None:
        title = (
            f"{degree_label} Drug–Target subgraph "
            f"({subgraph.number_of_nodes()} nodes, "
            f"{subgraph.number_of_edges()} edges)"
        )

    plt.title(title)
    plt.legend()
    plt.axis("off")
    plt.tight_layout()

    output_path = _prepare_output_path(title)

    # Save to the standard drug-gene output directory.
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Drug–Target subgraph saved to {output_path}")

    return subgraph


def visualize_similarity_subgraph(
    graph: Optional[nx.Graph] = None,
    title: Optional[str] = None,
    filename: Optional[str] = None,
    seed: Optional[int] = None,
    min_degree: int = 10,
    max_nodes: Optional[int] = None,
    community_membership: Optional[dict[str, str]] = None,
    output_dir: Optional[Path] = None,
    max_legend_items: int = 25,
    legend_columns: int = 2,
    legend_fontsize: int = 9,
    edge_list_path: Optional[Path] = None,
    node_metrics_path: Optional[Path] = None,
) -> nx.Graph:
    """Plot a random similarity-network subgraph from persisted CSV outputs.

    Parameters
    ----------
    graph : nx.Graph, optional
        Kept for backward compatibility. Visualization data is loaded from CSV files.
    title : str, optional
        Title for the figure; autogenerated if omitted.
    seed : int, optional
        Optional seed for reproducible node sampling and layout.
    min_degree : int, optional
        Minimum node degree required for inclusion. Default 10.
    max_nodes : int, optional
        Maximum number of randomly sampled nodes to visualize after degree filtering.
        Default ``None`` (no node cap).
    community_membership : dict, optional
        Mapping from node name to community label. When provided, nodes are colored
        by community instead of degree. If omitted, communities are loaded from
        ``node_metrics.csv`` when available.
    output_dir : Path, optional
        Directory where the figure should be saved. Defaults to the similarity folder.
    max_legend_items : int, optional
        Maximum number of communities to show in the legend. Useful to keep the legend
        readable when many communities exist. Set to 0 to hide the legend entries.
    legend_columns : int, optional
        Number of columns to use when drawing the legend (default 2).
    legend_fontsize : int, optional
        Font size for the legend labels (default 9). Increase for readability when
        the legend is narrow or long.
    edge_list_path : Path, optional
        Path to ``edge_list.csv``. Defaults to ``results/similarity/edge_list.csv``.
    node_metrics_path : Path, optional
        Path to ``node_metrics.csv``. Defaults to ``results/similarity/node_metrics.csv``.

    Returns
    -------
    nx.Graph
        The similarity subgraph that was visualized.
    """

    del graph  # Data source is the saved similarity CSV outputs.

    edge_csv = Path(edge_list_path) if edge_list_path else SIMILARITY_DIR / "edge_list.csv"
    node_csv = (
        Path(node_metrics_path)
        if node_metrics_path
        else SIMILARITY_DIR / "node_metrics.csv"
    )
    if not edge_csv.exists():
        raise FileNotFoundError(f"Similarity edge list not found: {edge_csv}")
    if not node_csv.exists():
        raise FileNotFoundError(f"Similarity node metrics not found: {node_csv}")

    edge_df = pd.read_csv(edge_csv)
    node_df = pd.read_csv(node_csv)

    required_edge_columns = {"drug_1", "drug_2", "similarity"}
    missing_edge_columns = required_edge_columns - set(edge_df.columns)
    if missing_edge_columns:
        raise ValueError(
            "edge_list.csv is missing required columns: "
            f"{sorted(missing_edge_columns)}"
        )

    required_node_columns = {"drug_id"}
    missing_node_columns = required_node_columns - set(node_df.columns)
    if missing_node_columns:
        raise ValueError(
            "node_metrics.csv is missing required columns: "
            f"{sorted(missing_node_columns)}"
        )

    loaded_graph = nx.Graph()
    for row in node_df.itertuples(index=False):
        loaded_graph.add_node(str(row.drug_id))

    for row in edge_df.itertuples(index=False):
        loaded_graph.add_edge(
            str(row.drug_1),
            str(row.drug_2),
            weight=float(row.similarity),
        )

    if loaded_graph.number_of_nodes() == 0:
        raise ValueError("Similarity graph is empty; cannot visualize.")

    nodes = sorted(n for n in loaded_graph.nodes() if loaded_graph.degree[n] >= min_degree)
    if not nodes:
        raise ValueError(
            f"Similarity graph has no nodes with degree >= {min_degree}; cannot visualize."
        )

    selected_nodes = nodes
    if max_nodes is not None and max_nodes <= 0:
        raise ValueError("max_nodes must be > 0 when provided.")

    if max_nodes is not None and len(nodes) > max_nodes:
        rng = random.Random(seed)
        selected_nodes = sorted(rng.sample(nodes, max_nodes))

    subgraph = loaded_graph.subgraph(selected_nodes).copy()

    node_count = subgraph.number_of_nodes()
    edge_count = subgraph.number_of_edges()

    print(f"Similarity nodes visualized: {node_count}")
    print(f"Similarity edges visualized: {edge_count}")

    layout_seed = seed if seed is not None else 42
    spacing = 3.0 / math.sqrt(max(1, node_count))
    pos = nx.spring_layout(subgraph, seed=layout_seed, k=spacing, iterations=200)

    # Scale node size by degree to emphasize local hubs.
    degrees = dict(subgraph.degree())
    max_degree = max(degrees.values()) or 1
    node_sizes = [
        40 + 80 * (degrees[node] / max_degree)
        for node in subgraph.nodes()
    ]

    edge_weights = [
        subgraph[u][v].get("weight", 1.0) for u, v in subgraph.edges()
    ]
    # Convert raw similarity weights to visible edge widths.
    normalized_widths = [
        0.5 + 2.0 * (weight - min(edge_weights, default=0))
        for weight in edge_weights
    ] if edge_weights else []

    fig, ax = plt.subplots(figsize=(12, 10))
    nx.draw_networkx_edges(
        subgraph,
        pos,
        alpha=0.6,
        width=normalized_widths if normalized_widths else 0.5,
        edge_color="gray",
        ax=ax,
    )
    legend_handles = None
    if community_membership is None and "community_id" in node_df.columns:
        community_membership = {
            str(row.drug_id): str(row.community_id)
            for row in node_df.itertuples(index=False)
        }

    if community_membership:
        # Color by community assignment and optionally build a compact legend.
        community_counts = Counter(
            community_membership.get(node, "Unassigned") for node in subgraph.nodes()
        )
        communities = [
            community for community, _ in sorted(
                community_counts.items(),
                key=lambda item: item[1],
                reverse=True,
            )
        ]
        cmap = plt.get_cmap("tab20")
        color_lookup = {
            community: cmap(idx % cmap.N) for idx, community in enumerate(communities)
        }
        node_colors = [
            color_lookup.get(community_membership.get(node, "Unassigned"), "tab:gray")
            for node in subgraph.nodes()
        ]
        node_collection = nx.draw_networkx_nodes(
            subgraph,
            pos,
            node_color=node_colors,
            node_size=node_sizes,
            alpha=0.9,
            edgecolors="black",
            linewidths=0.3,
            ax=ax,
        )
        if max_legend_items > 0:
            displayed = communities[:max_legend_items]
            legend_handles = [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=color_lookup[community],
                    markersize=7,
                    label=f"{community} (n={community_counts[community]})",
                )
                for community in displayed
            ]
            hidden = len(communities) - len(displayed)
            if hidden > 0:
                legend_handles.append(
                    Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        markerfacecolor="lightgray",
                        markersize=6,
                        label=f"{hidden} more communities",
                    )
                )
    else:
        # Fallback coloring: encode node degree with a continuous colormap.
        node_collection = nx.draw_networkx_nodes(
            subgraph,
            pos,
            node_color=list(degrees.values()),
            cmap="viridis",
            node_size=node_sizes,
            alpha=0.85,
            ax=ax,
        )
        sm = plt.cm.ScalarMappable(
            cmap="viridis",
            norm=plt.Normalize(vmin=0, vmax=max_degree),
        )
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, shrink=0.6)
        cbar.set_label("Node degree", rotation=270, labelpad=12)

    if title is None:
        title = (
            f"Drug Similarity Random Subgraph "
            f"({node_count} nodes, {edge_count} edges)"
        )

    ax.set_title(title)
    if legend_handles:
        ax.legend(
            handles=legend_handles,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            borderaxespad=0,
            fontsize=legend_fontsize,
            markerscale=1.0,
            handletextpad=0.4,
            columnspacing=0.8,
            ncol=max(1, legend_columns),
        )
    ax.axis("off")
    fig.tight_layout()

    target_dir = output_dir if output_dir is not None else SIMILARITY_DIR
    target_dir.mkdir(parents=True, exist_ok=True)
    file_stem = filename if filename is not None else _sanitize_title(title)
    output_path = target_dir / f"{file_stem}.png"
    # Persist figure to the selected directory.
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Similarity subgraph saved to {output_path}")

    return subgraph


def visualize_community_network(
    community_graph: nx.Graph,
    title: Optional[str] = None,
    seed: Optional[int] = None,
) -> nx.Graph:
    """Plot the community network produced by ``build_community_network``.

    Parameters
    ----------
    community_graph : nx.Graph
        Graph where each node represents a community and edges summarize
        cross-community connections.
    title : str, optional
        Title for the figure; autogenerated if omitted.
    seed : int, optional
        Seed forwarded to the ForceAtlas-style layout to keep runs reproducible.

    Returns
    -------
    nx.Graph
        The community graph, returned for convenience.
    """

    if community_graph.number_of_nodes() == 0:
        raise ValueError("Community graph is empty; cannot visualize.")

    # Layout the network of communities with edge weights influencing attraction.
    layout_seed = seed if seed is not None else 24
    node_count = community_graph.number_of_nodes()
    spacing = 3.0 / math.sqrt(max(1, node_count))
    pos = nx.spring_layout(
        community_graph,
        seed=layout_seed,
        k=spacing,
        weight="weight",
        iterations=200,
    )

    community_sizes = {
        node: data.get("size", 1) for node, data in community_graph.nodes(data=True)
    }
    # Scale node areas by community size.
    min_size = min(community_sizes.values()) if community_sizes else 1
    max_size = max(community_sizes.values()) if community_sizes else 1
    size_range = max(1, max_size - min_size)
    node_sizes = [
        200 + 300 * ((community_sizes[node] - min_size) / size_range)
        for node in community_graph.nodes()
    ]

    internal_weights = {
        node: community_graph.nodes[node].get("internal_weight", 0.0)
        for node in community_graph.nodes()
    }
    # Color communities by internal weight to highlight dense internal structure.
    color_values = list(internal_weights.values())
    color_min = min(color_values, default=0.0)
    color_max = max(color_values, default=1.0)
    if color_max == color_min:
        color_max = color_min + 1.0  # avoid zero range when building color map
    cmap = plt.cm.plasma

    edge_weights = [
        data.get("weight", 1.0) for _, _, data in community_graph.edges(data=True)
    ]
    # Normalize inter-community weights for consistent line thickness.
    if edge_weights:
        edge_min = min(edge_weights)
        edge_max = max(edge_weights)
        edge_range = max(1e-6, edge_max - edge_min)
        edge_widths = [
            0.3 + 2.7 * ((weight - edge_min) / edge_range)
            for weight in edge_weights
        ]
    else:
        edge_widths = []

    fig, ax = plt.subplots(figsize=(10, 8))
    nx.draw_networkx_edges(
        community_graph,
        pos,
        width=edge_widths if edge_widths else 0.5,
        alpha=0.3,
        edge_color="gray",
        ax=ax,
    )
    node_collection = nx.draw_networkx_nodes(
        community_graph,
        pos,
        node_color=color_values if color_values else "tab:blue",
        cmap=cmap if color_values else None,
        node_size=node_sizes,
        edgecolors="black",
        linewidths=0.5,
        alpha=0.9,
        ax=ax,
        vmin=color_min if color_values else None,
        vmax=color_max if color_values else None,
    )

    labels = {
        node: f"{node}\n(|C|={community_sizes.get(node, 0)})"
        for node in community_graph.nodes()
    }
    # Annotate each node with its label and cardinality.
    nx.draw_networkx_labels(
        community_graph,
        pos,
        labels=labels,
        font_size=7,
        font_weight="bold",
        ax=ax,
    )

    if color_values:
        norm = plt.Normalize(vmin=color_min, vmax=color_max)
        cbar = fig.colorbar(
            plt.cm.ScalarMappable(norm=norm, cmap=cmap),
            ax=ax,
            shrink=0.7,
        )
        cbar.set_label("Internal weight", rotation=270, labelpad=12)

    modularity = community_graph.graph.get("modularity")
    if title is None:
        title = f"Community Network ({node_count} communities)"
        if modularity is not None:
            title += f" – Modularity {modularity:.3f}"

    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()

    output_path = _prepare_community_output_path(title)
    # Save under the community visualization folder.
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Community network saved to {output_path}")

    return community_graph


def visualize_community_dag(
    dag: nx.DiGraph,
    community_id: int,
    output_dir: Path,
    title: Optional[str] = None,
    seed: Optional[int] = None,
) -> Path:
    """Render and save a DAG plot for a single community."""

    if dag.number_of_nodes() == 0:
        raise ValueError("DAG is empty; cannot visualize.")

    levels: dict[object, int] = {}
    is_dag = nx.is_directed_acyclic_graph(dag)
    if is_dag:
        # Place nodes by topological level to emphasize hierarchy depth.
        for node in nx.topological_sort(dag):
            preds = list(dag.predecessors(node))
            levels[node] = 0 if not preds else 1 + max(levels[p] for p in preds)
        groups: dict[int, list[object]] = {}
        for node, level in levels.items():
            groups.setdefault(level, []).append(node)

        pos: dict[object, tuple[float, float]] = {}
        for level, nodes in sorted(groups.items()):
            # Spread nodes horizontally inside each level band.
            nodes_sorted = sorted(nodes, key=lambda value: str(value))
            width = max(1, len(nodes_sorted) - 1)
            for idx, node in enumerate(nodes_sorted):
                x = (idx - width / 2.0) * 1.6
                y = -float(level) * 1.8
                pos[node] = (x, y)
    else:
        # Fallback layout if the graph is not acyclic.
        layout_seed = seed if seed is not None else 42
        pos = nx.spring_layout(dag, seed=layout_seed, k=1.2, iterations=200)
        levels = {node: 0 for node in dag.nodes()}

    out_degree = dict(dag.out_degree())
    max_out = max(out_degree.values(), default=1) or 1
    # Emphasize parent-like nodes with larger out-degree.
    node_sizes = [
        300 + 700 * (out_degree.get(node, 0) / max_out)
        for node in dag.nodes()
    ]
    level_values = [levels.get(node, 0) for node in dag.nodes()]
    edge_widths = [
        0.8 + 1.2 * float(data.get("set_difference", 1))
        for _, _, data in dag.edges(data=True)
    ]
    # Width can encode optional set-difference metadata when available.

    fig, ax = plt.subplots(figsize=(12, 8))
    nx.draw_networkx_edges(
        dag,
        pos,
        ax=ax,
        width=edge_widths if edge_widths else 0.8,
        alpha=0.5,
        arrows=True,
        arrowsize=14,
        arrowstyle="-|>",
        edge_color="gray",
        connectionstyle="arc3,rad=0.04",
    )
    nx.draw_networkx_nodes(
        dag,
        pos,
        ax=ax,
        node_size=node_sizes,
        node_color=level_values,
        cmap="viridis",
        alpha=0.95,
        edgecolors="black",
        linewidths=0.6,
    )
    nx.draw_networkx_labels(
        dag,
        pos,
        ax=ax,
        labels={node: str(node) for node in dag.nodes()},
        font_size=7,
    )

    if title is None:
        title = (
            f"Community {community_id} DAG "
            f"({dag.number_of_nodes()} nodes, {dag.number_of_edges()} edges)"
        )
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"community_{community_id}_dag.png"
    # Save the DAG image in the requested community folder.
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Community {community_id} DAG graph saved to {output_path}")
    return output_path
