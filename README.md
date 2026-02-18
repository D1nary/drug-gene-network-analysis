# Complex Network Analysis on ChG-InterDecagon

## Project objective

This project analyzes the **topological and functional organization** of the ChG-InterDecagon Chemical-Gene interaction dataset.

The workflow starts from the bipartite **drug-gene network**, then projects it into a **drug-drug similarity network** (Jaccard similarity on target sets) to study network-level properties such as density, connectivity, and modularity.

After that, the project detects **drug communities** with the Louvain method, evaluates whether communities are internally homogeneous or heterogeneous, and studies how communities interact.

Inside communities, the analysis shifts to genes via **co-occurrence networks** to detect co-targeting/redundancy patterns. A **spectral analysis** (normalized Laplacian) is also performed to characterize structural cohesion and possible modular substructures.

Finally, for selected communities, the project builds a **Directed Acyclic Graph (DAG)** based on target-set inclusion to reveal profile hierarchies (maximal profiles, incremental variants, and specific profiles).

## Repository structure

- `complex.py`: Main CLI pipeline. Loads and preprocesses data, builds all networks (drug-target, similarity, community, co-occurrence), runs spectral analysis, generates DAGs, and orchestrates saving + plotting outputs.
- `network.py`: Network construction logic:
  - bipartite drug-target graph
  - drug similarity graph (Jaccard threshold)
  - gene co-occurrence graph
  - Louvain-based community graph
  - target-inclusion DAG
  - transitive reduction and normalized Laplacian spectra.
- `saves.py`: Persistence and metrics utilities. Computes/saves global, node-level, edge-level, filtering, community, co-occurrence, and DAG parameters into JSON/CSV artifacts.
- `visualizzation.py`: Visualization utilities for:
  - drug-gene spotlight subgraphs
  - random similarity snapshots
  - community network plots
  - per-community DAG figures.
- `data/`: Input dataset(s). Currently includes `ChG-InterDecagon_targets.csv.gz`.
- `results/`: Generated analysis outputs.
  - `results/drug_gene/`: Drug-gene spotlight image and related parameter tables.
  - `results/similarity/`: Similarity network image and parameter exports.
  - `results/community/`: Community visualizations and community-level metrics (`community_network_metrics/`).
  - `results/co_occurence/`: Co-occurrence parameter summaries.
  - `results/dag/communities/<id>/`: Per-community DAG inputs/outputs (community nodes, DAG metrics, DAG image).
- `notes/`: Project notes/report material (markdown, PDF, and reference images used for documentation).
- `vir_env/`: Local Python virtual environment (dependencies and interpreter packages).
- `__pycache__/`: Python bytecode cache files.

## Python files and their role

### `complex.py` (entry point)
- Provides CLI arguments to control pipeline behavior (selected networks, thresholds, target communities, etc.).
- Reads raw data, cleans it, and prints data diagnostics.
- Builds and visualizes a mid-degree drug-target spotlight.
- Builds the drug similarity network and persists metrics.
- Builds Louvain communities and saves:
  - community network metrics
  - selected community member metadata
  - density-filtered community summaries
  - normalized Laplacian spectra.
- Builds per-community gene co-occurrence networks and saves aggregated parameters.
- Builds and saves inclusion DAGs for selected communities.

### `network.py` (graph builders)
- `build_drug_target_network(...)`: Creates bipartite graph with drug and gene node partitions.
- `build_drug_similarity_network(...)`: Builds weighted drug-drug network using Jaccard similarity of target sets.
- `build_gene_cooccurrence_network(...)`: Connects genes co-targeted by the same drugs, with edge weights = shared drug count.
- `build_community_network(...)`: Runs Louvain and returns a community-level aggregated graph + membership map.
- `build_target_inclusion_dag(...)`: Builds DAG from target-set inclusion relations.
- `reduce_transitive_edges(...)`: Applies transitive reduction to keep only informative DAG edges.
- `compute_community_normalized_laplacian_spectra(...)`: Exports normalized Laplacian matrices and eigenvalues by community.

### `saves.py` (metrics + export layer)
- Computes global/topological summaries, per-node statistics, per-edge tables, and filtering metadata.
- Saves standard network artifacts (`global_parameters.json`, `node_parameters.csv`, `edge_parameters.csv`, `filtering.json`).
- Computes/saves co-occurrence summaries (global or per-community).
- Computes/saves community network metrics and Louvain summary.
- Exports dense-community member lists and community target-profile summaries.
- Saves per-community node metadata for downstream DAG construction.
- Computes/saves DAG global metrics and node-level DAG parameters.

### `visualizzation.py` (plotting layer)
- Produces a filtered drug-target spotlight view (high/mid/low-degree focus).
- Produces random drug-similarity snapshots (optionally colored by community membership).
- Visualizes the community graph with node size/color encoding structural properties.
- Visualizes per-community DAGs with hierarchy-aware layout.

## Typical workflow

1. Place/update dataset in `data/`.
2. Run `complex.py` with desired CLI options.
3. Inspect artifacts in `results/` (JSON/CSV metrics + PNG figures).
4. Use `notes/` for report integration and interpretation material.
