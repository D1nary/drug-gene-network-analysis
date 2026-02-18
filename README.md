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

## Typical workflow

1. Place/update dataset in `data/`.
2. Run `complex.py` with desired CLI options.
3. Inspect artifacts in `results/` (JSON/CSV metrics + PNG figures).
4. Use `notes/` for report integration and interpretation material.
