# Complex Network Analysis on ChG-InterDecagon

## Project objective

This project analyzes the **topological and functional organization** of the ChG-InterDecagon Chemical-Gene interaction dataset.

The workflow starts from the bipartite **drug-gene network**, then projects it into a **drug-drug similarity network** (Jaccard similarity on target sets, or cosine similarity on learned embeddings) to study network-level properties such as density, connectivity, and modularity.

After that, the project detects **drug communities** with the Louvain method, evaluates whether communities are internally homogeneous or heterogeneous, and studies how communities interact.

Inside communities, the analysis shifts to genes via **co-occurrence networks** to detect co-targeting/redundancy patterns. A **spectral analysis** (normalized Laplacian) is also performed to characterize structural cohesion and possible modular substructures.

Finally, for selected communities, the project builds a **Directed Acyclic Graph (DAG)** based on target-set inclusion to reveal profile hierarchies (maximal profiles, incremental variants, and specific profiles).

A separate **CBP analysis** (Centrality-Based Profiling) computes betweenness centrality, closeness centrality, PageRank, and weighted degree on a user-supplied edge list.

---

## Repository structure

### Source files

- **`complex.py`** — Main CLI pipeline. Loads and preprocesses the dataset, builds all networks (drug-target, similarity, community, co-occurrence), runs spectral analysis, generates DAGs, and orchestrates saving and plotting outputs.
- **`network.py`** — Network construction and embedding logic:
  - Bipartite drug-target graph
  - All embedding algorithms (Node2Vec, Metapath2Vec variants — see below)
  - Drug-drug similarity graph (Jaccard threshold or cosine similarity on embeddings)
  - Gene co-occurrence graph
  - Louvain-based community graph
  - Target-inclusion DAG, transitive reduction, and normalized Laplacian spectra
- **`saves.py`** — Persistence and metrics utilities. Computes and saves global, node-level, edge-level, filtering, community, co-occurrence, and DAG parameters into JSON/CSV artifacts.
- **`visualizzation.py`** — Visualization utilities for:
  - Drug-gene spotlight subgraphs
  - Random similarity snapshots colored by community
  - Community network plots
  - Per-community DAG figures
- **`measurement.py`** — Standalone measurement module with its own CLI. Computes:
  - Per-community metrics (density, size, clustering coefficient, weighted degree) from Louvain assignments
  - CBP metrics (betweenness centrality, closeness centrality, PageRank, degree) from a generic edge list

### Data

- **`data/`** — Input dataset. Contains `ChG-InterDecagon_targets.csv.gz`.
- **`data_for_cbp/`** — Input edge list for CBP measurements (`edge_list.csv`). Must contain at least `drug_1` and `drug_2` columns; an optional `similarity` column is used as edge weight.

### Parameters

- **`parameters/node2vec/hyperparameters.json`** — Node2Vec hyperparameters.
- **`parameters/metapath2vec/hyperparameters.json`** — Metapath2Vec hyperparameters (shared by all metapath2vec variants).

### Results

- **`results/drug_gene/`** — Bipartite graph spotlight image (`bipartite_graph.png`) and graph statistics (`graph_stats.json`).
- **`results/embedding/`** — Embedding CSV files produced by the selected algorithm:
  - `node_embeddings.csv` (Node2Vec)
  - `metapath2vec_embeddings.csv` (Metapath2Vec pq-biased)
  - `metapath2vec_standard_embeddings.csv` (Metapath2Vec standard)
  - `metapath2vec_pp_embeddings.csv` (Metapath2Vec++)
  - `metapath2vec_pp_v0_embeddings.csv` (Metapath2Vec++ unoptimized reference)
- **`results/similarity/`** — Similarity network outputs:
  - `global_parameters.json` — Network-wide metrics (nodes, edges, density, modularity, diameter, etc.)
  - `node_metrics.csv` — Per-node metrics (degree, weighted degree, clustering, betweenness, closeness, community ID)
  - `edge_list.csv` — All similarity edges with weight and derived distance
  - `filtering.json` — Statistics on nodes/edges removed during threshold filtering
  - `community_measurement.csv` — Per-community aggregated metrics (density, size, clustering coefficient, weighted degree)
  - `similarity_network.png` — Visualization of the similarity network
- **`results/community/`** — Community visualizations and community-level metrics (`community_network_metrics/`).
- **`results/cbp_measurement/`** — CBP analysis outputs (one CSV per metric):
  - `betweenness_centrality.csv`
  - `closeness_centrality.csv`
  - `pagerank.csv`
  - `degree.csv`
- **`results/co_occurence/`** — Co-occurrence parameter summaries.
- **`results/dag/communities/<id>/`** — Per-community DAG inputs/outputs (community nodes, DAG metrics, DAG image).

---

## Embedding algorithms

Four embedding algorithms are available and selected via `--embedding-algorithm`:

| Algorithm | Flag | Description |
|---|---|---|
| Node2Vec | `node2vec` | Biased random walks (p/q) + Skip-gram via gensim |
| Metapath2Vec (pq) | `metapath2vec` | Metapath-guided walks with p/q bias + gensim |
| Metapath2Vec standard | `metapath2vec-standard` | Uniform metapath-guided walks + gensim |
| Metapath2Vec++ | `metapath2vec-pp` | Uniform metapath walks + heterogeneous negative sampling (custom NumPy SGD) |
| Metapath2Vec++ v0 | `metapath2vec-pp-v0` | Unoptimized reference implementation of Metapath2Vec++ (correctness baseline) |

All metapath2vec variants use the metapath defined in `parameters/metapath2vec/hyperparameters.json` (default: `["drug", "gene", "drug"]`). The `p` and `q` fields are only used by the `metapath2vec` (pq-biased) variant.

Metapath2Vec++ (`metapath2vec-pp`) differs from the standard variant in the training objective: negative samples are drawn exclusively from nodes of the same type as the context node (heterogeneous negative sampling), rather than from the full node vocabulary.

---

## Typical workflow

```bash
# 1. Run embedding generation
python complex.py --run-embedding --embedding-algorithm metapath2vec-pp-v0

# 2. Build the similarity network from embeddings
python complex.py --run-similarity --cosine-threshold 0.60

# 3. Visualize the similarity network
python complex.py --run-similarity-visualization

# 4. Compute per-community measurement metrics
python complex.py --run-measurement

# 5. Run CBP analysis on the CBP edge list
python complex.py --cbp

# Or run everything in a single call
python complex.py --run-embedding --run-similarity --run-similarity-visualization \
                  --run-measurement --cbp
```

`measurement.py` also exposes a standalone CLI:

```bash
# Community measurements
python measurement.py --node-metrics results/similarity/node_metrics.csv \
                      --edge-list results/similarity/edge_list.csv \
                      --output results/similarity/community_measurement.csv
```

---

## CLI reference (`complex.py`)

| Flag | Default | Description |
|---|---|---|
| `--embedding-algorithm` | `metapath2vec-pp-v0` | Embedding algorithm to use (`node2vec`, `metapath2vec`, `metapath2vec-standard`, `metapath2vec-pp`, `metapath2vec-pp-v0`) |
| `--data-path` | `data/ChG-InterDecagon_targets.csv.gz` | Path to the input dataset |
| `--cosine-threshold` | `0.60` | Cosine similarity threshold for the drug-drug network |
| `--run-embedding` / `--no-run-embedding` | disabled | Run embedding generation |
| `--run-graph-stats` / `--no-run-graph-stats` | disabled | Compute and save bipartite graph statistics |
| `--run-similarity` / `--no-run-similarity` | disabled | Build the similarity network from embeddings |
| `--run-similarity-visualization` / `--no-run-similarity-visualization` | disabled | Render the similarity network figure |
| `--similarity-path-centralities` / `--no-similarity-path-centralities` | disabled | Include betweenness/closeness centrality in node metrics |
| `--run-measurement` / `--no-run-measurement` | disabled | Compute per-community measurement metrics |
| `--cbp` / `--no-cbp` | enabled | Run CBP measurements from `data_for_cbp/edge_list.csv` |
| `--cbp-giant-component` / `--no-cbp-giant-component` | disabled | Restrict CBP metrics to the giant connected component |
| `--networks` | `similarity community cooccurence` | Which network types to build and save |
| `--community-density-threshold` | `0.99` | Density threshold for saving community member lists |
| `--community-min-size` | `10` | Minimum community size for member list export |
| `--community-ids` | `[21]` | Community IDs for DAG analysis |
| `--laplacian-community-min-size` | `4` | Minimum community size for Laplacian spectra computation |
| `--cooccurrence-min-drugs-per-gene` | `1` | Minimum drugs per gene in the co-occurrence network |
| `--cooccurrence-max-drugs-percentile` | `95.0` | Percentile cutoff for max drugs per gene |
| `--cooccurrence-community-min-size` | `15` | Minimum community size for per-community co-occurrence networks |
| `--similarity-min-degree` | `1` | Minimum node degree for similarity network visualization |

---

## Dependencies

Install required packages before running:

```bash
pip install networkx pandas numpy matplotlib gensim
```

`gensim` is required only for the `node2vec`, `metapath2vec`, and `metapath2vec-standard` algorithms. The `metapath2vec-pp` and `metapath2vec-pp-v0` variants use a custom NumPy-based SGD loop and do **not** require `gensim`.
