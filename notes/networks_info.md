# Drug-target network


# Similarity network
## Parametri
- Tipo di similarità: target (drug che colpiscono lo stesso target)
- Metrica di similarità: cosine
- Threshold di Similarità: Thresold assoluto. Esempi:
    - similarity > 0.3
    - similarity > 0.5
- Peso degli Archi: Rete weighted (consigliato) (il peso = valore della similarità)
- Normalizzazione e Preprocessing: 
    - Rimozione dei faramci con un solo target
    - Gestione di dati mancanti: Se alcuni farmaci hanno pochi target potresti o eliminarli oppure mantenere archi ma pesati meno (da vedere)

Modifiche:
- Layout: ForceAtlas2 (via Fa2 or Gephi)

## Interpretazione del grafico
- Colore dei nodi: grado di ogni nodo, cioè quante connessioni (edges) quel nodo ha all’interno del subgraph della similarity network.
- Edges: rappresentano la similarità tra due farmaci.

Il colore di ogni nodo è direttamente proporzionale al numero di edges ad esso collegati, cioè al suo grado nella similarity network.

## Problema visualizzazione grafico
Come mai ad ogni nuova esecuzione la visualizzazione cambia (anche non in maniera drastica) nonostante il seed sia sempre lo stesso (1):
Quello che può cambiare è solo come il viewer li mostra (zoom, margini, ecc.), quindi possono “sembrare” leggermente diversi anche se in realtà sono la stessa immagine.
Le piccole “differenze” che noti sono molto probabilmente dovute alla visualizzazione (zoom/aspect) o alla percezione, non a un cambiamento effettivo del sottografo.
Un forte indizio di questo è che nodi ed edges della similarity network rimangono sempre dello stesso numero durante diverse esecuzioni con lo stesso seed


# Parametri salvati
results/network_parameters
├─ drug_target_network
│  └─ global_parameters.json
│     • Statistiche aggregate: node_count, edge_count, density,
│       mean/median/min/max/std_weight, component_count, largest_component_size
├─ mid_degree_drug_spotlight
│  ├─ global_parameters.json   • Stesso set di metriche globali (node_count, …, largest_component_size)
│  ├─ filtering.json           • similarity_threshold, nodes_removed, edges_filtered,
│                               original_node_count, retained_node_count, potential_edges (tutti null)
│  ├─ node_parameters.csv      • Colonne: node, degree, weighted_degree,
│                               clustering_coefficient, betweenness_centrality, closeness_centrality
│  └─ edge_parameters.csv      • Colonne: source, target, weight
├─ random_ similarity_snapshot
│  ├─ global_parameters.json   • node_count, edge_count, density,
│                               mean/median/min/max/std_weight, component_count, largest_component_size
│  ├─ filtering.json           • similarity_threshold, nodes_removed, edges_filtered,
│                               original_node_count, retained_node_count, potential_edges
│  ├─ node_parameters.csv      • node, degree, weighted_degree,
│                               clustering_coefficient, betweenness_centrality, closeness_centrality
│  └─ edge_parameters.csv      • source, target, weight (coefficiente di similarità)
└─ similarity_snapshot_by_community
   ├─ global_parameters.json   • node_count, edge_count, density,
   │                             mean/median/min/max/std_weight, component_count, largest_component_size
   ├─ filtering.json           • similarity_threshold, nodes_removed, edges_filtered,
   │                             original_node_count, retained_node_count, potential_edges
   ├─ node_parameters.csv      • node, degree, weighted_degree,
   │                             clustering_coefficient, betweenness_centrality, closeness_centrality
   └─ edge_parameters.csv      • source, target, weight
