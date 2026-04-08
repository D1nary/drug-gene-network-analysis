# Schedule
SIMILARITY NETWORK:
Descrizione metodi +
	filtraggio
	densità
COMMIUNITY:
General parameters:
	Modularity
	size analysis
	
Analysis by community.

# Similarity network
La similarity network è stata costruita utilizzando sia un approccio embedding-based e sia tramite la Jaccard similarity con l'obiettivo di confrontarne analogie e differenze.

## Embedding method
Il primo metodo con cui è stata costruita la similarity network è l'algoritmo di embedding metapath2vec++ (Dong et al., 2017). Esso è simile all'algoritmo node2vec con alcune differenze importanti infatti, metapath2vec++ opera bene su grafi bipartiti, cioè grafi i cui nodi sono di due tipologie diverse come nel caso del dataset  ChG-InterDecagon. La prima differenza sta nel fatto che i random walk non sono liberi ma vincolati a seguire un path specifico  ovvero una sequenza ciclica di tipi di nodo. In questo contesto, il metapath adottato è Drug → Gene → Drug, il quale garantisce che l'embedding di ogni drug sia appreso esclusivamente attraverso la mediazione dei geni con cui interagisce.

Un'ulteriore differenza rispetto a node2vec è che, quando il modello aggiorna la rappresentazione di un drug (skip-gram), considera nel contesto solo gli altri drug e non i geni. I geni presenti nel walk hanno comunque un ruolo fondamentale, poiché fungono da "ponte" per mettere in relazione drug che interagiscono con geni simili, ma non contribuiscono direttamente all'aggiornamento dell'embedding dei drug. 

Con questo algoritmo è stato usato un threshold per la cosine similarity di 0,60. Come avverrà anche per la jaccard similarity, non si è voluto alzare troppo il threhsold in modo da non non lasciare solo connessioni ovvie eliminando quelle potenzialmente interessanti. In altre parole, l'obiettivo è cercare un equilibrio tra la cattura della similarità e l'inclusione di drug con profili genici non fortemente collegati ma comuque aventi collegamenti biologicamente interessanti. 

Siccome nelle seguenti analisi si vorrà analizzare la relativa community network, in cerca di famiglie di farmaci con porofili genetici simili, tenendo un thershold non troppo alto si riduce la frammentazione del dataset del dataset in nodi isolati con la conseguente perdita della possibile creazione di comunità biologicamente interessanti. 

Come verrà analizzato di seguito, la Jaccard similarity individua geni simili considerando esclusivamente il profilo genico. Di conseguenza, due farmaci che non condividono alcun gene avranno Jaccard = 0. Con l'algoritmo di embedding, invece, la similarità è determinata dalla topologia del grafo e non dalla semplice identità dei nodi. Due farmaci risulteranno vicini nello spazio vettoriale se tendono a co-occorrere frequentemente nei random walk. In questo contesto, due farmaci possono agire su geni diversi e condividere comunque una certa similarità, il che è biologicamente plausibile qualora geni diversi appartengano allo stesso pathway.

I seguenti sono i dati relativi al filtering dei nodi dopo aver applicato la cosine similarity:

\begin{table}[H]
\centering
\begin{tabular}{lr}
\toprule
\textbf{Parameter} & \textbf{Value} \\
\midrule
Similarity threshold      & 0.600 \\
Nodes removed             & 26 \\
Edges filtered            & 1\,447\,864 \\
Original node count       & 1\,774 \\
Retained node count       & 1\,748 \\
Potential edges           & 1\,526\,878 \\
\bottomrule
\end{tabular}
\caption{Similarity network parameters (threshold = 0.6).}
\label{tab:similarity_network_parameters_06}
\end{table}

In questo caso i nodi conservati sono il $98.5\%$. A total of 1,447,864 pairs are filtered out, while 79,014 out of 1,526,878 possible pairs are retained. The percentage of connections that manage to exceed the threshold is 5.17%.

Rispetto alla Jaccard, la cosine similarity con soglia 0.6 è molto più permissiva sui nodi (98.5% vs 81.2%). Questo è coerente con quanto detto prima infatti, anche se due nodi non hanno geni simili nel loro profilo, possono comunque essere essere simili a causa della topologia prodotta dall'embedding. Inoltre, la cosine similarity è più selettiva sugli edge. Questo può essere dovuto al threshold più elevato.

---
Nodi rimossi: un nodo viene eliminato se, dopo aver applicato la soglia, non ha nessun arco rimasto — cioè la sua similarità con tutti gli altri nodi è sotto la soglia. È un nodo isolato nella nuova rete.
Archi rimossi: un arco (u, v) viene eliminato semplicemente se sim(u, v) < threshold, indipendentemente da cosa succede agli altri archi.

Perché i numeri differiscono? Perché operano a livelli diversi:

- Un singolo nodo può perdere molti archi prima di diventare isolato.
- Un nodo con 1000 vicini potrebbe perdere 999 archi e restare comunque nella rete (basta uno sopra soglia).
- Un nodo con 2 vicini li perde entrambi e viene rimosso.

---

## Density
La density rappresenta  la frazione di archi esistenti rispetto al massimo possibile. È stata calcolata con la seguente formula:

density = 2 * m / (n * (n - 1))

dove m è il numero totale di archi dopo aver applicato il threshold e n il numero di nodi totale della rete. 

Ha un valore di 0.05 indicando che il $95\%$ dei drug pairs non supera la soglia di simialarità. This confirms a globally sparse network, consistent with the Jaccard-based result.

## Jaccard similarity
DESCRIZIONE IN BLU
In questo caso, si esegue un confronto diretto tra i profili genici di due farmaci. In particolare,

DESCRIZIONE IN GIALLO

Con questo metodo, è stato scelto un threshold di 0.4 per due ragioni principali. L'obiettivo iniziale era quello di preservare una struttura relativamente densa evitando di perdere connessioni deboli le quali potevano rappresentare interessanti legami farmacologici. La seconda ragione sta nel fatto che, già con questo valore, la rete risulta essere molto sconnessa e per questa ragione si è evitato di alzarlo ulteriormente, in modo da non perdere completamente il significato biologico.

I dati relativi al filtering dei nodi sono i seguenti:
TABELLA 5

Since the number of original nodes (original_node_count) is 1774 and the number of
"retained" nodes (retained_node_count) is 1441, it is observed that, after filtering, the
network retains 81.2% of the drugs.

A total of 963,962 pairs are filtered out, while 73,558 out of 1,037,520 possible pairs are retained. The percentage of connections that manage to exceed the threshold is 7.56%

### Density


Il Jaccard method ha prodotto una rete con 1,441 total nodes, 73,558 edges, and a density equal to 0.07, that is, 93% of drug pairs do not exceed the Jaccard threshold. This indicates a globally
sparse network. The global density is strongly influenced by the presence of giant
components, in particular by the one with a size of 359 nodes.

The higher density compared to the embedding-based network (0.05 vs 0.07) reflects the less selective nature of the jaccard threshold.





----------------------------------

# Community network analysis
After constructing the similarity network, the corresponding community network was built using the Louvain method. Below are the Louvain parameters, for both of the methods, saved after the execution of the analysis.

TABELLA

## Embedding method
La Louvain community detection applicata alla rete di similarità coseno sugli embedding prodotti da metapath2vec++ restituisce una struttura comunitaria con 143 comunità.

Il valore di modularità, pari a 0.428, si colloca in una fascia media suggerendo che gli embedding hanno catturato pattern significativi nella rete drug-gene. This suggests that the drugs do not have a random distribution but “organize” into communities, that is, groups characterized by partially shared gene profiles. As will be discussed later, genes belonging to different communities can be connected to each other, producing an interconnected structure.

Per quanto riguarda la distribuzione delle dimensioni delle comunità, si osserva una distribuzione fortemente asimmetrica a destra. La media (12.4) è circa sei volte superiore alla mediana (2.00), il che indica la presenza di poche comunità molto grandi che spostano il valore della media verso l'alto, mentre la maggior parte delle comunità è di piccole dimensioni. La comunità più grande raggiunge 377 nodi, mentre il minimo di 1 segnala la presenza di nodi isolati (singleton), ovvero farmaci o geni che non raggiungono la soglia di similarità coseno sufficiente per essere aggregati ad alcun altro nodo.

Dal punto di vista biologico, la struttura a coda lunga può indicare un nucleo grande di farmaci o geni hub con profili di interazione condivisi insieme a numerosi farmaci con target più specifici che formano piccoli cluster o rimangono isolati.


La distrubuzione delle dimensioni delle comunità può essere osservata anche dai dati nella seguente tabella:
\begin{table}[H]
\centering
\begin{tabular}{lrrrr}
\toprule
\textbf{Size range} 
    & \multicolumn{2}{c}{\textbf{Embedding-based}} 
    & \multicolumn{2}{c}{\textbf{Jaccard-based}} \\
\cmidrule(lr){2-3} \cmidrule(lr){4-5}
    & \textbf{Count} & \textbf{\%} 
    & \textbf{Count} & \textbf{\%} \\
\midrule
Size $< 5$              & 123 & 86.01 & 219 & 83.91 \\
$5 \leq$ Size $\leq 50$ & 11  & 7.69  & 38  & 14.56 \\
$50 <$ Size $\leq 100$  & 4   & 2.80  & 3   & 1.15  \\
Size $> 100$            & 5   & 3.50  & 1   & 0.38  \\
\bottomrule
\end{tabular}
\caption{Distribution of community sizes for the embedding-based (143 communities) and Jaccard-based (261 communities) similarity networks detected via the Louvain algorithm.}
\label{tab:community_size_distribution}
\end{table}

Tra le 5 comunità più grandi abbiamo la community_112 e la community_48 con rispettivamente 191 e 377 elementi. These may correspond to areas associated with widely studied targets or pathways. Another plausible reason for their presence may be the intrinsic nature of the dataset. Indeed, it contains data obtained through variants of the same experimental condition or under different conditions of the same compound.

Since these drugs act on very similar gene sets, more in-depth analyses could be con- ducted to investigate possible therapeutic combinations, including among drugs belonging to closely connected communities.

The high number of small pharmacological modules can indicate rare or highly specific profiles or may represent mechanisms of action that are poorly redundant or scarcely explored and, with further analysis, may serve as interesting starting points for the study of new drugs.

## Jaccard

La Louvain community detection applicata alla rete di similarità ottenuta tramite il jaccard-based method crea una struttura con 261 comunità.

La modurarity è più bassa con un valore di 0.20. Questo indica la presenza di una debole ma reale struttura modulare. Questo risultato è sensibilmente inferiore rispetto al valore ottenuto con il metodo embedding-based (0.428), e può essere ricondotto a diverse ragioni strutturali. Con la Jacacrd similarity, due farmaci sono considerati simili solo se condividono esplicitamente gli stessi geni target. Al contrario, gli embedding prodotti da metapath2vec++ proiettano i nodi in uno spazio in cui vengono preservate co-occorrenze anche indirette, catturando somiglianze tra farmaci che non condividono target diretti ma occupano posizioni topologicamente analoghe nella rete. Ne consegue che la rete embedding-based è intrinsecamente più informativa, in quanto preserva una maggiore quantità di informazione sulla similarità e separazione tra i farmaci, producendo così una modularità più elevata.


Anche in questo caso si ha la presenza di una distribuzione delle size delle community fortemente spostata verso destra con una mediana uguale a 2 indicando che più della metà delle comunity è composta da 2 farmaci o meno. Il valore medio è approssitivamente 5. In questo caso la comunità più grande possiede 359 elementi.

Analizzando la tabella community_size_distribution emerge che le differenze maggiori tra i due metodi per quanto riguarda la distribuzione delle size delle commnity avviene per dimenzioni più elevate. La quota di comunità di dimensione intermedia (5 ≤ size ≤ 50) è considerevolmente più alta nel caso Jaccard (14.56% vs 7.69%), mentre la presenza di comunità molto grandi risulta piuttosto ridotta. Infatti una sola comunità supera i 100 nodi (0.38%), contro le cinque del metodo embedding-based (3.50%).

Di seguito è riportato l'istogramma della distribuzione delle size delle comunità per entrambi i metodi

IMMAGINE

## Intra-community analisis
Per ogni community network creata (embedding-based e jaccard-based) e per ogni community identificata dall'algoritmo di Louvain, viene calcolata density e clustering coefficient del subgrafo indotto dalla community stessa. Il subgrafo indotto è definito come il sottografo del grafo di similarità che include esclusivamente i farmaci appartenenti alla community e gli archi che li connettono tra loro.
### Density
La density, come nel caso precedente, è definita come il rapporto tra il numero di archi presenti nel subgrafo e il numero massimo di archi possibili tra i nodi del subgrafo:
$$\text{density} = \frac{2|E|}{|V|(|V|-1)}$$

dove $|V|$ è il numero di farmaci nella community e $|E|$ è il numero di archi di similarità tra farmaci della stessa community.

In small modules, the number of possible connections is very limited. In fact, it is sufficient for a few nodes to be all connected to each other for the density to be high, often close to 1. This reflects a very strong similarity among the drugs in the group (for example, almost complete sharing of targets), but such values are not very robust from a statistical point of view, because they are strongly influenced by the low number of nodes. For this reason, they were removed from the analysis, considering only communities with size greater than or equal to 5.

#### EMBEDDING
Per quanto riguarda l'embedding-based, sono state calcolate le seguenti density:

GRAFICO SIZE-DENSITY

ISTOGRAMMA DENSITY MEDIA


Il range 5-50 mostra densità media intorno a 0.62 con varianza elevata. In questo range ci sono sia community quasi clique come i nodi 92 e 126 con density = 1, sia community più sparse ($\approx 0.23-0.26$). In questo range convivono quindi trutture qualitativamente diverse. Questo potrebbe essere dovuto al piccolo numero di nodi all'interno della community. Infatti the quadratic scaling of the density denominator fa si che la density di queste community sia molto sensibile alla variazione anche solo di un arco all'interno della community

Come si può osservare dal grafico sopra, nel range 50-100 la density ha una bassa varianza con un valore medio $\approx 0.20$. CI si può aspettare questo comportamento visto the quadratic scaling of the density denominator at larger sizes.




### Clustering coefficient
### degree?








