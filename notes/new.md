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

Le due similarity network ottenute per entrambi i metodi sono visualizzate di seguito:


## Intra-community analisis
Per ciascuna delle due reti di community (embedding-based e Jaccard-based), e per ogni community individuata dall'algoritmo di Louvain, vengono calcolati la densità e il clustering coefficient del sottografo indotto dalla community. Tale sottografo è definito come la porzione del grafo di similarità che include unicamente i farmaci appartenenti alla community e gli archi che li collegano reciprocamente.
### EMBEDDING BASED

Per quanto riguarda l'embedding-based method, sono state calcolate le seguenti density:

TABELLA

#### Density
La density, come nel caso precedente, è definita come il rapporto tra il numero di archi presenti nel subgrafo e il numero massimo di archi possibili tra i nodi del subgrafo:
$$\text{density} = \frac{2|E|}{|V|(|V|-1)}$$

dove $|V|$ è il numero di farmaci nella community e $|E|$ è il numero di archi di similarità tra farmaci della stessa community.

In small modules, the number of possible connections is very limited. In fact, it is sufficient for a few nodes to be all connected to each other for the density to be high, often close to 1. This reflects a very strong similarity among the drugs in the group (for example, almost complete sharing of targets), but such values are not very robust from a statistical point of view, because they are strongly influenced by the low number of nodes. For this reason, they were removed from the analysis, considering only communities with size greater than or equal to 5.

Il range 5-50 mostra densità media intorno a 0.62 con varianza elevata. In questo range ci sono sia community quasi clique come i nodi 92 e 126 con density = 1, sia community più sparse ($\approx 0.23-0.26$). Quindi qui convivono strutture qualitativamente diverse. Questo potrebbe essere dovuto al piccolo numero di nodi all'interno della community. Infatti the quadratic scaling of the density denominator fa si che la density di queste community sia molto sensibile alla variazione anche solo di un arco all'interno della community

Come si può anche osservare dal grafico sottostante, nel range 50-100 la density ha una bassa varianza con un valore medio $\approx 0.20$. Ci si può aspettare questo comportamento visto the quadratic scaling of the density denominator at larger sizes.

Per quanto riguarda le 5 community più grandi, si osservano bassi valori di density, come ci si aspetterebbe per comunità così grandi, tranne in due casi. Questi sono la community più grande ovvero la community_48 (size = 377) e la community_98 (size = 146) con densità rispettivamente di 0.57 e 0.68. Queste sono strutture molto dense e grandi, le quali rappresentano casi di particolare interesse.

---
Osservando la formula della density si capisce perchè per reti grandi ci si aspetta un valore di density basso per community grandi:

Il denominatore cresce quadraticamente con k (size), mentre il numero di archi reali ∣E∣ tende a crescere molto più lentamente — tipicamente in modo lineare o sub-quadratico nelle reti reali. Quindi anche se una community grande ha molti archi in valore assoluto, il rapporto density crolla inevitabilmente.

---

Tali comunità possono rappresentare hub farmacologici, ovvero insiemi di farmaci che condividono un elevato numero di target genici comuni, dando origine a zone della rete caratterizzate da alta densità di connessioni. Oppure, possono raggruppare farmaci che agiscono su pathway biologici centrali e condivisi.

Un'altra possibile ragione per la formazione di queste comunità grandi e dense potrebbe provenire dalla natura intrinseca del dataset. It does not only aggregate information coming from different sources (real screenings, experimental data, computationally generated results), but also information about the same compound obtained from different experimental conditions. This, combined with the fact that the drug IDs within the dataset do not distinguish the experimental context in which the data were produced, may lead to the generation of drugs with identical or nearly identical profiles

GRAFICO SIZE-DENSITY

ISTOGRAMMA DENSITY E CLUSTERING MEDIO

### Clustering coefficient (CC)
Il CC è strutturalmente più stabile rispetto alla density infatti i suoi valori medi rimangono nel range  0.62–0.74. Questo indica che anche quando sono globalmente sparse, le community mantengono una struttura locale triangolare robusta.

Nel range 5-50 il CC è alto e più variabile rispetto agli altri range. Qui abbiamo community molto piccole con struttura quasi-clique (CC $\approx$ 0.90–0.94) accanto a community più grandi e sparse (CC $\approx$ 0.50–0.58)

Negli altri due range, il CC è più omogeneo ed entrambi mostrano valori molto simili con standard deviation bassa ($0.624 \pm 0.049$ , $0.640 \pm 0.058$). Questo indica una struttura locale simile tra i due range nonostante la differenza di dimensione. 

Confrontando i valori di density e CC si osserva che, nel range piccolo (5–50), i due valori medi sono relativamente vicini (density 0.638, CC 0.741), il che è atteso. Nel range medio e grande invece, la density crolla metre il CC rimane comunque $\approx 0.60$. Questo disaccopiamento indica una struttura composta da zone dense separate dalle altre da connessioni deboli. Biologicalmente parlando, questo potrebbe indicare la presenza, all'interno di una stessa comunity, di gruppi di farmaci con target genici fortemente sovrapposti ma con poca sovrapposizione tra gruppi diversi.

---
community dense tendono ad avere anche CC alto, il che è atteso. Come mai?

Immagina una community con density molto alta — quasi tutti i possibili archi esistono. Prendi un nodo *v* qualsiasi con vicini *u* e *w*: poiché la density è alta, è molto probabile che anche l'arco *u*-*w* esista. Il triangolo *v*-*u*-*w* è quindi quasi certamente chiuso.

In altre parole: **se la rete è globalmente densa, ogni coppia di nodi ha alta probabilità di essere connessa, incluse le coppie di vicini comuni**. Il CC non fa altro che misurare localmente questa stessa proprietà — e in una rete densa non può che essere alto.

----

---
Ma se all'interno di una stessa community abbiamo gruppi di farmaci con target genici fortemente sovrapposti al proprio interno, ma con poca sovrapposizione tra gruppi diversi.
Come mai questi farmaci che sono separati all'interno di una stessa comunity sono stati messi nella stessa community?

Louvain non assegna i nodi a una community perché sono tutti densamente connessi tra loro in senso assoluto. Li assegna insieme perché hanno più connessioni tra loro di quante ne avrebbe un grafo casuale con la stessa distribuzione di gradi. È una misura relativa, non assoluta.
Quindi due sottogruppi di farmaci possono finire nella stessa community anche se la loro sovrapposizione reciproca è bassa, purché quella sovrapposizione sia comunque superiore all'atteso per caso.
---

È interessante notare che le due community 98 e 48 precedentemente analizzate mantengono un CC in linea con il resto del range.



### degree?

## Jaccard-based

Per quanto riguarda il Jaccard-based method, sono state calcolate le seguenti density e i seguenti clustering coefficients:

TABELLA

ISTOGRAMMA DENSITY E CLUSTERING MEDIO

Il range 50-100 mostra una density elevata con una deviazione standard più bassa rispetto allo stesso range nel embedding-based method. La presenza di numerose componenti con density = 1 è una caratteristica distintiva di questa rete e potrebbe essere dovuta alla natura della similarità prodotta dal Jaccard method. Infatti, operando direttamente sul profilo genico dei farmaci nel grafo bipartito, tende a produrre clique o near-clique per farmaci che condividono interamente o quasi il proprio target set. Come accade nella rete embedding-based il quadratic scaling del denominatore della density rende questo parametro molto sensibile alla variazione di singoli archi nelle community più piccole, il che spiega l'ampiezza dei valori di denità ($\approx 0.30 - 1.00$).

Il range 50-100 evidenzia la differenza più sostanziale in termini di valore medio e deviazione standard. In questo range sono presenti 3 community 2 delle quali simili in termini di densità e una diversa molto più sparsa.

Nel range size > 100 è presente una sola community con $\text{size} = 359$ e $\text{density} \approx 1.00$. Questo rappresenta un risultato strutturalmente molto diverso da quanto osservato nella rete embedding-based, dove le community più grandi mostravano density nell'intervallo 0.57–0.68. La Jaccard similarity, identifica qui un enorme insieme di farmaci a profilo genico quasi identico. Come osservato nel caso delle community di grandi dimensioni nel metodo embedding-based, le possibili ragioni alla base della formazione di questa community sono possono essere riconducibili alla presenza di farmaci che agiscono su pathway biologici centrali, oppure alla natura intrinseca del dataset.

GRAFICO DENSITY SCATTER

### CC
Nel range 5-50, in particolare per le community clique isolate (degree = 0) di piccole  dimensioni si nota un valore di CC nullo nonostante density unitaria. Escludendo queste utime, si ottiene un sottoinsieme di 16 community il cui CC medio è $\approx 0.661 \pm 0.131$, in linea con i valori osservati nell'analisi embedding-based. Quindi, strutture quasi-clique coesistono con strutture più sparse.


Nel range intermedio è da notare che tutti e tre i valori che non mostrano il fenomeno CC = 0.

La big community clique, ha un CC = 0.903. L'accoppiamento density-CC alto indica che tale community oltre ad essere densa mantiene anche una stuttura traingolare locale.

Precedentemente (embedding-based method) si era ossevato, un pattern di dicaccoppiamento tra density e CC in cui si aveva la presenza di community sparse ma localmente connesse. Nel contesto del Jaccard method, questo non accade. Nei range 50-100 e size>100 la maggiornanza delle community (ad eccezzione della community_188) presenta un accoppiamento tra i due valori. Si può concludere che nella rete Jaccard, farmaci con un profilo di target genici molto simile tendono a formare community in cui sia la densità globale che la connettività locale sono elevate mentre l'embedding genera strutture composte da zone più dense separate da ponti più deboli all'interno della stessa community.

GRAFICO CC SCATTER



## Considerazioni finali
Strutturalmente la jaccard similarity produce community compatte in cui density e CC hanno valori molto vicini fra loro. Al contrario la rete embedding-based genera strutture più eterogenee (globalmente sparse ma localmente coese) caratterizzate dal disaccoppiamento tra density e CC osservato nei range medio-grandi.

Biologicamente parlando, l'embedding-based meethod, lavorando nello spazio delle rappresentazioni vettoriali cattura similarità più latenti e sfumate permettendo di raggruppare farmaci che non condividono necessariamente gli stessi target ma che occupano posizioni analoghe nel grafo bipartito. Questo può essere interpretato come una maggiore capacità di indiivduare analogie indirette. Al contratio, il jaccard method opera direttamente sul profilo genico dei farmaci producendo una similarità di tipo binario la quale può non cogliere sfumature che l'embedding method riescie invece a trovre. In questo contesto l'embedding method potrebbe essere biologicamente più informativo meritando un'analisi più approfondita.

Entrambi i metodi sono soggetti alla natura del dataset che aggrega dati da contesti sperimentali diversi senza distinguere le condizioni di acquisizione. Anche se, il metodo jaccard potrebbe essere più sensibile a questo fattore producendo più comunità clique.

## Betweenness, Closeness and PageRank
Per caratterizzare il ruolo dei singoli farmaci all'interno della rete, sono stati calcolati tre parametri: PageRank, closeness e betweenness a partire dalla rete di similarità embedding-based.

Di seguito un'anticipazione dei valori ottenuti dall'analisi:

\begin{table}[H]
\centering
\begin{tabular}{lrr}
\toprule
\textbf{Metric} & \textbf{Mean} & \textbf{Median} \\
\midrule
PageRank               & $5.72 \times 10^{-4}$ & $5.82 \times 10^{-4}$ \\
Closeness centrality   & $1.741 \times 10^{-1}$ & $1.907 \times 10^{-1}$ \\
Betweenness centrality & $1.805 \times 10^{-3}$ & $1.31 \times 10^{-4}$ \\
\bottomrule
\end{tabular}
\caption{Mean and median values of centrality metrics computed on the drug--drug similarity network.}
\label{tab:centrality_summary}
\end{table}

### pagerank
Con il pagerank, si vuole analizzare la centralità di ciascun farmaco all'interno della drug-drug similarity network. Un pagerank alto indica un farmaco con profilo di interazione genica "hub" condiviso da molti altri. Questo lo potrebbe rendere un potenziale candidato per il drug repurposing.

Di seguito è riportato un istogramma rappresentante la distribuzione dei valori di pagerank. 

HISTOGRAMMA PAGERANK

La distribuzione è quasi gaussiana (con valore medio $\bar{PR} = 5.7 \times 10^{-4} \pm 1.96 \times 10^{-4}$) il che indica che la rete è abbastanza omogena in termini di popolarità di profilo genico con pochi farmaci che spiccano in maniera netta rispetto agli altri. In altre parole, ci sono pochi farmaci con un profilo genico altamente condiviso da altri (profilo genico "hub"). I tre farmaci con il più alto valore di pagerank sono:

\begin{table}[H]
\centering
\begin{tabular}{clr}
\toprule
\textbf{\#} & \textbf{Drug} & \textbf{PageRank} \\
\midrule
1 & Drug\_2451      & 0.001116 \\
2 & Drug\_100054454 & 0.001048 \\
3 & Drug\_151508717 & 0.001032 \\
\bottomrule
\end{tabular}
\caption{Top 3 drugs by PageRank score.}
\label{tab:top3_pagerank}
\end{table}



---
PageRank assegna a ogni nodo un punteggio scalare che riflette la sua importanza strutturale nella rete 
Nel tuo progetto hai una drug-drug similarity network costruita così:

- Nodi: farmaci (proiettati dal grafo bipartito ChG)
- Archi: cosine similarity ≥ 0.4 tra i vettori embedding (MP2Vec-pq o Metapath2Vec++)
- Pesi degli archi: il valore della cosine similarity stessa

In questo contesto PageRank calcola, per ciascun farmaco, quanto è "centrale" nella struttura di similarità globale, ovvero: Un farmaco ha PageRank alto se è simile a molti altri farmaci che a loro volta sono simili a molti farmaci importanti

- Pagerank alto: Farmaco con profilo di interazione genica "hub", condiviso da molti altri ---> Potenziale candidato per drug repurposing (molte relazioni di similarità)
- Pagerank basso: Farmaco con profilo genico periferico o molto specifico ---> Farmaco di nicchia o con target biologico unico

HISTOGRAMMA PAGERANK
L'istogramma mostra quanti farmaci hanno un valore di PageRank che cade in ciascun intervallo. Sull'asse orizzontale ci sono i possibili valori di PageRank (espressi in millesimi), sull'asse verticale il numero di farmaci con quel valore.
In sostanza risponde alla domanda: "quanti farmaci hanno un PageRank basso, quanti medio, quanti alto?"

Nel nostro caso la forma a campana indica che la maggior parte dei farmaci si addensa attorno ai valori centrali — né troppo bassi né troppo alti — con pochissimi nodi agli estremi. È proprio questa forma che ci aveva permesso di concludere che la rete è omogenea, senza farmaci particolarmente dominanti. Ci sono pochi farmaci con pagerank alto ovvero con un profilo dominante. 
---


### Closeness Centrality
Per analizzare quanto un farmaco sia ben posizionato topologicamente all'interno della similarity network si è deciso di calcolare la Closeness Centrality (CC). Un farmaco con CC bassa è strutturalmente lontano dalla maggiorparte della rete, rappresentando un nodo isolato o in un cluster di nicchia. Al contrario un farmaco con CC alta è mediamente separato da pochi passi di similarità da tutto il resto della rete, ovvero un nodo centrale accessibile da molti altri.

La distribuzione della Closeness centrality è la seguente. 

HISTOGRAMMA CLOSENESS CENTRALITY

La distribuzione ha una skewness più marcata con una coda a sinistra la quale indica che la maggiorparte dei farmaci è ben posizionata nella rete mentre un sottoinsieme più piccolo di nodi è più periferico. Il valore medio è $\bar{C} = 0.174, \quad \sigma = 0.074$.

I tre nodi con closeness più alta sono:

\begin{table}[H]
\centering
\begin{tabular}{clr}
\toprule
\textbf{\#} & \textbf{Drug} & \textbf{Closeness} \\
\midrule
1 & Drug\_100001978 & 0.2726 \\
2 & Drug\_9880      & 0.2721 \\
3 & Drug\_28864     & 0.2716 \\
\bottomrule
\end{tabular}
\caption{Top 3 drugs by closeness centrality.}
\label{tab:top3_closeness}
\end{table}

Da notare che molti nodi all'infurori della top 3 per CC hanno valori molto vicini a quello massimo. Questo fatto è confermato anche dalla distribuzione sopra mostrata indicando la presenza di molti nodi ben topologicamente posizionati all'interno della rete. 


---
Mentre PageRank chiede "chi è importante perché è vicino a nodi importanti?", la closeness chiede:

"Da questo farmaco, quanto velocemente posso raggiungere tutti gli altri tramite relazioni di similarità?"

Un farmaco con closeness alta è mediamente vicino a tutti gli altri nella rete: pochi "salti" di similarità lo separano dall'intera popolazione di farmaci. Un farmaco con closeness bassa è periferico, strutturalmente lontano dalla maggior parte della rete.

- Valore alto: Farmaco con profilo genico "intermedio", compatibile con molti cluster diversi --> Raggiungi l'intera rete con pochi salti di similarità
- Valore basso: Farmaco specializzato, simile solo a un sottoinsieme ristretto --> Periferico, magari isolato in un cluster di nicchia

L'interpretazione biologica è sottile: un farmaco con closeness alta non è necessariamente simile a tutti gli altri (quello lo direbbe il degree), ma è ben posizionato topologicamente — fa da intermediario potenziale tra gruppi di farmaci con profili genici diversi.

Differenza chiave
Closeness → posizione globale
(quanto sei vicino a tutti)
Betweenness → ruolo strutturale
(quanto sei importante come passaggio)


---
### Betweenness Centrality
La Betweenness Centrality (BC) è stata calcolata per verificare il ruolo strutturale dei nodi all'interno della similarity network. In particolare, la misura consente di verificare quali farmaci fungono da ponte tra altri farmaci o tra cluster di farmaci. Infatti, un farmaco con elevata betweenness centrality (BC) tende a connettere nodi o gruppi di nodi caratterizzati da profili genici differenti.

---
Cosa misura concretamente
"Quanto spesso questo farmaco si trova sul percorso più breve tra due altri farmaci?"

Un nodo con betweenness alta è un ponte strutturale: rimuoverlo frammenterebbe o allungherebbe significativamente i percorsi nella rete. Non deve avere necessariamente molti vicini — può anche avere degree basso, ma essere l'unico collegamento tra due zone della rete.

- Valore alto: Farmaco che connette cluster con profili genici diversi --> Potenziale candidato per drug repurposing cross-indicazione
- Valore Basso: Farmaco interno a un singolo cluster omogeneo --> Farmaco altamente specifico per un'unica classe terapeutica

Un farmaco con betweenness alta è simile a gruppi di farmaci che tra loro non sono simili — ha un profilo genico "trasversale" che fa da tramite tra domini biologici distinti. Questo è esattamente il tipo di farmaco interessante per ipotesi di repurposing.

Differenza chiave
Closeness → posizione globale
(quanto sei vicino a tutti)
Betweenness → ruolo strutturale
(quanto sei importante come passaggio)
---


La distribuzione che si osserva è la seguente, con valore medio di $\bar{B} = (1.805 \pm 5.281) \times 10^{-3}$

HISTOGRAMMA BETWEENESS

Si nota che la distribuzione è fortemente asimmetrica con coda a destra. Inoltre, il $28\%$ dei nodi ha betweeness uguale a 0 ovvero 489 farmaci non giacciono su alcun shortest path tra altri nodi,  il $72\%$ ha $BC < 0.001$, e solo 13 nodi hanno $BC > 0.04$. Quello che si evince è che pochi farmaci fungono da ponti tra farmaci con diversi profili di similarità. Di seguito i top 3 drug per betweeness:


\begin{table}[H]
\centering
\begin{tabular}{clr}
\toprule
\textbf{\#} & \textbf{Drug} & \textbf{Betweenness} \\
\midrule
1 & Drug\_9880       & 0.06804 \\
2 & Drug\_100002914  & 0.06487 \\
3 & Drug\_2763       & 0.06297 \\
\bottomrule
\end{tabular}
\caption{Top 3 drugs by betweenness centrality.}
\label{tab:top3_betweenness}
\end{table}s


---
DOMANDA
Nel mio contesto, se la CC analizza la vicinanza topologica ad altri nodi della rete e la BC analizza quanto i nodi sono in grado di colleggare altri nodi, stanno misurando sostanzialmente la stessa cosa?

🔴 Differenza concettuale forte (nel tuo caso)
🔹 Closeness Centrality (CC)

Risponde a:

“Quanto questo farmaco è globalmente vicino a tutti gli altri?”

Dipende dalla distanza media nella rete
È una proprietà globale e diffusa
Non richiede di stare su percorsi “critici”

👉 Nel tuo contesto:

Farmaci con CC alta =
ben immersi nella rete di similarità
Tipicamente stanno in regioni dense e ben connesse
🔹 Betweenness Centrality (BC)

Risponde a:

“Quanto questo farmaco è necessario per collegare altri farmaci?”

Dipende dai cammini minimi
È una proprietà strutturale e selettiva
Identifica colli di bottiglia / ponti

👉 Nel tuo contesto:

Farmaci con BC alta =
connettono cluster diversi di profili genici
Anche pochi nodi possono avere valori molto alti (come osservi)
⚖️ Perché NON sono la stessa cosa
1. Puoi avere CC alta ma BC bassa

➡️ Nodo centrale dentro un cluster denso

È vicino a tutti (alta CC)
Ma non è necessario per collegare altri nodi (bassa BC)

👉 Interpretazione biologica:

Farmaco “tipico” di un grande gruppo
Non particolarmente interessante per bridging tra meccanismi diversi
2. Puoi avere BC alta ma CC non altissima

➡️ Nodo ponte tra cluster

Non è vicinissimo a tutti (CC moderata)
Ma è cruciale per collegare regioni diverse (BC alta)

👉 Questo è il caso più interessante per te:

Farmaco con profilo genico trasversale
Forte candidato per drug repurposing cross-cluster
🔬 Collegamento diretto ai tuoi risultati

I tuoi dati confermano perfettamente la differenza:

✔ Closeness
Distribuzione abbastanza concentrata
Molti nodi con valori simili
➡️ rete ben connessa globalmente
✔ Betweenness
Distribuzione molto skewed
28% = 0
pochi nodi molto alti
➡️ pochi veri “ponti strutturali”

👉 Questo è un risultato molto informativo:

La rete è globalmente accessibile (alta CC diffusa)
Ma la comunicazione tra cluster passa da pochi nodi chiave (BC)

---

## Relazione CC - BB

Analizzando le relazioni tra CC e BC è possibile determinare il ruolo di ciascun farmaco all'interno della rete. I ruoli vengono definiti in base alla combinazione di valori alti o bassi di CC e BC, dove la soglia discriminante è fissata al 25° e 75° percentile.

Vengono identificati i seguenti ruoli:
- Hub globali (alta CC - alta BC): Essi si trovano in posizioni "centrali" della rete e, allo stesso tempo, possiedono un ruolo strutturale importante collegando fra loro farmaci con profili genici diversi. Biologicamente possono corrispondere a farmaci con profili di interazione genica molto ampi. Infatti, non solo interagiscono con molti geni, ma quei geni appartengono a pathway funzionali distinti.

- Core di cluster (alta CC - bassa BC): Farmaci ben posizionati nella rete globale, ovvero facili da raggiungere attraverso salti similarità ma che non fungono da ponte tra farmaci o comunità diverse. Probabilmente farmaci appartenenti a famiglie studiate e ben caratterizzate.

- Ponti periferici (bassa cc - alta BC): Questo è il gruppo più ricco. Sono farmaci ai margini della rete (bassa CC) ma indispensabili per l'interconnesione della rete stessa. Come detto prima, sono farmaci che fanno da ponte tra famiglie o farmaci con diverso profilo genico.  Biologimanete parlando, sono molecole con meccanismi d'azione misti, che probabilmetne interferiscono con pathway biologici trasversali.

- Farmaci perifierici (bassa CC e bassa BC): Farmaci alla perifieria della rete che non fungono neanche da ponte tra farmaci i cluster diversi. Possono corrispondere a farmaci altamente selettivi, poco caratterizzati o poco studiati all'interno del dataset

- Nodi misti (al di fuori dei 4 quadranti): Questi farmaci non sono né centrali né periferici in modo netto. Hanno una closeness nella fascia 0.138–0.231 e una betweenness nella facia 0–0.0011. Non fungono da ponti fondamentali per farmaci diversi (bassa BC) ma non sono isolati (CC non trascurabile).

\begin{table}[H]
\centering
\begin{tabular}{llrr}
\toprule
\textbf{Quadrant} & \textbf{Description} & \textbf{N} & \textbf{\%} \\
\midrule
High CC + High BC & Global hubs         & 121  & 6.9\%  \\
High CC + Low BC  & Cluster cores       & 9    & 0.5\%  \\
Low CC + High BC  & Peripheral bridges  & 40   & 2.3\%  \\
Low CC + Low BC   & Isolated peripherals & 302 & 17.3\% \\
Intermediate zone & (outside extreme quartiles) & 1276 & 73.0\% \\
\bottomrule
\end{tabular}
\caption{Distribution of drug nodes across centrality quadrants defined by closeness centrality (CC) and betweenness centrality (BC).}
\label{tab:centrality_quadrants}
\end{table}

---
I dati della similarity network, sono coerenti con il risultato ottenuto? ovvero che i nodi a bassa CC e bassa BC sono il gruppo di dimensione maggiore?

La domanda sulla coerenza con i parametri di rete è molto pertinente — e la risposta è sì, i dati sono coerenti in modo preciso.

Modularity 0.43 e pochissimi bridge
Una modularity di 0.43 indica comunità ben separate, con pochi legami inter-cluster. Questo spiega direttamente perché i bridge periferici siano solo 40: in una rete ad alta modularità i nodi che "attraversano" le comunità sono strutturalmente rari e preziosi. Se la modularity fosse 0.2, ti aspetteresti molti più bridge.

87 componenti e 302 periferici puri
La rete ha 87 componenti connesse di cui solo una gigante (1565 nodi su 1774). I restanti 209 nodi sono in componenti isolate o piccole. Questo spiega gran parte del gruppo LL: nodi disconnessi dal nucleo principale hanno per definizione CC bassissima e BC zero. Non è quindi un risultato biologico ma strutturale — sono farmaci con profili di embedding così distanti da tutti gli altri da non superare la soglia cosine 0.60.

Clustering medio 0.53 e quasi assenza di core di cluster puri (HL = 9)
Un clustering così alto significa che i vicini di un nodo tendono già ad essere connessi tra loro. In una tale struttura, un nodo con alta CC non riesce ad avere BC bassa perché i suoi shortest path vengono ridistribuiti su molti percorsi alternativi — non esiste un nodo che sia "al centro" senza che i suoi vicini si "scavalchino" a vicenda. Questo è esattamente perché HL è quasi vuoto: alta CC e BC nulla coesistono solo in reti con cluster molto isolati e densi, che qui la modularity moderata non permette completamente.

Diameter 17 e avg_path_length ~4.9
Il diametro elevato (17) rispetto al path length medio (4.9) indica una distribuzione asimmetrica dei cammini: la maggior parte delle coppie di nodi è relativamente vicina, ma esiste una coda lunga di percorsi che attraversano la rete da periferia a periferia. Questa struttura è coerente con i bridge periferici: sono proprio i nodi che, pur stando ai margini, accorciano quei percorsi lunghi tra cluster distanti — toglierli aumenterebbe il diametro in modo sproporzionato.

In sintesi, i parametri di rete non solo sono coerenti con la distribuzione dei quadranti, ma la spiegano meccanicisticamente. Il quadro che emerge è quello di una rete con comunità ben formate, bordi netti, e un nucleo connesso ma non dominante — esattamente il tipo di struttura in cui i bridge periferici sono rari, critici, e biologicamente i più interessanti da investigare.

SCATTER PLOT bc - cc



Scatter plot of closeness (CC) versus betweenness centrality (BC) for \(n = 1{,}748\) drug nodes in the embedding-based similarity network. The $y$-axis uses a $\log_{10}$ scale. Dashed lines mark the 25th/75th percentile thresholds of CC and the 75th percentile of BC

\caption{Distribution of PageRank scores across the 1{,}748 drug nodes in the embedding-based similarity network.}

\caption{Distribution of closeness centrality across the 1{,}748 drug nodes in the embedding-based similarity network.}

\caption{Distribution of betweenness centrality across the 1{,}748 drug nodes in the embedding-based similarity network.}
