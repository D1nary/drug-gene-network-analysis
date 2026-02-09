# Descrizoine del dataset
Il report analizza il dataset “Chemical-gene interaction network” (ID 10016-ChG-InterDecagon) che può essere trovato al seguente link: https://snap.stanford.edu/biodata/datasets/10016/10016-ChG-InterDecagon.html.
Esso rappresenta un network biologico in cui i nodi sono farmaci/composti chimici e geni/proteine. Gli archi della rete rappresentano le interazioni biologiche tra questi elementi. Tali interazioni sono associazioni funzionali o biomediche come il legame di un composto chimico ad una proteina target, l'attivazione o l'inibizione di un gene effetti osservati sperimentalmente e predizioni computazionali. Per le successive analisi, è importante evidenziare che il dataset non è limitato ad un singolo contesto ma aggrega anche dati provenienti da condizioni sperimentali differenti. Di seguito si trova una tabella con tutte le informazioni specifiche riguardanti il dataset

| Dataset statistic                         | Valore    |
|------------------------------------------|----------:|
| Nodes                                    | 9 569     |
| Drug nodes                               | 1 774     |
| Gene nodes                               | 7 795     |
| Edges                                    | 131 034   |
| Nodes in largest SCC                     | 9 538     |
| Fraction of nodes in largest SCC         | 1.000000  |
| Edges in largest SCC                     | 131 001   |
| Fraction of edges in largest SCC         | 0.999748  |
| Diameter (longest shortest path)         | 8         |
| 90-percentile effective diameter         | 3.864298  |
Label: Dataset statistic. Fonte: https://snap.stanford.edu/biodata/datasets/10016/10016-ChG-InterDecagon.html

## Preprocessing
Prima di ogni analisi e della creazione delle reti, è stato svolto un preprocessing dei dati nel dataset. Dopo una prima pulizia sintattica, ovvero rimozione di commenti, spazi vuoti e righe incomplete, è stata effettuata normalizzazione in modo da ottenere geni e farmaci ben identificati e facilmente utlizzabili per le successive analisi. In particolare gli identificativi sono stati ripuliti rimuovendo spazi vuoti e prefissi con una successiva conversione dell'identificativo da stringa a valore numerico. Inoltre, sono state rimosse tutte le righe duplicate in modo da non avere nodi "artificiali" nel dataset.

Dopo il prerocessing gli identificativi dei farmaci e gli ID dei geni hanno la seguente froma:
       
| Drug       | Gene   |
|------------|--------|
| 60752      | 3757.0 |
| 6918155    | 2908.0 |
| 103052762  | 3359.0 |
| 23668479   | 1230.0 |
| 28864      | 1269.0 |
Label: Esempio ID farmaci e geni dopo il preprocessing

# Methodology
## Grafo bipartito
Per una rappresentazione visuale della rete è stato creato un grafo bipartito drug-gene. Per una visualizzazione più chiara, non sono stati utilizzati tutti i nodi ma solo quelli con un grado compreso tra 5 e 15 imponendo un numero massimo di nodi drug a 50, un numero totale di nodi a 200. Il grafo ottenuto possiede:


| nodes | edges |
|---------------|
|  137  |  403  |
Label: Numero di nodi e edges nel grafo sopra rappresentato 

La visializzaazione del grafo è la seguente:

IMMAGINE GRAFO BIPARTITO
Label: grafo bipartito drug gene di esempio 

## Similarity network
È stata costruita una similarity network drug–drug utilizzando la Jaccard similarity, applicando inizialmente un threshold pari a 0.3 e successivamente uno pari a 0.4. I risultati ottenuti con le due soglie sono stati poi confrontati tra loro.

La Jaccard similarity misura il grado di sovrapposizione relativa tra due farmaci ed è definita come:
$$
J(A,B) = \frac{|A \cap B|}{|A \cup B|}
$$

Essa indica che due farmaci risultano simili solo se condividono un numero significativo di target e, allo stesso tempo, presentano pochi target differenti, cioè se hanno un profilo di bersagli complessivamente simile.

È stata scelta la Jaccard similarity perché opera nativamente su insiemi, penalizzando farmaci con un elevato numero di target non condivisi. Questo la rende particolarmente adatta all’analisi della similarità tra profili di target tipici del dataset analizzato

Ogni nodo del grafo della similarity network rappresenta un farmaco, caratterizzato da un profilo bersaglio, ovvero da un vettore contenente tutti i geni (target) con cui il farmaco interagisce. Un esempio di rappresentazione di un nodo è il seguente:

```bash
{
  "bipartite": "drug",
  "original_id": 12345,
  "targets": [1017, 1956, 7422]
}
```

### Threshold 0.3
I parametri relativi al filtraggio della rete di similarità con threshold 0,3 sono i seguenti:

| Parametro            | Valore    |
| -------------------- | --------- |
| similarity_threshold | 0.3       |
| nodes_removed        | 333       |
| edges_filtered       | 961,999   |
| original_node_count  | 1,774     |
| retained_node_count  | 1,441     |
| potential_edges      | 1,037,520 |
label: Parametri similarity network (threshold = 0,3)


Questo valore di threshold indica che un arco (drug–drug) è presente nella rete solo se la similarità tra due farmaci è ≥ 0.3. Essa è una soglia moderata che mantiene connessioni con similarità medio–bassa, quindi preserva una rete relativamente densa rispetto a soglie più stringenti. Genera una rete connessa ma non troppo, utile per il community detection riducendp il rumore ed eliminando similarità assolutamente deboli. 

Siccome il numero di nodi originali (original_node_count) è 1774 e il numero di nodi "sopreavvisuti" (retained_node_count) è 1441, si osserva che, dopo il filtraggio, la rete conserva l'81,2% dei farmaci.

Il parametro potential_edges rappresenta il totale delle coppie drug-drug possibili prima del filtraggio mentre edges_filtered indica quante coppie sono state eliminate perché avevano similarità < 0.3. Sottraendo questi due valori, si ottiene che la rete finale ha 85310 archi. Quindi solo circa l’8.22% delle possibili connessioni supera la soglia, indicando una rete relativamente sparsa. Questa sparsità è coerente con reti di similarità farmacologica infatti, generalmente pochi farmaci sono veramente simili.

### Threshold 0.4
I parametri relativi al filtraggio della rete di similarità con threshold impostato a 0,4 sono i seguenti:

| Parametro            | Valore    |
| -------------------- | --------- |
| similarity_threshold | 0.4       |
| nodes_removed        | 333       |
| edges_filtered       | 963,962   |
| original_node_count  | 1,774     |
| retained_node_count  | 1,441     |
| potential_edges      | 1,037,520 |

Con un threshold di 0.4 la soglia è più stringente mantenendo solo interazioni più forti rispetto a prima. Vengono filtrati 963,962 mantenendone 73,558 su 1,037,520 possibili coppie. La percentuale delle connessioni che riescono a superare il threshold è del 7.56%

Con tale threshold la rete perde $85.310 - 78.403 = 6.907 $ archi rispetto al caso precedente. Questo numero rappresenta una diminuzione del 
$$
\frac{6.907}{85.310} \approx 8.1 \%
$$
del numero di connessioni rispetto al threshold precedente.

Il numero di nodi rimossi rimane invariato a 333. Ciò indica che questi farmaci non presentano valori di similarità pari o superiori a 0.4 con nessun altro farmaco del dataset, risultando quindi isolati nella similarity network già al threshold più precedente.

Uno degli scopi di questa analisi, è quello di costruire una community network per ricercare comunità di farmaci con meccanismài d'azione simili. In questo contesto, non avrebbe senso alzare troppo del threshold potrebbe avere l'effetto di tenere solo le relazioni “ovvie” e perdere quelle deboli ma biologicamente interessanti. Infatti anche similarità più deboli rappresentano un caso interessante poichè, due farmaci che condividono pochi target possono essere candidati per drug repurposing. Un farmaco con profilo di target solo parzialmente simile può comunque avere effetti collaterali simili oppure opportunità di combinazione terapeutica. Imporre soglie troppo rigide rischia di semplificare eccessivamente una realtà biologica che è invece caratterizzata da relazioni sfumate e parziali.

Le seguenti analisi sono state eseguite considerando un threshold di 0.4.


### Similarity network global data

La rete prodotta possiede 1441 nodi totali, 73558 archi e una densità pari a 0.07. La densità è stata calcolata con:
\[
\text{density} = \frac{2E}{N(N - 1)}
\]

Per reti biologiche, il valore di density è un valore moderatamente alto. Rappresenta una rete connessa ma non prossima alla saturazione. Com'è possibile osservare dalla segunete rappresentazione della similarity network e come verrà discusso dopo, la density globale è fortemente influenzata dalla presenza di giant components e in particolare da quella avente dimensione di 359 nodi. 

Di seguito è mostrata una rappresentazione di 500 nodi scelti casualmente dalla similarity network

IMMAGINE SIMILARITY
label: Visaulizzation of 500 casual nodes from the similarity network


# Community network analysis
Dopo aver costruito la similarity nework, è stata costruita la relativa community network attraverso il louvian method.
Di seguito sono riportati i parametri di Louvian salvati nel file louvian_parameters.json dopo l'esecuzione dell'analisi.

| Parametro              | Valore |
|------------------------|--------|
| method                 | louvain |
| resolution             | 1.000  |
| modularity             | 0.204  |
| min_community_size     | 1      |
| max_community_size     | 359    |
| mean_community_size    | 5.521  |
| median_community_size  | 2.000  |

## Modularity
Una modularity di circa 0.20 indica la presenza di una struttura modulare debole ma reale. Questo ci suggerisce che i farmaci non possiedono una distruìibuzione casua ma si "organizzano" in community, ovvero gruppi caratterizzati da profili di geni parzialmente condivisi. 

La natuara dei sistemi biologici e faramcologici è complessa. Infatti in essi farmaci diversi possono avere effetti parzialmente simili, agendo su uno stesso set di geni condividendo pathway o target, pur appartendendo a a diverse classi funzionali. Questo comportamento è presente anche nel dataset analizzato in cui nodi appartenenti a community diverse sono connessi suggerendo la presenza di una struttura interconnessa. In questo contesto, tali farmaci possono avere potenziali effetti collaterali comuni e  possibili opportunità di combinazione terapeutica o drug repurposing.

## Other parameters
- Median:
Il valore mediano delle comunità risulta essere 2. Questo significa che più della metà delle comunità ha 2 farmaci o meno. Questo suggerisce gruppi di farmaci molto simili fra di loro i quali possono essere varianti strutturali di uno stesso composto oppure farmaci che condividono uno o pochissimi target molto specifici. In altre parole, queste comunità rappresentano farmaci di nicchia con bersagli rari oppure possono essere outliers interessanti per riposizionamento. Per esempio: se un singleton si collega debolmente a una grande comunità, potrebbe condividere qualche pathway con farmaci di un’altra indicazione.

- Mean, max community 
La media è rappresentata da una comunità di circa 5 farmaci. Questo valore insieme a quello della mediana, indicano una distribuzione fortemente sbilanciata con pochi cluster molto grandi e tanti piccoli. 

ISTOGRAMMA
label: The histogram rapresents the sizes of all the communities

Il valore massimo è rappresentato da una comunità con 359 elementi. Questa comunità verrà approfondita meglio nelle sezioni sucessive.


## Inforamazioni biologiche ricrea il file colab per la stampa del graficoavabili da questi valori

Anche in assenza di un’analisi dettagliata della composizione interna di ciascuna comunità, i parametri ottenuti tramite il metodo di Louvain consentono di trarre alcune considerazioni biologiche rilevanti. In primo luogo, la presenza di una modularità diversa da zero indica che la rete di farmaci non è organizzata in modo casuale, ma presenta una struttura modulare, suggerendo l’esistenza di gruppi di farmaci che condividono meccanismi e bersagli potenzialmente simili.

La presenza di farmaci appartenenti a comunità diverse suggerisce la possibilità che farmaci appartenenti a classi terapeutiche diverse, possano, poichè agenti su stessi geni, indurre degli effetti collaterali simili. Questo è interessante dal punto di vista del riposizionamento farmacologico. 



## Communities analisys
Le singole comunità vengono analizzate nel file community_parameters.csv. Il file è organizzato in colonne ciascuna con una caratteristica della comunità:
- community id: Identificativo della comunità
- size: Dimensione della community
- degree: Numero di arhci del nodo
- weighted degree: Somma dei pesi degli archi incidenti su un nodo
- clustering coefficient: Grado di chiusura locale delle comunità

Nella seguente tabella sono stati riportati i dati ottenuti dall'analisi solo delle comunità con size maggiore o uguale a 10.


community_id  size degree weighted_degree clustering_coefficient density
------------- ---- ------ --------------- ---------------------- -------
Community_5   359  0      0.000           0.903                  1.000  
Community_13  18   0      0.000           0.000                  0.791  
Community_20  22   0      0.000           0.593                  0.792  
Community_21  82   0      0.000           0.631                  0.778  
Community_22  83   1      12.990          0.851                  0.913  
Community_53  10   0      0.000           0.000                  0.733  
Community_54  22   0      0.000           0.574                  0.745  
Community_58  11   0      0.000           0.000                  1.000  
Community_59  10   0      0.000           0.000                  0.889  
Community_78  32   1      12.990          0.653                  0.905  
Community_81  20   0      0.000           0.964                  1.000  
Community_97  19   0      0.000           0.000                  0.754  
Community_100 13   0      0.000           0.000                  0.346  
Community_109 13   0      0.000           0.000                  1.000  
Community_127 14   0      0.000           0.000                  0.604  
Community_135 16   0      0.000           0.000                  0.425  
Community_187 10   0      0.000           0.000                  1.000  
Community_188 72   0      0.000           0.616                  0.324  
Community_192 39   0      0.000           0.706                  0.601  
Community_203 23   0      0.000           0.513                  0.482  


### Size

La distribuzione delle size delle community indica la presenza di un elevato numero di piccoli moduli farmacologici. Questi farmaci potrebbero avere profili rari o altamente specifici. Questi elementi del network, inoltre possono rappresentare meccanismi d'azione poco ridondanti o scarsamente esplorati e, con ulteriori analisi, possono essere interessanti punti di partenza per  l'analisi di nuovi farmaci.

Si nota la presenza di comunità molto grandi. Queste possono essere aree associate a target o pathway ampiamente studiati. Un'altra ragione plausibile per la presenza di queste ultime può essere, come vedremo più avanti, la natura intrinseca del dataset. Infatti esso contiene dat ottenuti attraverso varianti di una stessa condizione sperimentale o in condizioni differenti dello stesso composto. Siccome questi farmaci agiscono su set di geni molto simili, si potrebbero condurre delle analisi più approfondite per la ricerca di possibili combinazioni terapeutiche, anceh tra farmaci di communità strettamente connesse.

range           count percent
--------------- ----- -------
size < 5        219   83.91  
5 ≤ size ≤ 50   38    14.56  
50 < size ≤ 100 3     1.15   
size > 100      1     0.38   
label: Percentuali delle size delle community

### Density
La densità misura quanto i farmaci all’interno di una comunità siano effettivamente simili tra loro poichè queste ultime, benchè individuate dal metodo Louvian, possono essere tenute insieme da similarità parziali. In questo contensto la density aiuta a distinguere quali comunità sono realtente omogenee da quelle i cui componenti condividono solo alcuni bersagli.

$$
\text{density} = \frac{E_{\text{int}}}{\frac{n(n-1)}{2}}
$$
dove:

- density: frazione di connessioni interne presenti rispetto al massimo possibile nella comunità.
- E_int: numero di archi interni alla comunità (solo tra nodi della stessa comunità).
- n: numero di nodi nella comunità (size).
- n(n-1)/2: numero massimo di archi possibili in una comunità non orientata senza self‑loop.

Nei moduli piccoli, il numero di connessioni possibili è molto limitato, infatti basta che pochi nodi siano tutti connessi tra loro perché la density risulti elevata, spesso prossima a 1. Questo riflette una similarità molto forte tra i farmaci del gruppo (ad esempio condivisione quasi completa dei target), ma tali valori sono poco robusti dal punto di vista statistico, perché fortemente influenzati dal basso numero di nodi. Per questo sono stati rimossi dall'analisi considerando solo comunità con size maggiore o uguale di 10.

Dai dati riportati si osserva che comunità medio-grandi, presentano valori di densità variabile. Analizzando le comunità  con size compresa tra 39 e 83 si osservano densità molto diverse. Ovvero comprese tra 0.324 e 0.913. 
Le cause di questo comportamento possono essere sia matematiche che biologiche. 

Le cause matematiche sono da attribuirsi a come viene calcolata la density (formula precedente). Dal momento che il numero di possibili connessioni interne aumenta quadraticamente con la dimensione della comunità, anche piccole variazioni nei profili dei nodi tendono ad accumularsi in un numero crescente di archi mancanti. Questo effetto è amplificato nelle comunità più grandi, dove l’elevato numero di confronti tra profili rende più probabile l’emergere di differenze, riflettendosi in una maggiore variabilità dei valori di densità.

Inoltre, la presenza o l'assenza di un arco dipende da una soglia (0.4). Quindi, se due nodi sono appena sotto la soglia l'arco nella similarity network è assente mentre, se la similarity è appena sopra la soglia, l'arco è presente. Questo introduce variabilità artificiale e maggiore dispersione dei peorfili di densità. 

Da un punto di vista biologico, comunità grandi spesso aggregano famiglie farmacologiche ampie, pathway complessi o target parzialmente sovrapposti creando, di conseguenza cluster di più grandi dimensioni ma non completamente connessi.

Per le comunità di dimensione più contenuta (size compresa tra 10 e 32), si osserva in generale una densità elevata, con la maggior parte dei valori compresa tra ~0.7 e ~1.0. Sono tuttavia presenti alcune eccezioni con densità più moderata o bassa, come Community_135 (0.425), Community_203 (0.482) e Community_100 (0.346).

#### Comunità clique

Nella rete sono presenti due comunità quasi-clique. La prima, è la comunità con il maggior numero di elementi. Essa presenta le seguenti caratteristiche:

- density: 0.999 
- size: 359
- unique_profiles: 120
- shared_profiles: 81
- most_frequent_profile: 
   - drug_count: 125
   - gene_count: 158
   
Quindi più di un terzo dei componenti della comunità condivide lo stesso esatto insieme di geni. Questo implica che la Jaccard similarity è J = 1 per tutte le coppie di farmaci che condividono lo stesso profilo e J ≈ 1 per le coppie di nodi che differiscono tra loro per pochi geni. Di conseguenza il sottografo della comunità, ovvero il grafo contenente i nodi della comunità e tutti gli archi di similarità compresi, risulta essere quasi completamente connesso producendo density quasi unitaria.

Uno delle possibili cause di questo valore alto di density per una comunity cosi grande, può essere che la Jaccard diventa molto permissiva quando i set sono grandi e l'intersezione tra due set è molto ampia. Infatti avendo set di grandi dimensioni, anche con decine di geni diversi, la similarità resta alta. Avendo poi un threshold sufficientemente permissivo, come nel nostro caso, il valore di $E_{\text{int}} \approx \frac{n(n-1)}{2}$ con conseguente density $\approx 1$.

Un'ulteriore causa del valore quasi unitario di questa community, come accennato precedentemente, potrebbe provenire dalla natura intrinseca del dataset. Infatti esso, aggregando screening diversi producendo il risultato che, farmaci testati negli stessi screening sugli stessi pannelli genici con risultati simili diventano vettorialmente indistinguibili. Questo produce sottografi altamente densi nella rete di similarità.

Un'ulteriore causa del valore quasi unitario della density di questa community, può provenire dalla natura intrinseca del dataset. Esso non aggrega solo informazioni provenienti da fonti diverse (screening reali, sperimentali, prodotti computazionalmente) ma anche informazioni di uno stesso composto provenienti da screening diversi. Questo, unito al fatto che gli ID dei farmaci all'interno del dataset non distinguono il contesto sperimentale in cui sono stati prodotti i dati, può far si di generare farmaci con profili identici o quasi.

Senza nessuna ulteriore analisi, questa comunità rappresenta una famiglia di farmaci con meccanismo d'azione quali identico oppure il targetting di un grande pathway.

Sono presenti altre comunità con density unitaria per le quali può essere applicato lo stesso ragionamento appena fatto. Unito a questo fatto c'è da considerare che sono presenti clique banali come accade per la comunità 20. Essa continene 20 farmaci di cui 19 con lo stesso identico profilo composto da soli due farmaci. Questa rappresenta un clique più banale del precedente poichè il profilo è piccolissimo ed è replicato perfettamente.


### Clustering coefficient
Per ciascuna comunity, è stato calcolato il clustering coefficient:

community_id   size  density  clustering_coefficient
Community_5    359   1.000    0.903
Community_13   18    0.791    0.746
Community_20   22    0.792    0.593
Community_21   82    0.778    0.631
Community_22   83    0.913    0.851
Community_53   10    0.733    0.660
Community_54   22    0.745    0.574
Community_58   11    1.000    0.831
Community_59   10    0.889    0.629
Community_78   32    0.905    0.653
Community_81   20    1.000    0.964
Community_97   19    0.754    0.553
Community_100  13    0.346    0.442
Community_109  13    1.000    0.787
Community_127  14    0.604    0.671
Community_135  16    0.425    0.596
Community_187  10    1.000    0.660
Community_188  72    0.324    0.616
Community_192  39    0.601    0.706
Community_203  23    0.482    0.513


Si osservano diversi regimi di clustering coefficinet riassumibili in due: Clustering coefficient alto o moderatamente alto e clustering coefficient uguale a 0. Considerando il caso dei valori moderatamente elevati, si può dire che tali comunità  non presentano una struttura lineare o “a catena”, in cui i farmaci risultano simili solo a pochi vicini immediati, ma piuttosto costituiscono moduli fortemente coesi. Cioè, se un farmaco è simile ad altri farmaci all'interno della comuintà, è molto porbabile che, per come è stata definita la similarità e per la soglia utilizzata, esso sia siamile ad altri componenti della stessa.

Le community 13, 20, 21, 22, 53, 54, 59, 78, 97 presentano valori compresi tra::
| Size  | Density     | Clustering  |
| ----- | ----------- | ----------- |
| 10–83 | 0.733–0.913 | 0.553–0.851 |

quindi con densità alta ma minore di 1 e clustering moderatamente alto o alto. Da questi dati si può dire che non tutti i farmaci sono simili agli altri ma esistono sottogruppi molto coerenti all'interno della comunità. Biologicamente, è possibile la presenza di famiglie farmacologiche in cui esiste un nucleo di target comuni ed in cui ogni farmaco può introdurre variazioni marginali sul profilo.

Sono presenti inoltre, comunità con bassa density e clustering moderato

Tra le comunità con una density più bassa abbiamo:
| Community     | Size | Density | Clustering |
| ------------- | ---- | ------- | ---------- |
| Community_100 | 13   | 0.346   | 0.442      |
| Community_135 | 16   | 0.425   | 0.596      |
| Community_188 | 72   | 0.324   | 0.616      |
| Community_203 | 23   | 0.482   | 0.513      |

Queste mostrano poche connessioni globali ma connessioni locali ben strutturate evidenziando il fatto che la comunità non è un blocco compatto ma un insieme di cluster locali collegati indirettamente.

Un'ultimo caso che si può notare dalla tabella precedente sono comunità con density $\approx 1$ e clustering alto. 

| Community     | Size | Density | Clustering |
| ------------- | ---- | ------- | ---------- |
| Community_5   | 359  | 1.000   | 0.903      |
| Community_58  | 11   | 1.000   | 0.831      |
| Community_81  | 20   | 1.000   | 0.964      |
| Community_109 | 13   | 1.000   | 0.787      |
| Community_187 | 10   | 1.000   | 0.660      |

In queste comunità il clustering coefficient non aggiunge nuova informazione rispetto alla density, ma rafforza l’evidenza di omogeneità estrema.

## Weighted degree
Per descrivere il livello complessivo di interazione di una comunità con le altre è stato calcolato, il weighed degree. Esso è calcolato come la somma dei pesi degli archi inter-comunità (similarità Jaccard). 

Solo le comunità 23 e 78 risultano essere connesse tra di loro con un $\text{weighted degree} = 12.99$.

Siccome, un arco tra due comunità esiste solo se esistono farmaci appartenenti a comunità diverse ma connessi e che quindi condividono una porzione significativa di target genici (Jaccard maggiore uguale a 0.4), le due comunità rappresentano moduli farmacologici distinti ma non indipendenti. I farmaci "ponte" che connettono le due comunità possono agire su pathway diversi ma comunque interconnessi attraverso geni chiave.

Il valore di weighted degree legato alla connessione di queste due comunità è di 12.99.

# Co-occurence network
Per analizzare la frequenza con cui coppie di geni compaiono insieme nei profili dei farmaci, è sato deciso di costruire una serie di cooccurence network.
A causa dell'elevato costo computazionale relativo al calcolo di una co-occurence network globale, si è scelto di creare ed analizzare co-occurence network per ciascuna comunità precedentemente identificata attraverso la community analisys. Ci si è concentrati sulle comunità con una $\text{size} \geq 15$. Oltre al filtro appena citato, ne è stato applicato un altro in modo da rendere i dati di più facile interpretazione rimuovendo rumore e strutture quasi clique. In particolare sono stati rimossi dalle community geni super-frequenti ovvero i geni il cui numero di farmaci associati è maggiore del 95° percentile della distribuzione “farmaci per gene”. 

Di seguito sono riportati i parametri delle cooccurence network trovate per ciascuna community:

| Community ID | n_nodes | n_edges | density | component_count | giant_component_size | global_clustering_coefficient | community_size |
| ------------ | ------- | ------- | ------- | --------------- | -------------------- | ----------------------------- | -------------- |
| Community_5 | 357 | 47957 | 0.755 | 1 | 357 | 0.819 | 359 |
| Community_13 | 166 | 13626 | 0.995 | 1 | 166 | 0.997 | 18 |
| Community_20 | 17 | 136 | 1 | 1 | 17 | 1 | 22 |
| Community_21 | 26 | 280 | 0.862 | 1 | 26 | 0.91 | 82 |
| Community_22 | 77 | 2563 | 0.876 | 1 | 77 | 0.914 | 83 |
| Community_54 | 4 | 6 | 1 | 1 | 4 | 1 | 22 |
| Community_78 | 89 | 3394 | 0.867 | 1 | 89 | 0.921 | 32 |
| Community_81 | 2 | 1 | 1 | 1 | 2 | 0 | 20 |
| Community_97 | 215 | 22151 | 0.963 | 1 | 215 | 0.973 | 19 |
| Community_135 | 86 | 2863 | 0.783 | 1 | 86 | 0.883 | 16 |
| Community_188 | 82 | 2594 | 0.781 | 1 | 82 | 0.865 | 72 |
| Community_192 | 34 | 427 | 0.761 | 1 | 34 | 0.872 | 39 |
| Community_203 | 5 | 10 | 1 | 1 | 5 | 1 | 23 |
label: cooccurence network data

## component count e giant component size
Tutte le co-occurrence networks gene-gene hanno component count = 1. Questo vale indipendentemente dal numero di nodi (da 2 a 357), dalla densità (da ~0.75 a 1), dalla dimensione della community e dal clustering coefficient. Questo sta ad indicare che non esistono sottogruppi genetici isolati all’interno della community ovvero che la rete gene-gene è completamente connessa. In questo contesto, ogni gene è raggiungibile da ogni altro gene tramite almeno un cammino. 

Sempre tenendo conto dei limiti imopsti dall'analisi e dalla natura instrinseca del dataset (discussa prima), da un punto di vista biologico, un valore di component count di 1 suggerische che i geni all'interno di ciascuna community partecipano ad un sistema funzionale interconnesso. Non emergono sottoinsiemi genetici separati maun unico blocco funzionale.

Un dettaglio che si può notare dai risultati ottenuti è che in tutte le co-occurence network, la giannt component coincide con le dimensioni dell'intera rete. In questo scenario, il concetto di giant component perde la sua usuale accezione di “sottostruttura dominante” e diventa equivalente alla rete stessa.

## Density e global clustering coefficient
### density
Siccome la connettività globale della rete è già garantita (component conunt = 1) per costruzione, in questo caso la density ci dice quanto è ridondante la condivisione dei farmaci tra i geni della comunity.

La density, calcolata con la formula vista in precedenza, presenta valori elevati in quasi tutti i casi ovvero si ha assenza di reti sparse (density < 0.5). Possiamo identificare tre regimi diversi benchè tutti con valori elevati:

Nel primo regime di density, si hanno valori [0.95,1.00]. In questo caso siamo in assenza quasi totale di eterogeneità dei profili in cui tutte o quasi tutte le coppie di geni sono collegate (reti clique o quasi). Si ha una forte ridondanza farmacologica in cui i geni sono colpiti, nella maggioranza dei casi, dagli stessi farmaci (forte ridondanza farmacologica).

Nel seconod, la density compresa tra [0.75,0.90). Siamo in presenza di reti ancora pienamente connesse però con un numero più o meno significativo di archi mancanti. In queste strutture, meno clique e più articolate delle precedenti, ci possono essere geni coinvolti in processi distinti ma indirettamente interconnessi

### clustering coefficient
Anche nel caso del clustering coefficient, siamo in presenza di valori molot alti. Possiamo distinguere due regimi: valori massimi e valori alti ma < 1.
Nel primo caso siamo in presenza di strutture totalmente ridondanti geneticamente in cui ogni tripla forma forma un traingolo. 
Nel secondo caso (≈ 0.88 – 0.95) abbiamo strutture non completamente chiuse con moduli genetici più articolati rispetto al caso precedente con la presenza di geni ponte che collegano sottoinsiemi funzionali diversi. 

È presente una community con clustering coefficient = 0. Questa è una rete composta da 2 geni quindi assenza di trangoli per definizione.

## Wheight distribution and sparsity

Community	n_nodes	median	mean	max	weight_eq_1_pct	weight_eq_2_pct
Community_5	880	2,000	32,209	359	6,432	64,565
Community_13	453	2,000	3,016	18	14,154	65,961
Community_20	38	2,000	3,701	20	10,598	71,467
Community_21	52	2,000	3,812	78	31,768	35,368
Community_22	181	2,000	8,426	83	9,834	62,845
Community_54	14	2,000	2,588	21	0,000	94,118
Community_78	220	2,000	4,178	31	29,849	39,737
Community_81	3	1,000	7,333	20	66,667	0,000
Community_97	390	2,000	4,371	19	11,056	49,854
Community_135	175	2,000	2,573	9	46,854	30,391
Community_188	130	2,000	5,874	33	6,186	45,450
Community_192	63	2,000	5,691	31	1,700	71,000
Community_203	14	2,000	2,730	17	21,622	56,757
label: cooccurence network weight distribution and sparsity

Osservando i dati in tabella si osserva che la median è pari a 2 in tutte le comunità ad eccezzione della comunità 81 la quale è un caso limite avendo 3 nodi. Inoltre, generalmente, in tutte le community,  si osserva mean > median, indicativo di distribuzioni asimmetriche (con coda lunga a destra), caratterizzate dalla presenza di pochi archi con peso elevato e di molti archi con peso ridotto. In alcune community (ad esempio Community 5) questa differenza tra media e mediana è più marcata, mentre in altre risulta più contenuta (ad esempio Community 20). Ciò suggerisce la presenza di coppie di geni fortemente co-targettate, immerse in un contesto di co-targeting complessivamente più moderato.

Il parametro max misura peso massimo osservato nella community ovvero il numero massimo di farmci condivisi da una coppia di geni. In altre parole, è un indice di quanto può essere forte il co-targetting in quella community. Ci sono molte community con un max elevato indice che in tutte le reti esiste almeno una coppia di geni co-targettata da un numero alto di farmaci.

I parametri weight_eq_1_pct e weight_eq_2_pct sono stati introdotti come informazioni aggiuntive e rappresentano rispettivamente la quota di archi gene–gene sostenuti da uno e da due farmaci. Si osserva che weight_eq_1_pct è sempre ben al di sotto del 50% (ad eccezione di casi patologici di dimensione molto ridotta), mentre weight_eq_2_pct risulta spesso la classe dominante (≈ 40–70%). Questo comportamento è coerente con una mediana dei pesi pari a due, che implica che almeno la metà degli archi presenti nelle co-occurrence networks sia sostenuta da due o più farmaci.
