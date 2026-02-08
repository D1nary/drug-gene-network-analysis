# Descrizione del dataset

# Objective
L'obiettivo di questo progetto è quello di analizzare strutturalmente ...

# Similarity network
È stata costruita una similarity network drug–drug utilizzando la Jaccard similarity, applicando inizialmente un threshold pari a 0.3 e successivamente uno pari a 0.4. I risultati ottenuti con le due soglie sono stati poi confrontati tra loro (file filtering.json del programma).

La Jaccard similarity misura il grado di sovrapposizione relativa tra due farmaci ed è definita come:
$$
J(A,B) = \frac{|A \cap B|}{|A \cup B|}
$$

Ogni nodo del grafo della similarity network rappresenta un farmaco, caratterizzato da un profilo bersaglio, ovvero da un vettore contenente tutti i geni (target) con cui il farmaco interagisce. Un esempio di rappresentazione di un nodo è il seguente:
```bash
{
  "bipartite": "drug",
  "original_id": 12345,
  "targets": [1017, 1956, 7422]
}
```

La Jaccard similarity indica che due farmaci risultano simili solo se condividono un numero significativo di target e, allo stesso tempo, presentano pochi target differenti, cioè se hanno un profilo di bersagli complessivamente simile.

È stata scelta la Jaccard similarity perché opera nativamente su insiemi, penalizzando farmaci con un elevato numero di target non condivisi e riducendo l’impatto di overlap piccoli ma casuali. Questo la rende particolarmente adatta all’analisi della similarità tra profili di target nel contesto del dataset analizzato.




## Threshold 0.3

| Parametro            | Valore    |
| -------------------- | --------- |
| similarity_threshold | 0.3       |
| nodes_removed        | 333       |
| edges_filtered       | 961,999   |
| original_node_count  | 1,774     |
| retained_node_count  | 1,441     |
| potential_edges      | 1,037,520 |

Questo valore di threshold indica che un arco (drug–drug) è presente nella rete solo se la similarità tra due farmaci è ≥ 0.3. Essa è una soglia moderata non troppo permissiva ma nemmeno troppo stringente. Mantiene connessioni con similarità medio–bassa, quindi preserva una rete relativamente densa rispetto a soglie più stringenti. Genera una rete connessa ma non troppo, utile per il community detection riducendp il rumore ed eliminando similarità assolutamente deboli. 

Siccome il numero di nodi originali (original_node_count) è 1774 e il numero di nodi "sopreavvisuti" (retained_node_count) è 1441, si osserva che, dopo il filtraggio, la rete conserva l'81,2% dei farmaci. Oltre 1400 farmaci sono connessi tra di loro significa, con questo threshold, il dataset contiene informazioni ridondanti sui target farmacologici

Il parametro potential_edges rappresenta il totale delle coppie drug-drug possibili prima del filtraggio mentre edges_filtered indica quante coppie sono state eliminate perché avevano similarità < 0.3. Sottraendo questi due valori, si ottiene che la rete finale ha 85310 archi. Quindi solo circa l’8.22% delle possibili connessioni supera la soglia, indicando una rete relativamente sparsa. Questa sparsità è coerente con reti di similarità farmacologica poichè, in genere, pochi farmaci sono veramente simili.



## Threshold 0.4
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

Il numero di nodi rimossi rimane invariato a 333. Ciò indica che questi farmaci non presentano valori di similarità pari o superiori a 0.4 con nessun altro farmaco del dataset, risultando quindi isolati nella similarity network già al threshold più permissivo.

Uno degli scopi di questa analisi, è quello di costruire una community network per ricercare comunità di farmaci con meccanismài d'azione simili. In questo contesto, non avrebbe senso alzare troppo del threshold potrebbe avere l'effetto di tenere solo le relazioni “ovvie” e perdere quelle deboli ma biologicamente interessanti. Infatti anche similarità più deboli rappresentano un caso interessante poichè, due farmaci che condividono pochi target ma “giusti” possono essere candidati per drug repurposing. Un farmaco con profilo di target solo parzialmente simile può comunque avere effetti collaterali simili oppure opportunità di combinazione terapeutica. Imporre soglie troppo rigide rischia di semplificare eccessivamente una realtà biologica che è invece caratterizzata da relazioni sfumate e parziali.

Le seguenti analisi sono state eseguite considerando un threshold di 0.4.

---
NOTA:
**Come mai vengono rimossi nodi durante la costruzione della similarity network?**
Vengono rimossi nodi che non superano il threshold con nessuno

**Ma allora come mai, nella rappresentazione della similarity, vengono rappresentati nodi non collegati?**
Nella visualizzazione della similarity network sono stati mantenuti anche i nodi isolati, corrispondenti a farmaci che non superano la soglia di similarità con nessun altro composto. Sebbene tali nodi vengano esclusi dalle analisi di connettività e community detection (conteggiati in nodes_removed), la loro presenza nello snapshot consente di evidenziare la forte eterogeneità dei profili target e la natura selettiva del threshold adottato.
---




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
Una modularity di circa 0.20 indica la presenza di una struttura modulare debole ma reale, suggerendo che i farmaci non sono distribuiti in modo puramente casuale, ma mostrano una tendenza a organizzarsi in gruppi caratterizzati da profili di similarità parzialmente condivisi.

Pur non trattandosi di una modularità elevata, questo valore è coerente con la natura complessa dei sistemi biologici e farmacologici, nei quali molti farmaci presentano comportamenti polifarmacologici e condividono target o pathway tra differenti classi funzionali. Di conseguenza, la similarity network drug–drug non evidenzia una separazione netta in moduli fortemente isolati, ma una struttura interconnessa e continua, nella quale i collegamenti tra comunità riflettono sovrapposizioni funzionali, potenziali effetti collaterali comuni e possibili opportunità di combinazione terapeutica o drug repurposing.

## Other parameters
- Median:
Il valore mediano delle comunità risulta essere 2. Questo significa che più della metà delle comunità ha 2 farmaci o meno. Questo suggerisce gruppi di farmaci molto simili fra di loro i quali possono essere varianti strutturali di uno stesso composto oppure farmaci che condividono uno o pochissimi target molto specifici. In altre parole, queste comunità rappresentano farmaci di nicchia con bersagli rari oppure possono essere outliers interessanti per riposizionamento. Per esempio: se un singleton si collega debolmente a una grande comunità, potrebbe condividere qualche pathway con farmaci di un’altra indicazione.

- Mean, max community 
La media è rappresentata da una comunità di circa 5 farmaci. Questo valore insieme a quello della mediana, indicano una distribuzione fortemente sbilanciata con pochi cluster molto grandi e tanti piccoli. 

METTI HISTOGRAMMA DELLE COMUNITÀ IN FUNZIONE DEL NUMERO DI ELEMENTI.

Il valore massimo è rappresentato da una comunità con 359 elementi. In particolare, grandi famiglie farmacologiche come questa rappresentano famiglie che colpiscono la stessa famiglia di proteine (es. molte chinasi, GPCR, recettori nucleari…) oppure gli stessi pathway coinvolti in malattie comuni (es. segnali proliferativi tumorali, infiammazione, ecc.).  Questa comunità verrà approfondita meglio nelle sezioni sucessive.



## Inforamazioni biologiche ricavabili da questi valori

Anche in assenza di un’analisi dettagliata della composizione interna di ciascuna comunità, i parametri ottenuti tramite il metodo di Louvain consentono di trarre alcune considerazioni biologiche rilevanti. In primo luogo, la presenza di una modularità diversa da zero indica che la rete di farmaci non è organizzata in modo casuale, ma presenta una struttura modulare, suggerendo l’esistenza di gruppi di farmaci che condividono meccanismi molecolari, bersagli biologici o ambiti terapeutici affini.

La modularità moderata osservata riflette inoltre un’elevata polifarmacologia del sistema in molti farmaci risultano connessi a più comunità attraverso target o pathway condivisi. Questo comportamento è coerente con la possibilità che farmaci appartenenti a classi terapeutiche differenti possano indurre effetti collaterali simili o presentare opportunità di riposizionamento farmacologico, qualora condividano un sottoinsieme rilevante di bersagli biologici.

La distribuzione delle dimensioni delle comunità, caratterizzata dalla presenza di numerosi moduli di piccola dimensione e di nodi isolati, suggerisce l’esistenza di farmaci con profili di target rari o altamente specifici. Tali regioni del network possono rappresentare meccanismi d’azione poco ridondanti o ancora scarsamente esplorati potendo essere interassanti spunti di partenza per l'analisi o la scoperta di nuovi farmaci

Al contrario, la presenza di pochi moduli di grandi dimensioni può evidenziare aree del network associate a target e pathway ampiamente studiati e frequentemente sfruttati nello sviluppo farmacologico. In queste comunità si osserva un’elevata sovrapposizione farmacologica, che può essere utile per confrontare profili di efficacia e tossicità tra farmaci simili, nonché per individuare strategie di combinazione terapeutica all’interno della stessa comunità o tra comunità strettamente correlate.

## Communities analisys
Le singole comunità vengono analizzate nel file community_parameters.csv. Il file è organizzato in colonne ciascuna con una caratteristica della comunità:
- community id: Identificativo della comunità
- size: Dimensione della community
- degree: Numero di arhci del nodo
- weighted degree: Somma dei pesi degli archi incidenti su un nodo
- clustering coefficient: Grado di chiusura locale delle comunità

INSERISCI TABELLA CON TUTTI I DATI

### Size
Biologicamente, la dimensione di una comunità riflette il numero di farmaci che condividono un profilo di target genici simile, rappresentando quindi la densità di uno specifico spazio farmacologico nel grafo di similarità. Poiché il metodo di Louvain individua gruppi con elevata connettività interna, comunità di grandi dimensioni indicano insiemi di farmaci caratterizzati da profili di target altamente sovrapposti, compatibili con pathway biologici centrali o con una marcata ridondanza terapeutica. Al contrario, comunità di piccole dimensioni emergono in presenza di profili di target più specifici o poco condivisi, suggerendo spazi farmacologici più specializzati o meno ridondanti.

Mentre comunità grandi possono rappresentare meccanismi già esplorati (poichè diversi farmaci agiscono sullo stesso bersaglio o sullo stesso pattern), le comunità piccole individuano target non molto esplorati farmacologicamente oppure target presenti in pathway meno studiati e quindi famiglie proteiche non ancora saturate da altre molecole. Inoltre, se un farmaco non clusterizza con gli altri significa che non condivide un numero sufficiente di target da superare il threshold quindi non si inserisce nei pathway delle classi terapeutiche classiche o che potrebbe modulare il sistema biologico in modo differente. Cioè, essendo il sistema biologico una rete complessa di interazioni (proteine, pathway, feedback, ecc), un farmaco appartenente ad una piccola comunity, può alterare la rete in un punto con una dinamica diversa dalla maggiorparte dei farmaci esistenti. Future analisi di piccole comunità di farmaci possono essere utili in ambito farmacologico

Per semplicità di visualizzazione e per coerenza con l'analisy della density, nella tabella precedente sono riportate solo le comunità con una size maggiore o uguale di 20

MANCA L'ANALISI STATISTICA DELLE SIZE TROVATE ALL'INTERNO DELLA COMUNITY NETWORK. POSSIBILE ISTOGRAMMA DELLE SIZE DELLE COMMUNITY



### Density

La densità misura quanto i farmaci all’interno di una comunità siano effettivamente simili tra loro. Anche se il metodo di Louvain individua gruppi ben separati nel grafo, alcune comunità possono essere tenute insieme da similarità parziali o indirette. In questi casi la densità aiuta a distinguere moduli realmente omogenei da insiemi di farmaci che condividono solo alcuni bersagli o funzioni biologiche, senza costituire una vera classe farmacologica compatta. La densità è stata calcolata con:

$$
\text{density} = \frac{E_{\text{int}}}{\frac{n(n-1)}{2}}
$$
dove:

- density: frazione di connessioni interne presenti rispetto al massimo possibile nella comunità.
- E_int: numero di archi interni alla comunità (solo tra nodi della stessa comunità).
- n: numero di nodi nella comunità (size).
- n(n-1)/2: numero massimo di archi possibili in una comunità non orientata senza self‑loop.

Di seguito riportiamo in tabella solo i valori di densità per comunità con un numero maggiore o uguale a 15

ID              Size   density

Community_5     359   1
Community_13    18    0.791
Community_20    22    0.792
Community_21    82    0.778
Community_22    83    0.913
Community_54    22    0.745
Community_78    32    0.905
Community_81    20    1
Community_97    19    0.754
Community_135   16    0.425
Community_188   72    0.324
Community_192   39    0.601
Community_203   23    0.482

Nei moduli piccoli, il numero di connessioni possibili è molto limitato, infatti basta che pochi nodi siano tutti connessi tra loro perché la density risulti elevata, spesso prossima a 1. Questo riflette una similarità molto forte tra i farmaci del gruppo (ad esempio condivisione quasi completa dei target), ma tali valori sono poco robusti dal punto di vista statistico, perché fortemente influenzati dal basso numero di nodi. Per questo sono stati rimossi dall'analisi considerando solo comunità con size maggiore o uguale di 15. 

Dai dati riportati si osserva che comunità medio-grandi, presentano valori di densità variabile. Analizzando le comunità  con size compresa tra 39 e 72 si osservano densità molto diverse. Ovvero comprese tra 0.324 e 0.913. 
Le cause di questo comportamento possono essere sia matematiche che biologiche. 

Le cause matematiche sono da attribuirsi a come viene calcolata la density (formula precedente). Dal momento che il numero di possibili connessioni interne aumenta quadraticamente con la dimensione della comunità, anche piccole variazioni nei profili dei nodi tendono ad accumularsi in un numero crescente di archi mancanti. Questo effetto è amplificato nelle comunità più grandi, dove l’elevato numero di confronti tra profili rende più probabile l’emergere di differenze, riflettendosi in una maggiore variabilità dei valori di densità.

Inoltre, la presenza o l'assenza di un arco dipende da una soglia (0.4). Quindi, se due nodi sono appena sotto la soglia l'arco nella similarity network è assente mentre, se la similarity è appena sopra la soglia, l'arco è presente. Questo introduce variabilità artificiale e maggiore dispersione dei peorfili di densità. 

Da un punto di vista biologico, comunità grandi spesso aggregano famiglie farmacologiche ampie, pathway complessi o target parzialmente sovrapposti creando, di conseguenza cluster di più grandi dimensioni ma non completamente connessi.

Per le comunità di dimensione più contenuta (size compresa tra 16 e 32), si osserva in generale una densità mediamente elevata, con valori tipicamente compresi tra ~0.6 e ~0.9, e solo rare eccezioni con densità più bassa (es. Community_135 con density = 0.425).

Questo può essere spiegato matematicamente dal fatto che comunità piccole tendono, per costruzione, a mostrare valori di density più stabili e meno dispersi rispetto a quelle grandi. In queste comunità infatti, sebbene ogni singolo arco contribuisca in modo relativamente maggiore al valore della density nelle comunità di piccole dimensioni, una riduzione significativa della densità richiede l’assenza sistematica di molte connessioni interne. In assenza di eterogeneità strutturata nei profili dei nodi, le comunità piccole tendono quindi a mantenere valori di densità stabili.

Le rare comunità piccole con densità più bassa possono indicare quindi eterogeneità strutturale reale.

ANALIZZA IL DATASET, NON CERCARE UNA SPIEGAZIONE MATEMATICA O BIOLOGICA PER TUTTO

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

Un'ulteriore causa del valore quasi unitario di questa community, potrebbe provenire dalla natura intrinseca del dataset. Infatti esso, aggregando screening diversi producendo il risultato che, farmaci testati negli stessi screening sugli stessi pannelli genici con risultati simili diventano vettorialmente indistinguibili. Questo produce sottografi altamente densi nella rete di similarità (profili identici -> tutte le connessioni di similarità sono presenti -> con threshold basso tutti gli archi sono presenti -> density = 1)

--- 
Farmaci testati negli stessi screening:

- screening 1: testa A e B su geni {G1, G2, G3}
- screening 2: testa A e B sugli stessi geni {G1, G2, G3}
- screening 3: testa A e B sugli stessi geni {G1, G2, G3}

screening diversi = esperimenti diversi (utilizzando stessi farmaci e geni target)
---

All'interno del dataset, molti farmaci sono distinti per nome o contesto, ma non per profilo genetico producendo nodi distinti con vettori profilo identici e archi completi. Ogni farmaco è un nodo distinto identificato da un ID chimico con possibilità che provenga da fonti diverse, screening diversi o contesti sperimentali divesti. In altre parole due molecole possono essere chimicamente diverse o chimicamente molto simili o lo stesso composto annotato in contesti divesi. Si ha la produzione di un clique artificilae dal punto di vista topologico. 

Biologicamente, senza nessuna analisy più prodonda, questa comunità rappresenta una classe di composti con meccanismo d'azione quasi identico oppure il targeting di uno stesso grande pathway o complesso genico. Per capire meglio la natura di questa communiti andrebbero condotte analisi più approfondite. 

È inoltre presente una comunità con size 20 e density 1. In questo caso, il numero di archi possibili sono 190 e, similmente al caso sprecedente, il numero di archi presenti nella comunity è 190. In partiolare 19 farmaci su 20 hanno lo stesso identico profilo target il quale contiene solo 2 geni. Questo caso rappresenta un clique più banale del precedente poichè il profilo è piccolissimo ed essendo in preseza di una replica quasi perfetta. 

### Clustering coefficient
Per ciascuna comunity, è stato calcolato il clustering coefficient:

community_id size density clustering_coefficient

Community_5 359 1.000 0.903
Community_20 22 0.792 0.593
Community_21 82 0.778 0.631
Community_22 83 0.913 0.851
Community_54 22 0.745 0.574
Community_78 32 0.905 0.653
Community_81 20 1.000 0.964
Community_188 72 0.324 0.616
Community_192 39 0.601 0.706
Community_203 23 0.482 0.513

Globalmente si può osservare che le comunità hanno valori mediamente alti di clustering coefficient. Questo suggerisce che le comunità individuate non presentano una struttura lineare o “a catena”, in cui i farmaci risultano simili solo a pochi vicini immediati, ma piuttosto costituiscono moduli fortemente coesi. Cioè, se un farmaco è simile ad altri farmaci all'interno della comuintà, è molto porbabile che esso sia siamile ad altri componenti della stessa.

Le community 20, 21, 22, 54, 78 presentano valori compresi tra:
| Size  | Density   | Clustering |
| ----- | --------- | ---------- |
| 22–83 | 0.74–0.91 | 0.57–0.85  |

quindi con densità alta ma minore di 1 e clustering coefficient alto. Da questi dati si può dire che non tutti i farmaci sono simili agli altri ma esistono sottogruppi molto coerenti all'interno della comunità. In altre parole, la comunità è composta da blocchi locali (blocchi di traingoli) fortemente connessi. Biologicamente, è possibile la presenza di famiglie farmacologiche in cui esiste un nucleo di target comuni ed in cui ogni farmaco può introdurre cariazioni marginali sul profilo. Possono essere presenti sotto-meccanismi d'azione in cui i geni della comunity condividono un nucleo di geni bersaglio ma ciascuno interagice con geni aggiuntivi diversi. 

Sono presenti inoltre, comunità con bassa density e clustering moderato
---
Density bassa
→ molti farmaci non sono direttamente simili tra loro
→ la comunità non è un blocco compatto

Clustering moderato
→ quando un farmaco è simile a due altri, quei due tendono comunque a essere simili tra loro
→ esistono triangoli locali, cioè sottogruppi coerenti

👉 Questo implica che la comunità è composta da più sottogruppi locali (cluster densi) connessi da pochi nodi (farmaci ponte).
---
(questo è possibile poichè la density guarda la similarità globale mentre il clustering guarda solo ai vicini comuni)

Tra le comunità con una density più bassa abbiamo:
| Community | Size | Density | Clustering |
| --------- | ---- | ------- | ---------- |
| 188       | 72   | 0.324   | 0.616      |
| 203       | 23   | 0.482   | 0.513      |

Queste mostrano poche connessioni globali ma connessioni locali ben strutturate. In altre parole, la comunità non è un blocco compatto ma un insieme di cluster locali collegati indirettamente. Biologicamente, siamo in presenza di profili target eterogenei e farmci che condividono solo alcune componenti funzionali o pathway comuni. Quindi possibile presenza di pathway parzialmente condivisi farmaci ponte tra meccanismi diversi. In questo caso il clustering coefficient è molto informativo poichè rivela coerenza locale della community nascosta dalla bassa density.


Un'ultimo caso che si può notare dalla tabella precedente sono comunità con density $\approx 1$ e clustering alto. 
| Community | Size | Density | Clustering |
| --------- | ---- | ------- | ---------- |
| 5         | 359  | 1.000   | 0.903      |
| 81        | 20   | 1.000   | 0.964      |
In queste comunità il clustering coefficient non aggiunge nuova informazione rispetto alla density, ma rafforza l’evidenza di omogeneità estrema.

## Weighted degree
Per descrivere il livello complessivo di interazione di una comunità con le altre è stato calcolato, il weighed degree. Esso è calcolato come la somma dei pesi degli archi inter-comunità (similarità Jaccard). 

community_id	size	weighted_degree
Community_5	359	0.0
Community_13	18	0.0
Community_20	22	0.0
Community_21	82	0.0
Community_22	83	12.990404040404039
Community_54	22	0.0
Community_78	32	12.990404040404039
Community_81	20	0.0
Community_97	19	0.0
Community_135	16	0.0
Community_188	72	0.0
Community_192	39	0.0
Community_203	23	0.0


Solo le comunità 22 e 78 risultano essere connesse tra di loro con un $\text{weighted degree} = 12.99$.

Siccome, un arco tra due comunità esiste solo se esistono farmaci appartenenti a comunità diverse ma connessi e che quindi condividono una porzione significativa di target genici (Jaccard maggiore uguale a 0.4), le due comunità rappresentano moduli farmacologici distinti ma non indipendenti. Le due comunità possono essere associate a pathway diversi ma interconnessi con una sovrapposizione di geni chiave. I farmaci ponte che connettono le due comunità colpiscono geni presenti in entrambe le comunità.

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
 
## component count e giant component size
Tutte le co-occurrence networks gene-gene hanno component_count = 1
Questo vale indipendentemente da:
- numero di nodi (da 2 a 357),
- densità (da ~0.75 a 1),
- dimensione della community di farmaci di origine,
- clustering coefficient.

Questo indica che:
- la rete gene-gene è completamente connessa;
- ogni gene è raggiungibile da ogni altro gene tramite almeno un cammino;
- non esistono sottogruppi genetici isolati all’interno della community.

Coerente con la scelta metodologica di costruire le reti in comunità poichè aumenta la probabilità di profili genetici sovrapposti.

Dal punto di vista biologico, component_count = 1 suggerisce che:

- i geni all’interno di ciascuna community di farmaci partecipano a un sistema funzionale interconnesso;
- i farmaci della community colpiscono moduli genetici fortemente sovrapposti, non pathway indipendenti;
- non emergono sottosistemi genetici separati, ma un unico “blocco” funzionale.


Inoltre, in tutte le le co-occurrence networks gene-gene analizzate si osserva che: giant_component_size = n_nodes ovvero la componente gigante coincide sempre con l’intera rete. In questo scenario, il concetto di giant component perde la sua usuale accezione di “sottostruttura dominante” e diventa equivalente alla rete stessa. Questo parametro, non aggiunge nuova informazione ma verifica di coerenza strutturale della rete e conferma dell’elevata connettività genetica interna alle community.

## Density e global clustering coefficient
La connettività globale della rete è già garantita (component conunt = 1) per costruzione. In questo caso la density ci dice quanto è ridondante la condivisione dei farmaci tra i geni della comunity

Calcolata con la stessa formula del caso precedente
Valorielevati
assenza di reti sparse (<0.5)

Possiamo identificare tre regimi di density. Tutti e tre comunque con valori elevati
- regime 1:
density∈[0.95,1.00]
Tuttue o quasi tutte le coppie di geni sono collegate
reti clique (o quasi)
assenza quasi totale di eterogeneità dei profili
geni colpiti sempre o quasi dagli stessi faramci
forte ridondanza farmacologica
moduli genetici sovrapposti indicativi di stesso pathway, stesso complesso proteico, repliche o contesti sperimentali simili

- regime 2:
density∈[0.75,0.90)
Reti ancora pienamente connesse
numero più o meno significativo di archi mancanti
struttura meno clique leggermente più articolata della precedente
farmaci con target più specifici o combinazioni diverse di bersagli. 
geni coinvolti in processi distinti ma interconnessi



---
In questo caso tutte le co-occurrence networks hanno:
- component_count = 1
- giant_component_size = n_nodes

👉 quindi:
- la connettività globale è già garantita per costruzione,
- la density non sta più dicendo se la rete è connessa o no.

(anche una rete con density = 0.75 è topologicamente connessa quanto una con density = 1.0)

La density in generale ci dice quanto la rete è “piena” di archi rispetto al massimo possibile. Se una rete è piena di archi vuol dire che i geni al suo interno sono collegati. Essendo che, nella cooccurence network coppie di geni sono collegati se condividono almeno un farmaco, avere una density alta significa che geni condividono stessi farmaci e quindi che la condivisione dei farmaci è ridondante tra i geni della community
---

---
Un pathway dominante è un insieme di geni che partecipano allo stesso meccanismo biologico e che vengono frequentemente colpiti insieme dai farmaci di una community. In questo caso i farmaci non agiscono su singoli geni isolati, ma interferiscono in modo coordinato con un processo biologico specifico, che rappresenta il bersaglio funzionale principale della community.
---

---
Il global clustering coefficient misura la probabilità che due geni entrambi co-targettati con un terzo gene siano a loro volta co-targettati.
In altre parole quantifica quanto la rete è localmente “triangolare”, quindi quanto è vicina a una struttura clique-like.
---

Anche nel caso del clustering coefficient, siamo in presenza di valori molot alti
- regime di clustering massimo
ogni tripla di geni forma un triangolo
struttura completamente ridondanteGeneralmente, in tutte le community si osserva mean > median ovvero distribuzioni asimmetriche (con una coda lunga a destra) con pochi archi aventi un grande peso e molti con peso ridotto. In alcune (comunity 5) questa differenza tra mean e median è più marcata mentre in altre meno (community 20)
Questo vuol dire che esistono coppie di geni estremamente co-targettate immerse in un co-targetting relativamente più moderato
- regime di clustering alto ma < 1 (≈ 0.88 – 0.95)
struttura localmente compatta ma non completamente chiusa.
moduli genetici coerenti ma articolati;
pathway affini o parzialmente (anche se fortemente) sovrapposti;
presenza di geni “ponte” che collegano sottosistemi funzionali vicin

Siamo in presenza di un unico caso in cui il clustering = 0. Questa è una rete composta da 2 geni quindi assenza di trangoli per definizione

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

Medianan pari a 2 in tutte le comunità ad eccezzione della comunità 81 la quale è un caso limite avendo 3 nodi
Due o più farmaci che colpiscono la stessa coppia di geni può suggerire un legame funzionale reale e non un evento accidentale.

Generalmente, in tutte le community si osserva mean > median, indicativo di distribuzioni asimmetriche (con coda lunga a destra), caratterizzate dalla presenza di pochi archi con peso elevato e di molti archi con peso ridotto. In alcune community (ad esempio Community 5) questa differenza tra media e mediana è più marcata, mentre in altre risulta più contenuta (ad esempio Community 20). Ciò suggerisce la presenza di coppie di geni fortemente co-targettate, immerse in un contesto di co-targeting complessivamente più moderato.

Biologicamente si ha la possibile presenta di geni "hub", pathway centrali o target riccorrenti in più screening

Il parametro max misura peso massimo osservato nella community ovvero il numero massimo di farmci condivise da una coppia di geni nella community ovvero un indice di quanto può essere forte il co-targetting in quella community.
Molte community con valori elevati di max. 
In tutte le reti esiste almeno una coppia di geni co-targettata da un numero alto di farmaci.
Reti costituite da nuclei di co-targeting ridondandi

I parametri weight_eq_1_pct e weight_eq_2_pct sono stati introdotti come informazioni aggiuntive e rappresentano rispettivamente la quota di archi gene–gene sostenuti da uno e da due farmaci. Si osserva che weight_eq_1_pct è sempre ben al di sotto del 50% (ad eccezione di casi patologici di dimensione molto ridotta), mentre weight_eq_2_pct risulta spesso la classe dominante (≈ 40–70%). Questo comportamento è coerente con una mediana dei pesi pari a due, che implica che almeno la metà degli archi presenti nelle co-occurrence networks sia sostenuta da due o più farmaci.


---
I parametri di centralità di nodo, quali closeness, betweenness e PageRank, non sono stati inclusi nell’analisi in quanto l’obiettivo è la caratterizzazione globale delle co-occurrence networks. Inoltre, l’elevata densità e la connettività completa delle reti renderebbero tali misure poco discriminanti in questo contesto.
---

# Spettral analisys
The dug-drug similarity network was built with a threshold of...
In modo da codificare la struttura di connettività interna di ciascuna community di farmaci sono state costruite le normalizzed Laplacian matrix e di esse, lo spettro degli eigenvalues è stato analizzato. 

Sebbene, in una precedente analisi, siano stati calcolati parametri come density e clustering coefficient fornendo una prima caratterizzazione topologica delle community, l’analisi spettrale della Laplaciana normalizzata consente di valutare la coesione strutturale globale e l’eventuale presenza di sottostrutture modulari non rilevabili con le metriche appena citate. 

La density, per esempio, è una caratterizzazaione globale indicante quanti archi sono presenti nella community. Infatti, due community con lo stesso valore di density, possono essere strutturalmente diverse (comunità quai-clique oppure costituita da due blocchi densi separati). In questo contesto, è difficile distinguere comunità realmetne omogenee da aggregazioni indotte dalla misura di similarità.

Per analizzare la connectivity e la cohesivity delle reti, sono stati calcolati i Fiedler values ovveri il secondo eigenvalue più piccolo. Nella seguente tabella, sono riportati i Fiedler values per le comunità con size maggiore o uguale di 15.

Community name	Size	Fiedler value
Community_5	359	0.971
Community_13	18	0.270
Community_20	22	0.384
Community_21	82	0.514
Community_22	83	0.611
Community_54	22	0.672
Community_78	32	0.830
Community_81	20	1.036
Community_97	19	0.263
Community_135	16	0.029
Community_188	72	0.006
Community_192	39	0.045
Community_203	23	0.112

Per effettuare un'analisi dei Fiedler value individuati, introduciamo le seguenti classi operative
| Range λ₁      | Interpretazione strutturale                       |
| ------------- | ------------------------------------------------- |
| **λ₁ ≪ 0.1**  | Community **quasi separabile**, struttura fragile |
| **0.1 – 0.4** | Modularità interna marcata                        |
| **0.4 – 0.8** | Strutturalmente coesa con eterogeneità interna    |
| **0.8 – 1.1** | Community **fortemente coesa**                    |
| **> 1.1**     | Quasi-clique / struttura estremamente compatta    |

Nel dataset, esistono, sia comunità spettralmente deoli ((es. Community 100, 135, 188)) sia comunità spettralmente molto forti con λ₁ ≈ 1 o maggiore. (Questo non sarebbe distinguibile usando solo density e clustering. inserire?)


GRAFICO FIELDLER VALUE IN FUNZIONE DELLA SIZE?

Dal grafico, si osserva che, per comunità più piccole, il Fiedler value può assumere valori sia piccoli che grandi identificando comunità sia comunità con una grande connettività che altre con connettività ineriore. Per quanto riguada comunità più grandi (size > 50), a parte un solo caso, esse risultano essere più connesse.

Ricapitolando, la density ci da informazioni solo su quanti archi sono presenti e non come sono essi sono distribuiti o se sono presenti separazioni interne mentre il Fidelr value fornisce infomrazioni su quanto facile o costoso separare il grafo in due parte minimizzando il peso degli archi tagliati (tutorial “Algorithms for Graph Partitioning”). 
Sulla base di questo, se la density è alta (molti archi nella rete) e valuenè alto (difficoltà a separare il grafo) vuol dire che siamo in presenza di una ridondanza reale poichè ci sono archi uniformemente distribuiti, nessun sottogruppo separabile, ogni nodo è connesso “bene” con tutti.
Al contrario, se value è basso è "semplice" separare il grafo, siamo in presenza di una aggregazione artificiale ovvero una community che appare densa e compatta per costruzione matematica (similarità + threshold + algoritmo), ma che non rappresenta un insieme biologicamente omogeneo. Siamo in presenza di Due (o più) sottogruppi internamente molto densi e pochi archi tra i sottogruppi. 


Nel ChG-InterDecagon questo succede quando Jaccard è permissiva su set grandi, farmaci condividono parte dei target ma appartengono a pathway diversi. 
In questo caso il Louvain li mette insieme quindi la density resta alta ma spettralmente la community non regge



Community name	Size	density	Fiedler value
Community_5	357	0.755	0.971
Community_13	166	0.995	0.270
Community_20	17	1.000	0.384
Community_21	26	0.862	0.514
Community_22	77	0.876	0.611
Community_54	4	1.000	0.672
Community_78	89	0.867	0.830
Community_81	2	1.000	1.036
Community_97	215	0.963	0.263
Community_135	86	0.783	0.029
Community_188	82	0.781	0.006
Community_192	34	0.761	0.045
Community_203	5	1.000	0.112

Riprendendo l'analisi delle sezioni precedenti, dato che tutte le community elencate hanno density alta (≈ ≥ 0.75), la discriminante reale è il Fiedler value. Possiamo dividere le comunità (size >= 15) in 3 gruppi diversi:

Fiedler basso: λ₁ ≲ 0.1
coesione globale debole, quasi separabili
Di questo gruppo ne fanno parte le community 135, 188 e 192


Fiedler intermedio: 0.1 < λ₁ < 0.7
community dense ma strutturalmente eterogenee (community con una certa modularità interna)
Di questo gruppo ne fanno parte le community 13, 20, 21, 22, 97 e 203
Molte di queste sono comunità quasi-cliqie locali

Fiedler alto: λ₁ ≳ 0.7
community globalmente coese, ridondanza reale
Di questo gruppo ne fanno parte le community 5, 78, 54 e 2 (banale matematicamente)

In questo caso la community 5 precedentemente analizzata  non è un artefatto della Jaccard o del threshold, ma riflette una ridondanza reale dei profili target

SPETTRO SOLO COMMUNITY 5? CONFRONTO CON COMUNITÀ NON CLIQUE?

# DAG
All'interno delle community, abbiamo farmaci simili per costruzione (Jaccard similarity). Un'informazione che le community non forniscono, sono delle relazioni di "generalità"/"specifità" all'interno di community con farmaci simili. 
Se due farmaci sono somili all'interno della community, allora è prbalile che condividono un core di geni e differiscono per una piccola periferia (pochi geni in più/in meno). Le DAG costruite e analizzate di seguito (in particolare le loro orientation rule) vogliono analizzare proprio questo fatto. 

## Community analizzata
La seguente community è stata scelta per l'analisi perchè è una comunità informativa ovvero non banale (size piccola) non triviale (clique o quasi-clique) ma strutturamente eterogenea e biologicamente interpretabile nodo per nodo
Quello che la segunete analisi vuole fare è capire che ruolo ha ogni farmaco dentro la community.

Cerca comunità con:
- size 30 ≤ size ≤ 150
- density alta ma non satura: 0.7 ≤ density ≤ 0.9
   - Questo significa:
      - esiste coerenza interna
      - ma non tutti sono simili a tutti
      - quindi c’è eterogeneità strutturale reale
- Clustering coefficient alto ma non estremo: 0.6 ≤ clustering ≤ 0.9
   - Interpretazione:
      - sottogruppi locali ben definiti
      - possibilità di famiglie di farmaci o sottopercorsi biologici
- Fiedler value non estremo (chiave!): Fiedler intermedio (né ≪ 0.1, né ≫ 0.8)
   - Questo identifica comunità:
      - ben connesse
      - ma non rigidamente compatte
     - con colli di bottiglia strutturali biologicamente interessanti
     
La migliore candidata per questa analisi è la community 21 la quale possiede i seguenti parametri:
- size = 82
- density ≈ 0.78
- clustering ≈ 0.63
- Fiedler ≈ 0.51

## Orientation rule
La DAG è stata costruita con la seguente orientation rule:

- Si ordina per dimensione dei target. 
- Si crea un arco dal set più grande al set più piccolo se il set piccolo è sottoinsieme del grande: A -> B se T(B) ⊂ T(A). 
- C’è anche un vincolo sulla differenza di cardinalità: si considerano solo coppie con 0 < |T(A)| - |T(B)| <= min_set_difference (di default 3). 
- Se la differenza è maggiore, il ciclo si interrompe per quel node_small.

In questo modo:
La regola  𝐴 → 𝐵 se  T(B)⊂T(A):
- A = profilo più “generale” (include tutto ciò che fa B + extra)
- B = profilo più “specifico” (un sottoinsieme del generale)

Inoltre, il nodo viene creato, solo se Quindi A può avere al massimo 3 target in più rispetto a B. Questo è fondamentale per una rappresentazione gerarchica locale e per non perdere di interpretabilità biologica

---
Se A ha molti più target di B, le possibili cause sono indistinguibili:
- A è stato testato in molte più condizioni
- B è incompleto (missing data)
- A ha effetti aspecifici / tossici
- aggregazione di screening diversi
- annotazioni ridondanti nel dataset

Con una differenza piccola:
- queste cause sono meno plausibili
- o comunque più distinguibili
---

In questo modo analizziamo la rete individuando le seguenti tipologie di nodi:
- Sorgenti: Farmaci massimali
- Sink: Farmaci minimali
- Nodi intermedi: varianti incrementali


## DAG global parameters
+-------------+--------+
| parametro   | valore |
+-------------+--------+
| n_nodes     | 82     |
| n_edges     | 681    |
| density_dag | 0.103  |
| n_sources   | 35     |
| n_sinks     | 21     |
| max_depth   | 6      |
+-------------+--------+

Il valore di density osservato, suggerisce che la gerarchia di inclusione tra i nodi è presente ma non troppo fitta con una conseguente rappresentazione parziale delle relazioni di dominanza tra tutte le possibili. Questo è coerente con l'uso della regola di orientamento selettiva utilizzata , che tende a collegare tra loro solo profili di target molto simili in termini di dimensione. In altre parole, la density rappresenta una struttura interpretabile biologicamente senza indurre troppi collegamenti ridondanti o poco informativi.


Il numero elevato di sources (circa il 43% dei nodi totali) è coerente con il vincolo imposto sulla differenza tra due insiemi "successivi" il quale limita fortemente la connessione tra insieme molto grandi e molto piccoli. Di conseguenza, molti nodi con set di target grandi non trovano sottoinsiemi “abbastanza vicini” e rimangono senza archi entranti, diventando sorgenti della DAG.

Il numero più basso di sink, suggerisce una struttura a ventaglio più che a imbuto.

Il valore di max_depth suggerisce l'esistenza di strutture a catena del tipo: T1​⊃T2​⊃T3​⊃T4​⊃T5​⊃T6​⊃T7​. Questo suggerisce la possibile presenza di processi di specializzazione progressiva da farmaci “broad-target” a farmaci sempre più selettivi coerenti con la natura del dataset (effetti condizione-specifici o versioni sperimentali dello stesso composto). 

Analizzando il numero elevato di archi si può dire che molti farmaci differiscono per pochi geni e quindi spesso condividono un core comune. 

I risultati ottenuti, sono fortemente influenzati dall'orentetion rule utilizzata. È quindi importante tenere conto delle proprietà introdotte da essa. Tra esse abbiamo la gerarchia locale (e non globale) imposta dal vincolo Δ∣T∣≤3 il quale produce molte gerarchie locali e molte sorgenti. La DAG risulta quindi essere composta da molte gerarchie locali indipendenti, ciascuna con la propria sorgente e pochi livelli di profondità (struttura a più alberi).

## DAG node parameters
I seguenti parametri sono stati calcolati per ciascun nodo:

- in_degree: numero di archi entranti 
- out_degree: numero di archi uscenti 
- degree_ratio: out_degree / (in_degree + 1) 
- topological_level: livello topologico (0 per sorgenti, cresce lungo i predecessori)



Dai file salvati si mostra che:

| Categoria grado | Nodi `in_degree` | Nodi `out_degree` |
|---|---:|---:|
| Tra 1 e 3 | 8 | 14 |
| Uguale a 0 | 35 | 21 |
| Maggiore di 5 | 39 | 47 |

Ovvero, il 25% dei nodi ha in_degree = 0. Questi non sono sottoinsieme di nessun altro e quindi rappresentano i profili massimali. Questo faramci rappresentano una firma comune per molti altri farmaci diversi.

39 farmaci hanno un in_degree > 5. Sono nodi contenuti in molti altri profili. Questi farmaci possono essere di varianti sperimetali o in condizioni differenti dello stesso composto. Proprio come accade per farmaci nel dataset ChG-InterDecagon.

I nodi con outdegree alto rappresentano gli hub gerarchici mentre quelli con outdegree nullo sono nodi specifici che non contengono nessun altro. Di Questa ultima categoria ne fanno parte sia nodi che non si sono mai collegati ad altri (con nessun sottoinsieme) che i due nodi sink mostrati in giallo nel seguente grafico.

Il degree ratio è stato calcolato come segue:
$$
\mathrm{degree\_ratio} = \frac{\mathrm{out\_degree}}{\mathrm{in\_degree} + 1}
$$

Esso è la misura dell'assimetria direzionale di un nodo. È inoltre un indice del ruolo gerarchico di un nodo poichè può distinguere un nodo in:
- Ruolo dominante: degree ratio alto
- Nodo di transizione: degree_ratio ≈ 1
- Nodo foglia/ specializzato: degree ratio basso

I valori calcolati sono:
| Statistica (degree_ratio) | Valore |
|---|---:|
| Media ± SEM | 5.151 ± 1.181 |
| Mediana ± MAD | 0.375 ± 0.375 |
| Minimo | 0.000 |
| Massimo | 31.000 |
| % nodi con ratio < 1 | 73.171% |

Si osserva una distribuzione fortemente sbilanciata in cui la maggiorparte dei nodi hanno un valore << 1 (nodi nodi foglia)
Pochi nodi con valori molto alti dominano la media. 
Osservando il grafico sottostante, si osserva la presenza di diversi ruoli:
- Sorgenti isolate (viola in alto, degree_ratio = 0): Questi sono profili non ordinabili dall'orientation rule scelta
- Sorgenti dominanti (viola in alto, degree_ratio ≫ 1): out degree elevato. Questi sono profili target massimali e quindi radici di grandi sotto gerarchie. 
- Nodi intermedi di smistamento (blu scuro / azzurri centrali): questi hanno:
    - in_degree > 0
    - out_degree > 0
    - degree_ratio ≈ 1
Essi collegano grandi profili a profili più specifici
- Sottoinsiemi condivisi (turchesi, livello medio-basso) e pre-foglie (verdi chiari, penultimo livello)
Essi rappresentano i colli di bottiglia della gerarchia e sono gli ultimi passaggi prima della specializzazione massima.
Di conseguenza con una degree ratio molto basso
- Foglie pure (gialli in fondo, degree_ratio = 0). Essi sono i terminali veri della dag ovvero i profili massimamente specifici. Rappresentano farmaci molto selettivi o condizioni sperimentali estremamente specifiche.

| Categoria                  | # nodi | % sul totale |
| -------------------------- | -----: | -----------: |
| **Sorgenti isolate**       |     18 |       21.95% |
| **Sorgenti dominanti**     |     17 |       20.73% |
| **Nodi intermedi**         |     37 |       45.12% |
| **Sottoinsiemi condivisi** |      6 |        7.32% |
| **Pre-foglie**             |      2 |        2.44% |
| **Foglie**                 |      2 |        2.44% |
| **Totale**                 | **82** |     **100%** |



