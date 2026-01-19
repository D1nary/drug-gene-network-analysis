# Descrizione del dataset

# Similarity network
È stata costruita una similarity network drug-drug attraverso la Jaccard similarity utilizzando prima un threshold di 0.3 e successivamente uno di 0.4. Sono stati successivamente comparati i risultati (file filtering.json nel programma). La Jaccard similarity analizza la sovrapposizione relativa tra due farmaci
$$
J(A,B) = \frac{|A \cap B|}{|A \cup B|}
$$

Ogni nodo nel grafico della similarity network è rappresenta un famaco, il quale è caratterizzato da un profilo bersaglio ovvero da un vettore contenente tutti i farmaci al quale tale farmaco è legato.
```bash
{
  "bipartite": "drug",
  "original_id": 12345,
  "targets": [1017, 1956, 7422]
}
```

La Jsccard similatiry ci sice che due farmaci risultano similli solo se condividono molti target e hanno pochi target diversi, ovvero se hanno un profilo di target simile. È stata scelta la jaccard similarity poichè essa lavora nativamente su insiemi e penalizza farmaci con tanti target non condivisi e overlap piccoli ma casuali. 




## Threshold 0.3

| Parametro            | Valore    |
| -------------------- | --------- |
| similarity_threshold | 0.3       |
| nodes_removed        | 333       |
| edges_filtered       | 961,999   |
| original_node_count  | 1,774     |
| retained_node_count  | 1,441     |
| potential_edges      | 1,037,520 |

Questo valore di threshold indica che un arco (drug–drug) è presente nella rete solo se la similarità tra due farmaci è ≥ 0.3. Essa è una soglia moderata non troppo permissiva ma nemmeno troppo stringente. Mantiene connessioni con similarità medio–bassa, quindi preserva una rete relativamente densa rispetto a soglie più stringenti. Genera una rete connessa ma non troppo, utile per il community detection. Riduce il rumore eliminando similarità assolutamente deboli. 

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

Uno degli scopi di questa analisi, è quello di costruire una community network per ricercare comunità di farmaci con meccanismài d'azione simili. In questo contesto, non avrebbe senso alzare troppo del threshold potrebbe avere l'effetto di tenere solo le relazioni “ovvie” e perdere quelle deboli ma biologicamente interessanti. Infatti anche similarità più deboli rappresentano un caso interessante poichè, due farmaci che condividono pochi target ma “giusti” possono essere candidati per drug repurposing. Un farmaco con profilo di target solo parzialmente simile può comunque avere effetti collaterali simili oppure opportunità di combinazione terapeutica. Soglie troppo rigide vanno contro questa natura “sfumata”.

Le seguenti analisi sono state eseguite considerando un threshold di 0.4.

# Community network analysis
Dopo aver costruito la similarity nework, è stata costruita la relativa comunity network attraverso il louvian method.
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
Il valore mediano delle comunità risulta essere 2. Questo significa che più della metà delle comunità ha 2 farmaci o meno. Questo suggerisce gruppi di farmaci molto simili fra di loro i quali possono essere varianti strutturali di uno stesso composto oppure farmaci che condividono uno o pochissimi target molto specifici. In altre parole, queste comunità rappresentano farmaci di nicchia con bersagli rari oppure possono essere outliers interessanti per riposizionamento (se un singleton si collega debolmente a una grande comunità, potrebbe condividere qualche pathway con farmaci di un’altra indicazione).

- Mean, max community 
La media è rappresentata da una comunità di circa 5 farmaci. Questo valore insieme a quello della mediana, indicano una distribuzione fortemente sbilanciata con pochi cluster molto grandi e tanti piccoli. 

METTI HISTOGRAMMA DELLE COMUNITÀ IN FUNZIONE DEL NUMERO DI ELEMENTI.

Il valore massimo è rappresentato da una comunità con 359 elementi. In particolare, grandi famiglie farmacologiche come questa rappresentano famiglie che colpiscono la stessa famiglia di proteine (es. molte chinasi, GPCR, recettori nucleari…) oppure gli stessi pathway coinvolti in malattie comuni (es. segnali proliferativi tumorali, infiammazione, ecc.). Nella sezione successiva approffondiamo meglio questo valore nelle sezioni successive.



## Inforamazioni biologiche ricavabili da questi valori

Anche in assenza di un’analisi dettagliata della composizione interna di ciascuna comunità, i parametri ottenuti tramite il metodo di Louvain consentono di trarre alcune considerazioni biologiche rilevanti. In primo luogo, la presenza di una modularità diversa da zero indica che la rete di farmaci non è organizzata in modo casuale, ma presenta una struttura modulare, suggerendo l’esistenza di gruppi di farmaci che condividono meccanismi molecolari, bersagli biologici o ambiti terapeutici affini.

La modularità moderata osservata riflette inoltre un’elevata polifarmacologia del sistema: molti farmaci risultano connessi a più comunità attraverso target o pathway condivisi. Questo comportamento è coerente con la possibilità che farmaci appartenenti a classi terapeutiche differenti possano indurre effetti collaterali simili o, viceversa, presentare opportunità di riposizionamento farmacologico, qualora condividano un sottoinsieme rilevante di bersagli biologici.

La distribuzione delle dimensioni delle comunità, caratterizzata dalla presenza di numerosi moduli di piccola dimensione e di nodi isolati, suggerisce l’esistenza di farmaci con profili di target rari o altamente specifici. Tali regioni del network possono rappresentare meccanismi d’azione poco ridondanti o ancora scarsamente esplorati, e risultano quindi potenzialmente interessanti dal punto di vista della scoperta e dello sviluppo di nuovi farmaci.

Al contrario, la presenza di pochi moduli di grandi dimensioni può evidenziare aree del network associate a target e pathway ampiamente studiati e frequentemente sfruttati nello sviluppo farmacologico, come avviene tipicamente per alcuni pathway oncogeni o infiammatori ben caratterizzati. In queste comunità si osserva un’elevata sovrapposizione farmacologica, che può essere utile per confrontare profili di efficacia e tossicità tra farmaci simili, nonché per individuare strategie di combinazione terapeutica all’interno della stessa comunità o tra comunità strettamente correlate.

## Communities analisys
Le singole comunità vengono analizzate nel file community_parameters.csv. Il file è organizzato in colonne ciascuna con una caratteristica della comunità:
- community id: Identificativo della comunità
- size: Dimensione della community
- degree: Numero di arhci del nodo
- weighted degree: Somma dei pesi degli archi incidenti su un nodo
- clustering coefficient: Grado di chiusura locale delle comunità

### Size
Biologicamente rappresenta quanti farmaci condividono un profilo di target simile oppure, in altre parole, la popolarità di un certo spazio farmacologico. Comunità grandi rappresentano target molto studiati e sfruttati, pathway centrali e oppure alta ridondanza terapeutica. Al contrario, comunità piccole si manifestano quando siamo in presenza di target rari o specifici. Mentre comunità grandi rappresentano meccanismi già esplorati (poichè diversi farmaci agiscono sullo stesso bersaglio o sullo stesso pattern), le comunità piccole individuano target non molto esplorati farmacologicamente oppure target presenti in pathway meno studiati e quindi famiglie proteiche non ancora saturate da altre molecole. Inoltre, se un farmaco non clusterizza con gli altri significa che non condivide un numero sufficiente di target da superare il threshold, non si inserisce nei pathway delle classi terapeutiche classiche o che potrebbe modulare il sistema biologico in modo differente. Cioè, essendo il sistema biologico una rete complessa di interazioni (proteine, pathway, feedback, ecc), un farmaco appartenente ad una piccola comunity, altera la rete in un punto con una dinamica diversa dalla maggiorparte dei farmaci esistenti. Questo è molto utile in farmacologia. 

### Density

Mentre il Louvian method identifica comunità di farmaci strutturamente disinte, il calcolo della densità permette di valutare il grado di coerenza interna di ciascun modulo permettendo di distinguere vere classi farmacologiche omogenee da insiemi di farmaci che sono stati raggruppati nella stessa comunità non perchè siano tutti altamente simili tra loro, ma perchè condividono alcune funzioni biologiche o bersagli. La densità è stata calcolata con:

$$
\text{density} = \frac{E_{\text{int}}}{\frac{n(n-1)}{2}}
$$
dove:

- density: frazione di connessioni interne presenti rispetto al massimo possibile nella comunità.
- E_int: numero di archi interni alla comunità (solo tra nodi della stessa comunità).
- n: numero di nodi nella comunità (size).
- n(n-1)/2: numero massimo di archi possibili in una comunità non orientata senza self‑loop.

Nei moduli piccoli, il numero di connessioni possibili è molto limitato, infatti basta che pochi nodi siano tutti connessi tra loro perché la density risulti elevata, spesso prossima a 1. Questo riflette una similarità molto forte tra i farmaci del gruppo (ad esempio condivisione quasi completa dei target), ma tali valori sono poco robusti dal punto di vista statistico, perché fortemente influenzati dal basso numero di nodi.

Un'ulteriore osservazione è che comunità medio-grandi, presentano valori di densità variabile. Analizzando le comunità presenti nella tabella con size compresa tra 39 e 72 si osservano densità molto diverse. Ovvero comprese tra 0.324 e 0.913. 

Le cause di questo comportamento possono essere sia matematiche che biologiche. 

Le cause matematiche sono da attribuirsi a come viene calcolata la density (formula precedente). Dal momento che il numero di possibili connessioni interne aumenta quadraticamente con la dimensione della comunità, anche piccole variazioni nei profili dei nodi tendono ad accumularsi in un numero crescente di archi mancanti. Questo effetto è amplificato nelle comunità più grandi, dove l’elevato numero di confronti tra profili rende più probabile l’emergere di differenze, riflettendosi in una maggiore variabilità dei valori di densità.

Inoltre, la presenza o l'assenza di un arco dipende da una soglia (0.4). Quindi, se due nodi sono appena sotto la soglia l'arco nella similarity network è assente mentre, se la similarity è appena sopra la soglia, l'arco è presente. Questo introduce variabilità artificiale e maggiore dispersione dei peorfili di densità. 

Da un punto di vista biologico, comunità grandi spesso aggregano famiglie farmacologiche ampie, pathway complessi o target parzialmente sovrapposti. Di conseguenza vengono creati cluster di più grandi dimensioni ma non completamente connessi. Inoltre comunità grandi sono spesso caratterizzate da una struttura core-periferia. Il core è molto denso mentre le periferie lo sono meno. Tale struttura porta ad una diminuzione della densità.

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

Uno delle possibili cause di questo valore alto di density per una comunity cosi grande, può essere che la Jaccard diventa molto permissiva quando i set sono grandi e l'intersezione tra due set è molto ampia. Infatti avendo set di grandi dimensioni, anche con decine di geni diversi, la similarità resta alta. Avendo poi un threshold sufficientemente permissivo, come nel nostro caso, il valore di $E_{\text{int}} \approx \frac{n(n-1)}{2}$ con conseguente density $\approx 1$. In questo caso non siamo difronte ad una ad una causa biologica dei farmaci ma piuttosto ad una causa dovuta da una proprietà matematica della Jaccard similarity. 


Un'ulteriore causa del valore quasi unitario di questa community, potrebbe provenire dalla natura intrinseca del dataset. Infatti esso, aggregando screening diversi producendo il risultato che farmaci testati negli stessi screening, sugli stessi pannelli genici con risultati simili diventano vettorialmente indistinguibili. Questo produce sottografi altamente densi nella rete di similarità (profili identici -> tutte le connessioni di similarità sono presenti -> con threshold basso tutti gli archi sono presenti -> density = 1)

--- 
Farmaci testati negli stessi screening:

- screening 1: testa A e B su geni {G1, G2, G3}
- screening 2: testa A e B sugli stessi geni {G1, G2, G3}
- screening 3: testa A e B sugli stessi geni {G1, G2, G3}

screening diversi = esperimenti diversi (utilizzando stessi farmaci e geni target)
---

All'interno del dataset quindi, molti farmaci sono distinti per nome o contesto, ma non per profilo genetico producendo nodi distinti con vettori profilo identici e archi completi. Ogni farmaco è un nodo distinto identificato da un ID chimico con possibilità che provenga da fonti diverse, screening diversi o contesti sperimentali divesti. In altre parole due molecole possono essere chimicamente diverse o chimicamente molto simili o lo stesso composto annotato in contesti divesi. Si ha la produzione di un clique artificilae dal punto di vista topologico. 

Biologicamente, senza nessuna analisy più prodonda, questa comunità rappresenta una classe di composti con meccanismo d'azione quasi identico oppure il targeting di uno stesso grande pathway o complesso genico. Per capire meglio la natura di questa communiti andrebbero condotte analisi più approfondite. 

È inoltre presente una comunità con size 20 e density 1. In questo caso, il numero di archi possibili sono 190 e, similmente al caso sprecedente, il numero di archi presenti nella comunity è 190. In partiolare 19 farmaci su 20 hanno lo stesso identico profilo target il quale contiene solo 2 geni. Questo caso rappresenta un clique più banale del precedente poichè il profilo è piccolissimo ed essendo in preseza di una replica quasi perfetta. 

Biologicamente la community può essere interpretata come lo stesso evento biologico ripetuto in condizioni sperimentali diverse o come replica tecnica. Non rappresenta una famiglia farmacologica ma è una ridondanza del dataset

### Clustering coefficient
Per ciascuna comunity, è stato calcolato il clustering coefficient.
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

Globalmente di può osservare che le comunità hanno valori mediamente alti di clustering coefficient. Questo suggerisce che le comunità individuate non presentano una struttura lineare o “a catena”, in cui i farmaci risultano simili solo a pochi vicini immediati, ma piuttosto costituiscono moduli fortemente coesi. Cioè, se un farmaco è simile ad altri farmaci all'interno della comuintà, è molto porbabile che esso sia siamile ad altri componenti della stessa. In questi casi i famaci hanno profili bersaglio molto sovrapposti.

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
| Community | Size | Density | Clustering |
| --------- | ---- | ------- | ---------- |
| 188       | 72   | 0.324   | 0.616      |
| 203       | 23   | 0.482   | 0.513      |

Queste comunità mostrano poche connessioni globali ma connessioni locali ben strutturate. In altre parole, la comunità non è un blocco compatto ma un insieme di cluster locali collegati indirettamente. Biologicamente, siamo in presenza di profili target eterogenei e farmci che condividono solo alcune componenti funzionali o pathway comuni. Quindi possibile presenza di pathway parzialmente condivisi farmaci ponte tra meccanismi diversi. In questo caso il clustering coefficient è molto informativo poichè rivela coerenza locale della community nascosta dalla bassa density.


Un'ultimo caso che si può notare dalla tabella precedente sono comunità con density $\approx 1$ e clustering alto. 
| Community | Size | Density | Clustering |
| --------- | ---- | ------- | ---------- |
| 5         | 359  | 1.000   | 0.903      |
| 81        | 20   | 1.000   | 0.964      |
In queste comunità il clustering coefficient non aggiunge nuova informazione rispetto alla density, ma rafforza l’evidenza di omogeneità estrema.

## Weighted degree
Per descrivere il livello complessivo di interazione di una comunità con le altre è stato calcolato, per ciascuna comunittà il parametro weighed degree. Esso è calcolato come la somma dei pesi degli archi inter-comunità (similarità Jaccard). 

Solo le comunità 22 e 78 risultano essere connesse tra di loro con un $\text{weighted degree} = 12.99$. Il fatto che solo due comunità risultano collegate è coerente con le scelte fatte per l'individuazzione della similarità e delle comunità (jaccsrd similarity e Lousvian  method). 

Siccome, un arco tra due comunità esiste solo se esistono farmaci appartenenti a comunità diverse ma connessi e che quindi condividono una porzione significativa di target genici (Jaccard maggiore uguale a 0.4), le due comunità rappresentano moduli farmacologici distinti ma non indipendenti. Le due comunità possono essere associate a pathway diversi ma interconnessi con una sovrapposizione non trascurabile di geni chiave. Alcuni esempi possono essere pathway di segnalazione che convergono su nodi comuni, pathway a monte / a valle oppure pathway cellulari trasversali. I farmaci che connettono le due comunità vengono chiamati farmaci ponte. Essi sono multi-target colpendo geni presenti in entrambe le comunità avendo inoltre, profili più "larghi" rispetto a farmaci altamente specifici. 

Il valore di weighted degree legato alla connessione di queste due comunità è di 12.99.

# Co-occurence network



