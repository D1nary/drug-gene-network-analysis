# Similarity network
È stata costruita una similarity network drug-drug attraverso la cosine similarity utilizzando prima un threshold di 0.3 e successivamente uno di 0.4. Sono stati successivamente comparati i risultati (file filtering.json nel programma).

## Threshold 0.3

| Parametro            | Valore    |
| -------------------- | --------- |
| similarity_threshold | 0.3       |
| nodes_removed        | 333       |
| edges_filtered       | 952,210   |
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
| edges_filtered       | 959,117   |
| original_node_count  | 1,774     |
| retained_node_count  | 1,441     |
| potential_edges      | 1,037,520 |

Con un threshold di 0.4 la soglia è più stringente mantenendo solo interazioni più forti rispetto a prima. Vengono filtrati 959,117 mantenendone 78,403 su 1,037,520 possibili coppie. La percentuale delle connessioni che riescono a superare il threshold è del 7.56%

Con tale threshold la rete perde $85.310 - 78.403 = 6.907 $ archi rispetto al caso precedente. Questo numero rappresenta una diminuzione del 
$$
\frac{6.907}{85.310} \approx 8.1 \%
$$
del numero di connessioni rispetto al threshold precedente.

Il numero di nodi rimossi rimane invariato a 333. Ciò indica che questi farmaci non presentano valori di similarità pari o superiori a 0.3 con nessun altro farmaco del dataset, risultando quindi isolati nella similarity network già al threshold più permissivo.

Uno degli scopi di questa analisi, è quello di costruire una community network per ricercare comunità di farmaci con meccanismài d'azione simili. In questo contesto, non avrebbe senso alzare troppo del threshold potrebbe avere l'effetto di tenere solo le relazioni “ovvie” e perdere quelle deboli ma biologicamente interessanti. Infatti anche similarità più deboli rappresentano un caso interessante poichè, due farmaci che condividono pochi target ma “giusti” possono essere candidati per drug repurposing. Un farmaco con profilo di target solo parzialmente simile può comunque avere effetti collaterali simili oppure opportunità di combinazione terapeutica. Soglie troppo rigide vanno contro questa natura “sfumata”.

Le seguenti analisi sono state eseguite considerando un threshold di 0.4.

# Community network analysis
Dopo aver costruito la similarity nework, è stata costruita la relativa comunity network attraverso il louvian method.
Di seguito sono riportati i parametri di Louvian salvati nel file louvian_parameters.json dopo l'esecuzione dell'analisi.

| Parametro              | Valore |
|------------------------|--------|
| method                 | louvain |
| resolution             | 1.000  |
| modularity             | 0.258  |
| min_community_size     | 1      |
| max_community_size     | 359    |
| mean_community_size    | 13.220 |
| median_community_size  | 2.000  |


## Modularity
Una modularity di circa 0.26 indica che esiste una struttura a moduli reale cioè che i farmaci non sono distribuiti a caso, ma si raggruppano in insiemi che condividono pattern di similarità. 

Non si tratta di una modularità “estrema” (ad esempio >0.4–0.5), ma di un valore intermedio che riflette la presenza simultanea di cluster farmacologici ben riconoscibili e di numerosi collegamenti tra comunità diverse. Questo risultato è coerente con la natura intrinsecamente complessa dei sistemi biologici, nei quali molti farmaci presentano profili polifarmacologici e possono agire su target o pathway condivisi tra differenti classi terapeutiche. Di conseguenza, la similarity network drug–drug non mostra una separazione netta tra moduli completamente isolati, ma una struttura modulare interconnessa, realistica e biologicamente plausibile, in cui i “ponti” tra comunità rappresentano potenziali sovrapposizioni funzionali, effetti collaterali comuni o opportunità di combinazione terapeutica e drug repurposing.

## Other parameters
- Median:
Il valore mediano delle comunità risulta essere 2. Questo significa che più della metà delle comunità ha 2 farmaci o meno. Questo suggerisce gruppi di farmaci molto simili fra di loro i quali possono essere varianti strutturali di uno stesso composto oppure farmaci che condividono uno o pochissimi target molto specifici. In altre parole, queste comunità rappresentano farmaci di nicchia con bersagli rari oppure possono essere outliers interessanti per riposizionamento (se un singleton si collega debolmente a una grande comunità, potrebbe condividere qualche pathway con farmaci di un’altra indicazione).

- Mean, max community 
La media è rappresentata da una comunità di circa 13 farmaci. Questo valore insieme a quello della mediana, indicano una distribuzione fortemente sbilanciata con pochi cluster molto grandi e tanti piccoli. 

METTI HISTOGRAMMA DELLE COMUNITÀ IN FUNZIONE DEL NUMERO DI ELEMENTI.

Il valore massimo è rappresentato da una comunità con 359 elementi. In particolare, grandi famiglie farmacologiche come questa rappresentano famiglie che colpiscono la stessa famiglia di proteine (es. molte chinasi, GPCR, recettori nucleari…) oppure gli stessi pathway coinvolti in malattie comuni (es. segnali proliferativi tumorali, infiammazione, ecc.). Oppure



## Inforamazioni biologiche ricavabili da questi valori

Anche in assenza di un’analisi dettagliata della composizione interna di ciascuna comunità, i parametri ottenuti tramite il metodo di Louvain consentono di trarre alcune considerazioni biologiche rilevanti. In primo luogo, la presenza di una modularità diversa da zero indica che la rete di farmaci non è organizzata in modo casuale, ma presenta una struttura modulare, suggerendo l’esistenza di gruppi di farmaci che condividono meccanismi molecolari, bersagli biologici o ambiti terapeutici affini.

La modularità moderata osservata riflette inoltre un’elevata polifarmacologia del sistema: molti farmaci risultano connessi a più comunità attraverso target o pathway condivisi. Questo comportamento è coerente con la possibilità che farmaci appartenenti a classi terapeutiche differenti possano indurre effetti collaterali simili o, viceversa, presentare opportunità di riposizionamento farmacologico, qualora condividano un sottoinsieme rilevante di bersagli biologici.

La distribuzione delle dimensioni delle comunità, caratterizzata dalla presenza di numerosi moduli di piccola dimensione e di nodi isolati, suggerisce l’esistenza di farmaci con profili di target rari o altamente specifici. Tali regioni del network possono rappresentare meccanismi d’azione poco ridondanti o ancora scarsamente esplorati, e risultano quindi potenzialmente interessanti dal punto di vista della scoperta e dello sviluppo di nuovi farmaci.

Al contrario, la presenza di pochi moduli di grandi dimensioni evidenzia aree del network associate a target e pathway ampiamente studiati e frequentemente sfruttati nello sviluppo farmacologico, come avviene tipicamente per alcuni pathway oncogeni o infiammatori ben caratterizzati. In queste comunità si osserva un’elevata sovrapposizione farmacologica, che può essere utile per confrontare profili di efficacia e tossicità tra farmaci simili, nonché per individuare strategie di combinazione terapeutica all’interno della stessa comunità o tra comunità strettamente correlate.

## Communities analisys
Le singole comunità vengono analizzate nel file community_parameters.csv. Il file è organizzato in colonne ciascuna con una caratteristica della comunità:
- community id: Identificativo della comunità
- size 
- degree
- weighted degree
- clustering coefficient

### Size
Biologicamente rappresenta quanti farmaci condividono un profilo di target simile oppure, in altre parole, la popolarità di un certo spazio farmacologico. Comunità grandi rappresentano target molto studiati e sfruttati, pathway centrali e oppure alta ridondanza terapeutica. Al contrario, comunità piccole si manifestano quando siamo in presenza di target rari o specifici. Mentre comunità grandi rappresentano meccanismi già esplorati (poichè diversi farmaci agiscono sullo stesso bersaglio o sullo stesso pattern), le comunità piccole individuano target non molto esplorati farmacologicamente oppure target presenti in pathway meno studiati e quindi famiglie proteiche non ancora saturate da altre molecole. Inoltre, se un farmaco non clusterizza con gli altri significa che non condivide un numero sufficiente di target da superare il threshold, non si inserisce nei pathway delle classi terapeutiche classiche e che potrebbe modulare il sistema biologico in modo differente. Cioè, essendo il sistema biologico una rete complessa di interazioni (proteine, pathway, feedback, ecc), un farmaco appartenente ad una piccola comunity, altera la rete in un punto con una dinamica diversa dalla maggiorparte dei farmaci esistenti. Questo è molto utile in farmacologia. 

### Density

AGGIUNGI FORMULA DELLA DENSITÀ

Mentre il Louvian method identifica comunità di farmaci strutturamente disinte, il calcolo della densità permette di valutare il grado di coerenza interna di ciascun modulo permettendo di distinguere vere classi farmacologiche omogenee da insiemi di farmaci che sono stati raggruppati nella stessa comunità non perchè siano tutti altamente simili tra loro, ma perchè condividono alcune funzioni biologiche o bersagli.

---
Esempio concettuale
Immagina una comunità che include:
- inibitori diretti di una chinasi A;
- modulatori upstream che regolano l’attivazione di A;
- farmaci che colpiscono un pathway parallelo ma convergente.
Questi farmaci:
- non condividono tutti gli stessi target,
- ma agiscono sulla stessa funzione biologica finale (es. proliferazione cellulare).
---

Si osserva che comunità di piccole dimensioni mostrino spesso una density elevata, mentre comunità più grandi presentino density più basse, è un comportamento atteso sia dal punto di vista matematico dei network sia da quello biologico dei farmaci.

Nei moduli piccoli, il numero di connessioni possibili è molto limitato: basta che pochi nodi siano tutti connessi tra loro perché la density risulti elevata, spesso prossima a 1. Questo riflette una similarità molto forte tra i farmaci del gruppo (ad esempio condivisione quasi completa dei target), ma tali valori sono poco robusti dal punto di vista statistico, perché fortemente influenzati dal basso numero di nodi.

Al contrario, nelle comunità più grandi il numero di connessioni possibili cresce quadraticamente con la size. In questi casi è biologicamente e matematicamente improbabile che tutti i farmaci siano simili a tutti gli altri: la density tende quindi a diminuire, pur rimanendo sufficientemente alta da indicare una coerenza funzionale del modulo. Queste density moderate sono più affidabili e descrivono classi farmacologiche ampie, caratterizzate da meccanismi d’azione o pathway comuni ma con una naturale eterogeneità interna.

Di seguito riportiamo in tabella solo i valori di densità per comunità con un numero maggiore o uguale a 15

Community id	Size	Density
Community_1	359	1.000
Community_68	225	0.205
Community_67	211	0.096
Community_17	202	0.055
Community_37	135	0.494
Community_49	18	0.490
Community_81	16	0.667

