# Similarity network
Costruito una similarity network drug-drug attraverso la cosine similarity con un threshold impostato a 0,3 quindi considerando farmaci che condividono lo stesso target in modo significativo

È stata costruita una similarity nework attraverso la cosine similarity. Prima è stato utilizzato un threshold di 0.3 e successivamente uno di 0.4 e sono stati comparati i risultati (file filtering.json nel programma). 

## Threshold 0.3

| Parametro            | Valore    |
| -------------------- | --------- |
| similarity_threshold | 0.3       |
| nodes_removed        | 333       |
| edges_filtered       | 952,210   |
| original_node_count  | 1,774     |
| retained_node_count  | 1,441     |
| potential_edges      | 1,037,520 |

Questo valore di threshold indica che un arco (drug–drug) è presente nella rete solo se la similarità tra due farmaci è ≥ 0.3. Essa è una soglia moderata non troppo permissiva ma nemmeno troppo stringente. Mantiene connessioni con similarità medio–bassa, quindi preserva una rete relativamente densa. Genera una rete connessa ma non troppo, utile per il community detection. Riduce il rumore eliminando similarità assolutamente deboli. 

Siccome il numero di nodi originali (original_node_count) è 1774 e il numero di nodi "sopreavvisuti" (retained_node_count) è 1441, si osserva che, dopo il filtraggio, la rete conserva l'81,2% dei farmaci. Oltre 1400 farmaci sono connessi tra di loro significa, con questo threshold, il dataset contiene informazioni ridondanti sui target farmacologici

Il parametro potential edges rappresenta il totale delle coppie drug-drug possibili prima del filtraggio mentre edges filtered indica quante coppie sono state eliminate perché avevano similarità < 0.3. Sottraendo questi due valori, si ottiene che la rete finale ha 85310 archi. Quindi solo circa l’8.22% delle possibili connessioni supera la soglia, indicando una rete relativamente sparsa. Questa sparsità è coerente con reti di similarità farmacologica poichè, in genere, pochi farmaci sono veramente simili.



## Threshold 0.4
| Parametro            | Valore    |di tenere solo le relazioni “ovvie” e perdere quelle deboli ma biologicamente interessanti.
| -------------------- | --------- |
| similarity_threshold | 0.4       |
| nodes_removed        | 333       |
| edges_filtered       | 959,117   |
| original_node_count  | 1,774     |
| retained_node_count  | 1,441     |
| potential_edges      | 1,037,520 |

Con un threshold di 0.4 la soglia è più stringente mantenendo solo interazioni più forti rispetto a prima. Vengono filtrati 959,117 mantenendone 78,403. La percentuale delle connessioni che riescono a superare il threshold è del 7.56%

Con tale threshold la rete perde $85.310 - 78.403 = 6.907 $ archi rispetto al caso precedente. Questo numero rappresenta una diminuzione del 
$$
\frac{6.907}{85.310} \approx 8.1 \%
$$
del numero di connessioni rispetto al threshold precedente.

Il numero di nodi rimossi rimane invariato a 333. Questo significa che questi 333 farmaci avevano similarità basse con qualunque altro farmaco già prima.

Uno degli scopi di questa analisi, è quello di costruire una community network per ricercare comunità di farmaci con meccanismài d'azione simili. In questo contesto, non avrebbe senso alzare troppo del threshold potrebbe avere l'effetto di tenere solo le relazioni “ovvie” e perdere quelle deboli ma biologicamente interessanti. Infatti anche la similarità ovvia rappresenta un caso interessante poichè, due farmaci che condividono pochi target ma “giusti” possono essere candidati per drug repurposing. Un farmaco con profilo di target solo parzialmente simile può comunque avere effetti collaterali simili oppure opportunità di combinazione terapeutica. Soglie troppo rigide vanno contro questa natura “sfumata”.

Le seguenti analisi sono state eseguite considerando un threshold di 0.4.

# Community network analysis
Dopo aver costruito la similarity nework, è stata costruita la relativa comunity network attraverso il louvian method.
Di seguito sono riportati i parametri di Louvian creati nel file louvian_parameters.json

---
NON INSERIRE: La modularity misura quanto una rete è suddivisibile in comunità reali, ovvero gruppi di nodi più densamente connessi tra loro rispetto a quanto ci si aspetterebbe per caso in una rete casuale con lo stesso grado dei nodi. Più alta è la modularity, più la rete presenta una vera struttura a comunità.
---

- method: "louvain"
- resolution: 1.0
- modularity: 0.25754534108516103
- min_community_size: 1
- max_community_size: 359
- mean_community_size: 13.220183486238533
- median_community_size: 2.0

## Modularity
Una modularity di circa 0.26 indica che esiste una struttura a moduli reale: i farmaci non sono distribuiti a caso, ma si raggruppano in insiemi che condividono pattern di similarità. 

Non è una modularità “estrema” (tipo >0.4–0.5), quindi:

ci sono molti legami anche tra comunità diverse (farmaci polifarmacologici, target condivisi tra classi diverse, ecc.);

è coerente con un sistema biologico complesso, dove gli stessi bersagli o pathway possono essere presi di mira da più classi di farmaci.

In termini biologici:
la tua similarity network drug-drug mostra cluster farmacologici riconoscibili, ma con parecchi “ponti” tra cluster, come ci si aspetta da un network di farmaci/target realistico.

Non è una modularità “estrema” quindi ci sono molti legami anche tra comunità diverse (farmaci polifarmacologici, target condivisi tra classi diverse, ecc.)


In altri termini, la similarity network drug-drug mostra cluster farmacologici riconoscibili, ma con parecchi “ponti” tra cluster, come ci si aspetta da un network di farmaci/target realistico.

## Other parameters
- Median:
Il valore mediano delle comunità risulta essere 2. Questo significa che più della metà delle comunità ha 2 farmaci o meno. Questo suggerisce gruppi di farmaci molto simili fra di loro i quali possono essere varianti strutturali di uno stesso composto oppure farmaci che condividono uno o pochissimi target molto specifici. In altre parole, queste comunità rappresentano farmaci di nicchia con bersagli rari oppure possono essere outliers interessanti per riposizionamento (se un singleton si collega debolmente a una grande comunità, potrebbe condividere qualche pathway con farmaci di un’altra indicazione).

- Mean, max community 
La media è rappresentata da una comunità di circa 13 farmaci. Questo valore insieme a quello della mediana, indicano una distribuzione fortemente sbilanciata con pochi cluster molto grandi e tanti piccoli. 
METTI HISTOGRAMMA DELLE COMUNITÀ IN FUNZIONE DEL NUMERO DI ELEMENTI.
Il valore massimo è rappresentato da una comunità con 359 elementi. In particolare, grandi famiglie farmacologiche come questa rappresentano famiglie che colpiscono la stessa famiglia di proteine (es. molte chinasi, GPCR, recettori nucleari…) oppure gli stessi pathway coinvolti in malattie comuni (es. segnali proliferativi tumorali, infiammazione, ecc.). Oppure

---
Cosa significa biologicamente “farmaci che colpiscono gli stessi pathway”
Quando due farmaci hanno profili di similarità simili (nella tua rete, cosine similarity ≥ 0.4 basata sui target), questo non significa necessariamente che colpiscano gli stessi geni identici, ma spesso che colpiscono set di geni appartenenti allo stesso pathway.

Esempi:
- pathway proliferativi (MAPK, PI3K/AKT)
- pathway infiammatori (NF-κB, JAK/STAT)
- pathway metabolici (mTOR)
- pathway apoptotici

Un pathway non è un gruppo casuale di proteine, ma una catena di azioni in cui ogni molecola:
- riceve un segnale,
- lo elabora,
- e lo trasmette alla molecola successiva.
È coordinata perché:
- le reazioni avvengono in ordine,
- con tempi specifici,
- e con relazioni di attivazione/inibizione.

Un pathway è costruito da interazioni tra:
- recettori (che percepiscono un segnale esterno),
- enzimi (che modificano altre proteine mediante fosforilazione, taglio, ecc.),
- proteine adattatrici (che collegano diversi moduli),
- trasduttori di segnale (che propagano l’informazione),
- fattori di trascrizione (che cambiano l’espressione genica).

L’obiettivo finale è una risposta cellulare definita.
Esempi di funzioni finali:
- proliferazione cellulare
- morte cellulare (apoptosi)
- risposta infiammatoria
- risposta allo stress
- regolazione del metabolismo
- differenziamento cellulare

Perché vengono generate comunità grandi nella similarity network relativamente ad uno stesso pathway?

Perché i pathway rilevanti nelle malattie comuni (tumori, infiammazione, malattie cardiovascolari, metaboliche…) sono:
- molto studiati,
- molto popolati da target farmacologici,
- e quindi molti farmaci diversi finiscono per agire su nodi dello stesso circuito molecolare.
Di conseguenza:
- anche se due farmaci non hanno esattamente gli stessi target,
- se entrambi agiscono su nodi dello stesso pathway,
- la loro cosine similarity sulle liste dei target risulta comunque elevata.

Nel mio dataset i nodi rappresentano geni e farmaci. Se costruisco una rete di similarità drug–drug, in che modo posso giustificare l’ipotesi che le comunità individuate nella rete corrispondano a gruppi di farmaci che agiscono sugli stessi pathway biologici associati a malattie comuni, piuttosto che limitarsi a condividere un singolo target molecolare?

Nel dataset:
- ogni farmaco è collegato a più geni
- nella similarity network drug-drug, due farmaci sono connessi se condividono un pattern simile di target

farmaci con profili di target simili avranno effetti biologici simili, perché “agitano” gli stessi pezzi del sistema cellulare.

Dal punto di vista biologico:
- Molti farmaci non agiscono su un solo bersaglio, ma su più proteine dello stesso pathway o di pathway fortemente collegati (polifarmacologia funzionale).
- Le malattie comuni (tumori, infiammazione, malattie cardiovascolari, metaboliche) sono spesso dovute a disregolazioni di interi pathway, non di un singolo gene isolato.
- La progettazione (o selezione) dei farmaci tende a convergere su nodi e rami critici di questi pathway (recettori, chinasi, fattori di trascrizione…).

Se la similarità drug–drug è basata su interi vettori di target (cosine sui vettori gene-farmaco), allora:
- ogni arco non dice solo “condividono un target”,
- ma “hanno un pattern globale di target simile” (anche molti target diversi ma nello stesso contesto). 

POSSIBILE ANALISI
Per ogni comunità Louvain:
1. prendi l’unione dei geni target dei farmaci in quel modulo;
2. fai un’analisi di arricchimento (enrichment) su:
   - KEGG
   - Reactome
   - GO Biological Process
3. confronti la lista di target del modulo con:
   - tutti i target della rete, oppure
   - tutti i geni del genoma (a seconda dell’approccio).
Se una comunità:
- ha un arricchimento forte per pochi pathway specifici (es. PI3K–AKT, JAK–STAT, NF-κB, mTOR, angiogenesi…),
- con p-value corretti (FDR) molto bassi,
allora puoi dire con buona confidenza:
“Questa comunità non è solo un gruppo di farmaci che condividono qualche target,
ma è un modulo di farmaci che nel complesso convergono su questo pathway biologico specifico”.

Se più community mostrano arricchimenti in pathway classici di malattie comuni (tumori, infiammazione, diabete, ecc.), hai un’evidenza diretta che la modularità della rete riflette la struttura dei pathway di malattia.
---

## Riassumendo
### Tipi di informazioni biologiche che puoi derivare già da questi numeri.

Anche senza ancora guardare il contenuto di ogni comunità, i parametri di Louvain ti permettono di dire che:

1. La rete di farmaci è modularmente organizzata
   - → Esistono gruppi di farmaci che probabilmente condividono meccanismi molecolari o destinazioni terapeutiche.

2. La modularità moderata riflette un’elevata polifarmacologia
   - → Molti farmaci hanno target condivisi con più comunità, quindi:
      - potenziali effetti collaterali simili tra comunità diverse;
      - possibilità di riposizionamento (un farmaco di tumore che si avvicina a una comunità di anti-infiammatori, ecc.).

3. La presenza di molte comunità piccole e singleton
→ Indica:
   - bersagli rari o specifici;
   - possibili meccanismi “di frontiera”;
   - regioni del network ancora poco esplorate o poco ridondanti → interessanti per sviluppare nuovi farmaci.
4. La presenza di pochi moduli molto grandi
→ Evidenzia:
   - target e pathway molto “popolari” nello sviluppo di farmaci (es. pathway oncogeni classici);sw
   - aree dove esiste sovrapposizione farmacologica alta → utile per:
      - confrontare efficacia e tossicità;
      - scegliere combinazioni terapeutiche all’interno della stessa comunità o tra comunità vicine.
      
## Communities analisys
Le singole comunità vengono analizzate nel file community_parameters.csv. Il file è organizzato in colonne ciascuna con una caratteristica della comunità:
- community id: Identificativo della comunità
- size 
- degree
- weighted degree
- clustering coefficient

### Size
Biologicamente rappresenta quanti farmaci condividono un profilo di target simile oppure, in altre parole, la popolarità di un certo spazio farmacologico. Comunità grandi rappresentano target molto studiati e sfruttati, pathway centrali e oppure alta ridondanza terapeutica. Al contrario, comunità piccole si manifestano quando siamo in presenza di target rari o specifici. Mentre comunità grandi rappresentano meccanismi già esplorati (poichè diversi farmaci agiscono sullo stesso bersaglio o sullo stesso pattern), le comunità piccole individuano target non molto esplorati farmacologicamente oppure target presenti in pathway meno studiati e quindi famiglie proteiche non ancora saturate da altre molecole. Inoltre, se un farmaco non clusterizza con gli altri significa che
