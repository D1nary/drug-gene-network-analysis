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

Nella seguente tabella, per motivi di visualizzazione, sono stati riportati i dati ottenuti dall'analisi solo delle comunità con size maggiore o uguale a 10.


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


