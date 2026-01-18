# Forma del dataset
Il dataset di partenza è il file compresso data/ChG-InterDecagon_targets.csv.gz: una tabella a 2 colonne (Drug, Gene) dove ogni riga è un’interazione farmaco‑gene. Se un farmaco è associato a più geni, lo stesso ID drug appare su più righe (una per gene).

## Rappresentazione vettoriale geni
Un farmaco è rappresentato come un vettore perché può essere descritto dall’insieme delle sue interazioni con tutti i geni del dataset, e ogni gene diventa una dimensione dello spazio.

Supponiamo che il dataset abbia **5 geni totali**:

| Gene | G1 | G2 | G3 | G4 | G5 |
|------|----|----|----|----|----|

**Farmaco A** interagisce con **G1, G3, G5**

\[
\vec{A} = (1, 0, 1, 0, 1)
\]

**Farmaco B** interagisce con **G1, G3**

\[
\vec{B} = (1, 0, 1, 0, 0)
\]

Questi vettori **vivono nello stesso spazio** e quindi possono essere **confrontati**.

Formalmente, nel caso dell'esempio, il vettore è un punto (o una freccia dall’origine) nello spazio a 5 dimensioni. Non esiste una rappresentazione grafica diretta in uno spazio a 5 dimensioni.

## Cosine similarity
La cosine similarity misura quanto due vettori “puntano nella stessa direzione”, cioè quanto è simile il loro profilo di interazioni genetiche, indipendentemente dalla quantità assoluta di interazioni.

Con threshold = 0.4, stai dicendo: “Considero due farmaci biologicamente simili solo se condividono una porzione significativa del loro profilo di interazioni genetiche"

Più concretamente:
- i due farmaci agiscono su insiemi di geni in parte sovrapposti
- non è una coincidenza casuale
- la sovrapposizione è sufficientemente strutturata da suggerire un meccanismo comune

0.4 NON significa “40% degli stessi geni”. Significa che la geometria dei loro vettori è abbastanza allineata, cioè:
- colpiscono geni simili
- con pattern simili
- anche se uno dei due ha più target dell’altro

## Modularity
La modularity misura quanto una rete è suddivisibile in comunità reali, ovvero gruppi di nodi più densamente connessi tra loro rispetto a quanto ci si aspetterebbe per caso in una rete casuale con lo stesso grado dei nodi. Più alta è la modularity, più la rete presenta una vera struttura a comunità.

## Pathway come bersaglio
### Cosa significa biologicamente “farmaci che colpiscono gli stessi pathway”?

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

### Perché vengono generate comunità grandi nella similarity network relativamente ad uno stesso pathway?

Perché i pathway rilevanti nelle malattie comuni (tumori, infiammazione, malattie cardiovascolari, metaboliche…) sono:
- molto studiati,
- molto popolati da target farmacologici,
- e quindi molti farmaci diversi finiscono per agire su nodi dello stesso circuito molecolare.

Di conseguenza:
- anche se due farmaci non hanno esattamente gli stessi target,
- se entrambi agiscono su nodi dello stesso pathway,
- la loro cosine similarity sulle liste dei target risulta comunque elevata.

### Come mai parliamo di similarità su pathway e non su singoli target molecolari?
Nel mio dataset i nodi rappresentano geni e farmaci. Se costruisco una rete di similarità drug–drug, in che modo posso giustificare l’ipotesi che le comunità individuate nella rete corrispondano a gruppi di farmaci che agiscono sugli stessi pathway biologici associati a malattie comuni, piuttosto che limitarsi a condividere un singolo target molecolare?

Nel dataset:
- ogni farmaco è collegato a uno o più geni
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

## Modulazione di un sistema biologico da parte di un farmaco
Quando diciamo che un farmaco “modula il sistema biologico” intendiamo come e dove interviene nella rete molecolare cellulare cioè quali proteine (target) colpisce, in quali pathway si inserisce e con che intensità e specificità altera il flusso di segnali biologici. Il sistema biologico non è una catena lineare, ma una rete complessa di interazioni (proteine, pathway, feedback, cross-talk).

## Densità
### DIFFERENZA TRA L'IDENTIFICAZIONE DELLE COMUNITÀ E IL CALCOLO DELLA DENSITÀ PER CIASCUNA COMUNITÀ
La risposta breve è:

> **perché Louvain ti dice _chi_ sta insieme, la density ti dice _quanto bene_ sta insieme.**

Ora ti spiego bene il perché, passo per passo, senza scorciatoie.

1. Cosa fa davvero Louvain (e cosa NON fa)

Il **Louvain method**:

- non cerca gruppi “molto densi” in senso assoluto 
- cerca gruppi che hanno **più connessioni interne di quelle attese per caso**

Formalmente, massimizza la **modularity**:

Q = (edge interni osservati) − (edge interni attesi)

👉 Quindi una comunità può essere:
- poco densa,
- ma comunque una comunità valida,
- se è più connessa internamente che verso l’esterno.

💡 Traduzione biologica:

Louvain trova gruppi **relativamente simili**, non necessariamente **fortemente simili**.


2. Due comunità possono essere entrambe “valide” ma molto diverse

Immagina due comunità trovate da Louvain:

**Comunità A**
- 20 farmaci 
- 150 archi interni 
- density ≈ 0.79 

👉 quasi tutti simili a quasi tutti

**Comunità B**
- 20 farmaci 
- 25 archi interni 
- density ≈ 0.13 

👉 pochi legami, ma più interni che esterni

⚠️ Louvain accetta entrambe come comunità 
ma biologicamente **non raccontano la stessa storia**.

Qui entra in gioco la density.


3. A cosa serve davvero la density (risposta concettuale)

La density risponde a una domanda diversa:

> **Quanto è omogeneo il meccanismo biologico all’interno della comunità?**

**Louvain**
- “Questi farmaci stanno insieme?”
- criterio relativo
- separazione dal resto
- struttura globale

**Density**
- “Quanto sono davvero simili?”
- criterio assoluto
- coesione interna
- qualità biologica del cluster

4. Perché questo è cruciale biologicamente (non solo tecnicamente)

Nel tuo caso (farmaci ↔ target):

**Caso 1 – Density alta**

Farmaci:
- colpiscono gli stessi target 
- o la stessa famiglia proteica 

Biologicamente:
- stessa classe farmacologica 
- effetti simili, tossicità simile 

👉 Qui puoi parlare di **classe funzionale ben definita**

**Caso 2 – Density bassa**

Farmaci:
- condividono solo alcuni target 
- oppure convergono su uno stesso pathway 

Biologicamente:
- meccanismi diversi, stesso effetto finale 
- possibili combinazioni terapeutiche 

👉 Qui **NON** hai una classe omogenea, ma un **modulo funzionale**

Senza la density, questi due casi sembrerebbero uguali, perché Louvain li etichetta entrambi come “comunità”.

5. Il punto chiave (molto importante)

Louvain **NON** garantisce che una comunità sia internamente compatta. 
Garantisce solo che è **più compatta di quanto ci si aspetterebbe per caso**.

La density serve a:
- distinguere:
  - comunità forti 
  - comunità deboli 
- evitare over-interpretazioni biologiche 
- capire quanto fidarti del cluster

### Come mai per comunità più grandi la density è più variabile?
Dal momento che il numero di possibili connessioni interne aumenta quadraticamente con la dimensione della comunità, anche piccole variazioni nei profili dei nodi tendono ad accumularsi in un numero crescente di archi mancanti. Questo effetto è amplificato nelle comunità più grandi, dove l’elevato numero di confronti tra profili rende più probabile l’emergere di differenze, riflettendosi in una maggiore variabilità dei valori di densità.

####Esempio numerico: effetto della dimensione della comunità sulla density

Consideriamo due comunità, una piccola e una grande, e osserviamo come la densità
reagisce alla presenza di archi mancanti.

Comunità piccola
- Numero di nodi: n = 10
- Numero massimo di archi possibili:
  
  n(n−1)/2 = 10·9/2 = 45

- Supponiamo che 5 coppie di nodi non superino la soglia di similarità.
  Gli archi interni osservati sono quindi:

  E = 40

- La densità della comunità è:

  density = 40 / 45 ≈ 0.89

Comunità grande
- Numero di nodi: n = 100
- Numero massimo di archi possibili:

  n(n−1)/2 = 100·99/2 = 4950

- Supponiamo che solo il 5% delle coppie di nodi non superi la soglia di similarità.
  Gli archi osservati sono quindi:

  E = 0.95 · 4950 = 4702

- La densità della comunità è:

  density = 4702 / 4950 ≈ 0.95

Se la frazione di coppie che non supera la soglia aumenta leggermente, ad esempio
al 15%:

- E = 0.85 · 4950 = 4207
- density ≈ 0.85

Osservazione chiave
In una comunità grande, anche piccole percentuali di coppie non connesse producono
variazioni apprezzabili della densità. La densità riflette quindi l'effetto cumulativo
di una modesta eterogeneità interna, mentre nelle comunità piccole la stessa
eterogeneità può rimanere nascosta.


## Weighted degree
Nel file community_parameters.csv, il weighted degree di una comunità non descrive la sua coesione interna, ma il livello complessivo di interazione con le altre comunità.

Questo perché il valore è calcolato sul community graph, un grafo in cui:
- ogni nodo rappresenta una comunità individuata con il metodo di Louvain;
- un arco tra due comunità esiste se nel grafo originale di similarità sono presenti connessioni tra farmaci appartenenti alle due comunità;
- il peso dell’arco tra due comunità è definito come la somma dei pesi di tutte le connessioni farmaco–farmaco che collegano nodi appartenenti alle due comunità (come implementato in network.py).

Di conseguenza, il weighted degree di una comunità è la somma dei pesi di tutti gli archi che la collegano alle altre comunità nel community graph, e quantifica quanto intensamente quella comunità interagisce, in termini di similarità complessiva, con il resto della rete.
In altre parole, questo parametro misura:
- la forza delle connessioni esterne di una comunità,
- il suo possibile ruolo di ponte o hub tra moduli diversi,

## Clustering coefficient
Il clustering coefficient misura il grado di chiusura locale delle comunità, quantificando la tendenza dei farmaci a formare triangoli di similarità. Valori elevati indicano moduli altamente coerenti, in cui farmaci simili a un terzo sono anche reciprocamente simili, mentre valori più bassi suggeriscono strutture eterogenee o la presenza di farmaci con ruolo di ponte tra profili funzionali differenti.

Il clustering coefficient risponde alla domanda: Se il farmaco A è simile a B e a C, quanto è probabile che B e C siano simili tra loro?

Nel nostro caso (archi pesati con Jaccard similarity), si usa una versione pesata che tiene conto anche della forza delle connessioni: triangoli con archi forti contribuiscono più di triangoli con archi deboli


Clustering coefficient alto:
Indica che:
- i farmaci della community sono mutuamente simili
- condividono insiemi di target molto sovrapposti
- la community è funzionalmente coerente

Interpretazione biologica:
- stessa classe farmacologica
- stesso meccanismo d’azione
- possibili duplicati funzionali

Clustering coefficient basso
Indica che:
- i farmaci sono collegati indirettamente
- uno stesso farmaco è simile a molti altri che non sono simili tra loro

Interpretazione biologica:
- farmaci ponte
- polifarmacologia
- target comuni “generici” (hub biologici)



### Interpretazione del clustering coefficient e delle connessioni indirette

**1. Connessioni dirette e indirette**
In una drug–drug similarity network, una connessione diretta tra due farmaci indica che la loro similarità (ad esempio Jaccard) supera una soglia prefissata. Due farmaci possono invece risultare collegati in modo indiretto quando non sono direttamente simili tra loro, ma risultano entrambi simili a un terzo farmaco che funge da intermediario.

**2. Clustering coefficient come misura di transitività**
Il clustering coefficient quantifica la probabilità che, dati tre farmaci A, B e C, se A è simile a B e a C, allora anche B e C siano simili tra loro. In termini topologici, misura la frazione di triplette che risultano chiuse in triangoli.

**3. Clustering basso e triplette aperte**

Un valore basso di clustering coefficient indica che la maggior parte delle triplette è aperta: esistono molte configurazioni del tipo A–B–C senza il collegamento diretto A–C. Questo implica che le relazioni di similarità non sono transitive e che i farmaci sono spesso connessi solo tramite percorsi indiretti.

Analizziamo il seguente esempio

   B
   |
A--X--C
   |
   D

- X è simile a molti farmaci
- B, C, D non sono simili tra loro
- quasi tutte le triplette sono aperte
- triangoli ≈ 0 → clustering basso

Interpretazione:
- B e C sono nella stessa community
- ma non sono direttamente simili
- sono collegati indirettamente tramite X

**4. Community individuate dal Louvain method**
Il Louvain method massimizza la modularity globale della rete e può quindi individuare comunità che sono topologicamente connesse, ma non necessariamente dense o caratterizzate da un’elevata chiusura locale. Di conseguenza, una community può presentare un clustering coefficient basso pur risultando coerente dal punto di vista della modularità.

**5. Interpretazione biologica**
Dal punto di vista biologico, una community con clustering coefficient basso suggerisce la presenza di farmaci con ruolo di “ponte”, caratterizzati da profili target ampi o polifarmacologici. In questo scenario, i farmaci appartenenti alla stessa community non sono tutti direttamente simili tra loro, ma condividono sovrapposizioni funzionali parziali mediate da uno o più nodi centrali.


# Problemi
## Density = 1 nella community con 359 farmaci (threshold impostato a 0.4)
La densità = 1 in una grande comunity è sbaglaita.

---
Facciamo un chiarimento:
Nella similarity network, i nodi sono rappresentati semplicemente da farmaci, e due farmaci sono collegati se la loro similarità supera la soglia di 0.4. Quella similarità, non nasce dal nulla. È calcolata a partire da un profilo bersaglio (vettore binario rappresentante quali geni sono collegati a tale farmaco) che esiste prima della rete. In altre parole:
- Nella rete il nodo è un solo farmaco
- Mentre nel calcolo, il farmaco è rappresentato da un profilo di target.

In altre parole, ogni farmaco è rappresentato da un profilo target
---

Con un threshold di 0.4, nella community analysis, si ha una grande comunità di 359 farmaci con density = 1. In queste condizioni, nella comunity, i 359 farmaci non hanno tutti un profilo diverso ma condividono solo 120 configurazioni distinte di target.

In media (359/120 = 3), ogni profilo viene riutilizzato da ~3 farmaci.

Un singolo profilo bersaglio contiene 158 geni ed è identico per 125 farmaci diversi. Quindi:
- 125 nodi su 359 hanno esattamente lo stesso vettore
- per ogni coppia di questi 125 farmaci:
   - intersezione = 158
   - norma = identica
   - cosine similarity = 1
Considerando questi 125 farmaci abbiamo 7750 coppie.

Oltre ai **125 identici**, gli altri profili:
- differiscono magari per **1–5 geni**
- condividono comunque **150+ target su 158**

La **cosine similarity** tra due vettori così è:
\[
\cos(\theta) \approx \frac{150}{\sqrt{158 \cdot 158}} \approx 0.95\text{–}0.99
\]
**molto sopra la soglia 0.4**

Quindi:
- anche farmaci **non identici** sono comunque **collegati**
- il grafo diventa **densissimo** 

Siccome la density misura la frazione di coppie di farmaci all’interno della community che risultano connesse da un arco nella similarity network, e in questo caso ogni coppia di farmaci supera la soglia di similarità, il sottografo indotto dalla community risulta completamente connesso e la density assume valore 1.

### Primas soluzione provata
Algoritmo TF-IDF:

Nel dataset drug–gene, l’algoritmo TF-IDF pesa ogni gene in modo da rendere simili due farmaci solo se condividono target rari e specifici, riducendo invece l’influenza dei geni molto comuni che interagiscono con molti farmaci.

ESEMPIO:

Farmaco A → geni {G1, G2}

Farmaco B → geni {G1, G3}

Se G1 è un gene molto comune (presente in molti farmaci), il suo peso TF-IDF è basso e contribuisce poco alla similarità.
Se G2 e G3 sono geni rari, il TF-IDF li valorizza: di conseguenza A e B risultano poco simili, perché condividono solo un gene comune e poco informativo.

👉 In questo modo la similarità riflette meccanismi d’azione specifici, non la semplice condivisione di geni “onnipresenti”.

TF-IDF viene descritto separando TF e IDF perché queste due componenti misurano aspetti concettualmente distinti e indipendenti dell’informazione, ed è proprio questa separazione che rende l’algoritmo interpretabile ed efficace.

La TF (Term Frequency) misura la rilevanza locale: risponde alla domanda “quanto questo gene è importante per questo specifico farmaco?”. È una quantità che guarda esclusivamente all’interno del singolo farmaco e non tiene conto del resto del dataset. In termini biologici, la TF quantifica quanto un gene caratterizzi il profilo di interazione di un farmaco.

La IDF (Inverse Document Frequency) misura invece la rilevanza globale: risponde alla domanda “quanto questo gene è informativo nel dataset complessivo?”. Valuta quanto un gene sia raro o specifico rispetto all’insieme di tutti i farmaci. Dal punto di vista biologico, l’IDF penalizza i geni ubiquitari (hub), che compaiono in moltissimi farmaci e quindi hanno scarso potere discriminante, ed enfatizza i geni più selettivi.

Queste due componenti vengono trattate separatamente perché non sono intercambiabili: usare solo la TF porta a una similarità dominata da geni molto comuni, mentre usare solo la IDF scollega il peso del gene dal suo ruolo effettivo nel farmaco specifico. La loro combinazione tramite prodotto,

TF-IDF(d,g)=TF(d,g)⋅IDF(g),

ha un significato preciso: un gene è considerato importante per un farmaco solo se è rilevante localmente per quel farmaco e informativo a livello globale nel dataset. Se una delle due componenti è bassa, il contributo complessivo del gene viene ridotto. In questo modo TF-IDF implementa un principio biologicamente sensato, secondo cui la similarità tra farmaci deve basarsi su target specifici e discriminanti, non sulla condivisione di geni generici e onnipresenti.

### Seconda soluzione provata
Introduzione della jaccard similarity. Non ha risolto il problma. La Jaccard risolve sovrapposizioni spurie dovute a geni hub quando i profili sono diversi ma non può distinguere farmaci con profili identici o quasi identici (nel nostro caso, 125 farmaci hanno lo STESSO set e decine di altri hanno set che differiscono di pochissimi geni). Le possibili soluizoni possono essere: collassare i profili duplicati, porre un vincolo minimo sull'intersezione o rimuovere gli hub cioè i target troppo frequenti.

Ha senso introdurre queste modifiche? dipende dall'obiettivo. NON LE INTRODUCIAMO VISTO LO SCOPO SEMPLICEMENTE DESCRITTIVO DEL PROGETTO




