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
DIFFERENZA TRA L'IDENTIFICAZIONE DELLE COMUNITÀ E IL CALCOLO DELLA DENSITÀ PER CIASCUNA COMUNITÀ
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
---
