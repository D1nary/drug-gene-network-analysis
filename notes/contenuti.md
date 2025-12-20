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
