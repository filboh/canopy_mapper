# CanopyMapper

Script per la **Processing Toolbox di QGIS** (3.x) che individua automaticamente
le chiome di alberi e arbusti a partire da un'ortofoto, con calcolo di area,
perimetro, compattezza e — se disponibili un DSM e un DTM — altezza di ogni
chioma e separazione delle chiome adiacenti.

**Stato: sperimentale.** Funziona bene su dati di buona qualità e risoluzione,
ma richiede occhio critico sui risultati — vedi [Limiti](#limiti-e-cosa-sapere-prima-di-usarlo)
prima di usarlo per qualsiasi decisione che conta.

## Cosa fa

1. Calcola un indice di vegetazione (NDVI se l'ortofoto ha una banda NIR,
   altrimenti ExG dalle sole bande RGB)
2. Sogliatura automatica (Otsu) o manuale
3. Se fornisci DSM e DTM: calcola il CHM (Canopy Height Model) e lo usa per
   escludere vegetazione bassa (prati, colture) prima ancora di generare i
   poligoni
4. Separa le chiome adiacenti/fuse con un algoritmo di tipo watershed, invece
   di trattarle come un unico blob da accettare o scartare in massa
5. Filtra il rumore residuo per area e compattezza
6. Restituisce un layer vettoriale con area, perimetro, compattezza, altezza
   media/massima (se disponibile), e una classificazione dimensionale
   automatica incrociata con l'altezza

## Perché

Pensato per chi ha bisogno di una stima rapida della copertura arborea di
un'area — pianificazione, analisi territoriale preliminare, studi di
architettura del paesaggio — partendo da ortofoto e modelli altimetrici
pubblici (geoportali regionali/nazionali), senza dover digitalizzare le
chiome a mano una per una.

**Non è, e non vuole essere, un sostituto di rilievi certificati o censimenti
ufficiali.** È uno strumento di analisi preliminare/esplorativa: utile per
farsi un'idea rapida prima di rivolgersi a un professionista o all'ente
competente quando serve un dato con valore legale/amministrativo (permessi,
contenziosi, atti pubblici, pratiche catastali).

## Installazione

1. Scarica `canopy_mapper_algorithm.py`
2. In QGIS: **Processing → Toolbox** → icona Python (in alto) → **Add Script
   to Toolbox...** → seleziona il file
3. Lo strumento compare sotto **Scripts → CanopyMapper → Individua canopy da
   ortofoto**

In alternativa, copia il file nella cartella
`processing/scripts/` del tuo profilo QGIS (**Impostazioni → Profili utente
→ Apri cartella profilo attivo**) per caricarlo automaticamente ad ogni avvio.

### Requisiti

- QGIS 3.x
- `numpy` (incluso di serie in ogni installazione QGIS)
- `scipy` (incluso di serie in QGIS 3.4x; su versioni più vecchie potrebbe
  mancare — se assente, la sola separazione watershed viene saltata
  automaticamente con un avviso, il resto dello strumento funziona comunque)

## Limiti e cosa sapere prima di usarlo

Emersi da test estensivi su dati reali, non da considerazioni teoriche:

- **La risoluzione dell'ortofoto conta moltissimo.** Sotto i 20-30 cm/pixel
  la segmentazione delle singole chiome diventa inaffidabile, soprattutto per
  arbusti e alberi piccoli/giovani
- **Un indice di vegetazione da solo non distingue un albero da una coltura
  agricola.** Un campo di grano o un vigneto in pieno vigore possono avere un
  indice di vegetazione alto quanto una chioma vera. Fornire DSM e DTM aiuta
  molto (l'altezza separa bene alberi da colture basse), ma non risolve tutto
  — vigneti/frutteti con strutture (pali, tralicci) superano comunque la
  soglia di altezza minima
- **DSM e DTM devono essere temporalmente coerenti tra loro e con
  l'ortofoto.** Se provengono da rilievi fatti in periodi diversi, un campo
  che oggi appare arato/nudo nell'ortofoto può risultare "alto" nel CHM se al
  momento del rilievo LiDAR c'era ancora una coltura in piedi — un problema
  di dati, non risolvibile via software
- **La separazione automatica delle chiome (watershed) è un compromesso.**
  Riduce ma non elimina il rischio di frammentare chiome molto vicine o,
  viceversa, di non separare bene chiome molto ravvicinate — il parametro
  "distanza minima tra cime" va tarato sul contesto (specie, densità,
  risoluzione del dato)
- **La classe `sospetto_area_ampia_bassa`** è un campanello d'allarme
  automatico per poligoni con area estesa ma altezza media bassa — un buon
  punto di partenza per un controllo visivo manuale, non una garanzia di
  correttezza per tutto il resto

In sintesi: **verifica sempre i risultati su un campione noto** (un'area di
cui conosci il numero reale di alberi) prima di fidarti dei numeri su
un'area nuova.

## Output

Layer vettoriale poligonale con i campi:

| Campo | Significato |
|---|---|
| `id_pianta` | numerazione progressiva |
| `area` | area della chioma in m² |
| `perimeter` | perimetro in m |
| `compattez` | indice di compattezza (0-1) |
| `h_mean`, `h_max` | altezza media/massima (solo con DSM+DTM) |
| `classe_dim` | arbusto / albero_piccolo / albero_medio / albero_grande / sospetto_area_ampia_bassa |
| `specie` | campo vuoto, per compilazione manuale |

## Licenza

GPL v2 o successiva — vedi [LICENSE](LICENSE). Coerente con la licenza di
QGIS stesso, come richiesto per qualsiasi strumento distribuito nel suo
ecosistema.

## Contribuire

Segnalazioni di bug, casi di test su altri dataset, e miglioramenti sono
benvenuti — apri una issue o una pull request. Se segnali un problema, un
estratto piccolo di ortofoto/DSM/DTM su cui riproducibile aiuta moltissimo
a diagnosticarlo.

## Storia del progetto

Nato da un video su un workflow simile per QGIS, sviluppato e testato
iterativamente su dati reali (ortofoto + DSM + DTM LiDAR di un'area in
Emilia-Romagna) fino a raggiungere un comportamento affidabile sui casi
comuni. Vedi [CHANGELOG.md](CHANGELOG.md) per la cronologia delle versioni.
