# MAGeCK Method Comparison - Comprehensive Guide

## Übersicht

Dieses Modul ermöglicht einen systematischen Vergleich zwischen verschiedenen MAGeCK-Analysemethoden (RRA vs. MLE) um die robusteste und zuverlässigste Methode für deinen Datensatz zu identifizieren.

## Die vier Säulen der Methodenvergleichsanalyse

### 1. Leave-One-Replicate-Out Analyse (Replikat-Konsistenz)

**Was wird gemessen:**
- Stabilität der Ergebnisse über Replikate hinweg
- Reproduzierbarkeit der Top-Hits
- Robustheit gegenüber einzelnen Replikaten

**Wie es funktioniert:**
```python
# Beispiel mit 3 Replikaten:
# Run 1: Rep1 + Rep2 (ohne Rep3)
# Run 2: Rep1 + Rep3 (ohne Rep2)  
# Run 3: Rep2 + Rep3 (ohne Rep1)
# 
# Vergleiche dann die Ergebnisse paarweise
```

**Metriken:**
- **Spearman Correlation**: Rang-Korrelation der Gen-Scores zwischen Runs
  - > 0.8: Sehr stabil
  - 0.6 - 0.8: Moderat stabil
  - < 0.6: Instabil
  
- **Jaccard Index (Top-N)**: Überlapp der Top-N Gene zwischen Runs
  - > 0.7: Hoher Überlapp
  - 0.4 - 0.7: Moderater Überlapp
  - < 0.4: Geringer Überlapp

**Interpretation:**
```
Methode A:
  mean_spearman = 0.85
  mean_jaccard_top_100 = 0.72

Methode B:
  mean_spearman = 0.65
  mean_jaccard_top_100 = 0.45

→ Methode A ist deutlich stabiler und reproduzierbarer
```

---

### 2. sgRNA Kohärenz-Analyse

**Was wird gemessen:**
- Konsistenz mehrerer sgRNAs pro Gen
- Richtungsübereinstimmung (alle positiv oder alle negativ)
- Stärke des Effekts

**Wie es funktioniert:**
Für jedes Top-Gen:
1. Sammle alle sgRNAs die dieses Gen targetieren
2. Prüfe: Zeigen die sgRNAs in die gleiche Richtung?
3. Prüfe: Ist die Effektstärke ähnlich?
4. Zähle: Wie viele sgRNAs sind signifikant?

**Metriken:**
- **Direction Consistency**: Anteil der sgRNAs mit gleichem Vorzeichen
  - > 0.8: Sehr kohärent
  - 0.6 - 0.8: Moderat kohärent
  - < 0.6: Inkohärent (möglicher Jackpot-Effekt)

- **Fraction Significant**: Anteil signifikanter sgRNAs
  - Höher = besser (nicht nur eine sgRNA verantwortlich)

- **CV (Coefficient of Variation)**: log-fold-change Variabilität
  - Niedriger = konsistenter

**Interpretation:**
```
Gen BRCA1:
  n_sgrnas = 6
  direction_consistency = 1.0  (alle 6 negativ)
  fraction_significant = 0.83  (5 von 6 signifikant)
  cv_lfc = 0.3                 (ähnliche Effektstärke)
  
→ Sehr robustes Signal, nicht nur ein Jackpot-Guide

Gen XYZ:
  n_sgrnas = 6
  direction_consistency = 0.67  (4 negativ, 2 positiv)
  fraction_significant = 0.17   (nur 1 signifikant)
  cv_lfc = 2.5                  (sehr variable Effekte)
  
→ Wahrscheinlich Artefakt oder Jackpot-Guide
```

**Wann ist eine Methode besser?**
- Höhere durchschnittliche Direction Consistency bei Top-Genes
- Höhere durchschnittliche Fraction Significant
- Niedrigere CV bei wahren Hits

---

### 3. Control sgRNA False-Positive Analyse

**Was wird gemessen:**
- Wie viele Non-Targeting Controls erscheinen in den Top-N Hits?
- Falsch-Positiv-Rate der Methode

**Wie es funktioniert:**
```python
# Zähle Control-sgRNAs in:
# Top 50 Genes
# Top 100 Genes
# Top 200 Genes
# Top 500 Genes

# Erwartung: Sollten 0 oder sehr wenige sein
```

**Metriken:**
- **n_controls**: Absolute Anzahl Controls in Top-N
- **fraction_controls**: Anteil Controls an Top-N

**Interpretation:**
```
Methode A (Top 100):
  n_controls = 2
  fraction = 0.02
  
Methode B (Top 100):
  n_controls = 15
  fraction = 0.15
  
→ Methode A ist besser kalibriert
→ Methode B hat zu viele Falsch-Positive
```

**Typische Werte:**
- Gut kalibriert: 0-2% Controls in Top 100
- Akzeptabel: 2-5% Controls
- Schlecht: > 5% Controls

**Was tun bei hohen False-Positives?**
1. Andere Normalisierungsmethode versuchen
2. Strengere FDR-Cutoffs verwenden
3. sgRNA-Kohärenz als zusätzliches Filter

---

### 4. Permutations-Test

**Was wird gemessen:**
- Wie verhält sich die Methode bei kaputten/zufälligen Daten?
- Baseline False-Positive Rate

**Wie es funktioniert:**
```python
# Typ 1: Within-Sample Permutation
# → Mische Counts innerhalb jeder Sample-Spalte
# → Zerstört echte biologische Signale
# → Behalte technische Artefakte

# Typ 2: Label-Swap
# → Vertausche Gene-Labels
# → Komplette Randomisierung
```

**Erwartung:**
- Echte Signale sollten verschwinden oder schwächer werden
- Technische Artefakte bleiben bestehen
- Eine robuste Methode findet wenige "Hits" in permutierten Daten

**Metriken:**
- **mean_perm_sig_genes**: Durchschnittliche Anzahl signifikanter Gene in Permutationen
- **std_perm_sig_genes**: Variabilität über Permutationen

**Interpretation:**
```
Methode A:
  Original: 500 signifikante Gene
  Permutiert: mean = 10, std = 5
  → Gut! Nur 2% der Hits sind False Positives

Methode B:
  Original: 500 signifikante Gene
  Permutiert: mean = 150, std = 50
  → Schlecht! 30% könnten False Positives sein
```

---

## Comprehensive Comparison Workflow

### Schritt 1: Setup

```python
import pypipegraph2 as ppg
from pathlib import Path
from crispr_screens import mageck_method_comparison_job
from crispr_screens.core.mageck import mageck_test, mageck_mle

ppg.new()

# Definiere zu vergleichende Methoden
methods = {
    "RRA_paired_median": {
        "run_func": mageck_test,
        "params": {"paired": True, "norm_method": "median"},
        "gene_col": "id",
    },
    "MLE_batch_median": {
        "run_func": mageck_mle,
        "params": {
            "design_matrix": "design_batch.tsv",
            "norm_method": "median"
        },
        "gene_col": "Gene",
    },
}
```

### Schritt 2: Run Comparison

```python
comparison_job = mageck_method_comparison_job(
    count_table="results/counts.txt",
    control_ids=["Total_Rep1", "Total_Rep2", "Total_Rep3"],
    treatment_ids=["Sort_Rep1", "Sort_Rep2", "Sort_Rep3"],
    output_dir="results/method_comparison",
    control_sgrnas="controls.txt",
    methods=methods,
    top_n_list=[50, 100, 200],
    run_leave_one_out=True,
    run_coherence=True,
    run_control_fp=True,
    run_permutation=True,
    n_permutations=5,
)

ppg.run()
```

### Schritt 3: Ergebnisse interpretieren

```python
import pandas as pd

# Lade Summary
summary = pd.read_csv(
    "results/method_comparison/method_comparison_summary.tsv",
    sep="\t"
)

print(summary)
```

**Output Beispiel:**
```
method               mean_spearman  mean_jaccard_top_100  mean_direction_consistency  controls_in_top_100  mean_perm_sig_genes
RRA_paired_median    0.82           0.68                  0.85                        3                    12
MLE_batch_median     0.88           0.75                  0.91                        1                    8
```

**Interpretation:**
- MLE hat höhere Spearman Korrelation → stabilere Rankings
- MLE hat höheren Jaccard Index → mehr konsistente Top-Hits
- MLE hat höhere Direction Consistency → kohärentere sgRNAs
- MLE hat weniger Controls in Top 100 → weniger False Positives
- MLE hat weniger Hits in Permutationen → besser kalibriert

**→ MLE ist für diesen Datensatz die bessere Methode!**

---

## Detaillierte Analyse pro Methode

Für jede Methode gibt es ein eigenes Verzeichnis:

```
results/method_comparison/
├── method_comparison_summary.tsv        # Gesamtübersicht
├── RRA_paired_median/
│   ├── leave_one_out/
│   │   ├── RRA_paired_median_leave_one_out_comparison.tsv
│   │   ├── leave_out_rep1/
│   │   │   └── RRA_paired_median_leave_out_rep1.gene_summary.tsv
│   │   ├── leave_out_rep2/
│   │   └── leave_out_rep3/
│   ├── full_run/
│   │   ├── RRA_paired_median_full.gene_summary.tsv
│   │   └── RRA_paired_median_full.sgrna_summary.txt
│   ├── RRA_paired_median_sgrna_coherence.tsv
│   ├── RRA_paired_median_control_false_positives.tsv
│   └── permutations/
│       ├── permutation_within_sample_1/
│       ├── permutation_within_sample_2/
│       └── ...
└── MLE_batch_median/
    └── ... (gleiche Struktur)
```

### Leave-One-Out Details

```python
loo = pd.read_csv(
    "results/method_comparison/MLE_batch_median/leave_one_out/"
    "MLE_batch_median_leave_one_out_comparison.tsv",
    sep="\t"
)

# Paarweise Vergleiche
print(loo[["run1", "run2", "spearman_correlation", "top_100_jaccard"]])
```

### sgRNA Coherence Details

```python
coherence = pd.read_csv(
    "results/method_comparison/MLE_batch_median/"
    "MLE_batch_median_sgrna_coherence.tsv",
    sep="\t"
)

# Top Genes mit bester Kohärenz
best_genes = coherence.nlargest(20, "direction_consistency")
print(best_genes[["gene", "n_sgrnas", "direction_consistency", "fraction_significant"]])
```

---

## Entscheidungsbaum

```
START: Welche Methode ist besser?
│
├─ Schritt 1: Leave-One-Out Analyse
│  │
│  ├─ Hohe Korrelation (> 0.8) & Jaccard (> 0.7)?
│  │  └─ JA → Methode ist stabil ✓
│  │  └─ NEIN → Methode instabil, andere Normalisierung testen
│  │
│  └─ Große Unterschiede zwischen Methoden?
│     └─ JA → Wähle stabilere Methode
│     └─ NEIN → Beide ähnlich, weiter mit Schritt 2
│
├─ Schritt 2: sgRNA Kohärenz
│  │
│  ├─ Hohe Direction Consistency (> 0.8)?
│  │  └─ JA → sgRNAs konsistent ✓
│  │  └─ NEIN → Viele Jackpot-Guides, Library-Problem?
│  │
│  └─ Unterschiede zwischen Methoden?
│     └─ JA → Wähle Methode mit höherer Konsistenz
│     └─ NEIN → Beide ähnlich, weiter mit Schritt 3
│
├─ Schritt 3: Control False-Positives
│  │
│  ├─ Wenige Controls in Top-N (< 5%)?
│  │  └─ JA → Gut kalibriert ✓
│  │  └─ NEIN → Zu liberal, strengere Cutoffs oder andere Norm
│  │
│  └─ Unterschiede zwischen Methoden?
│     └─ JA → Wähle Methode mit weniger FP
│     └─ NEIN → Beide ähnlich, weiter mit Schritt 4
│
└─ Schritt 4: Permutations-Test
   │
   ├─ Wenige Hits in Permutationen (< 5% von Original)?
   │  └─ JA → Gut kalibriert ✓
   │  └─ NEIN → Zu viele False Positives
   │
   └─ Unterschiede zwischen Methoden?
      └─ JA → Wähle Methode mit weniger Perm-Hits
      └─ NEIN → Beide ähnlich
      
ENDE: Wähle Methode mit besten Overall-Metriken
```

---

## Häufige Szenarien

### Szenario 1: MLE vs. RRA bei Batch-Effekten

**Symptome:**
- RRA zeigt niedrige Replicate Consistency
- Große Unterschiede zwischen Runs aus verschiedenen Batches

**Lösung:**
```python
# Teste MLE mit Batch-Design-Matrix
methods = {
    "RRA_paired": {
        "run_func": mageck_test,
        "params": {"paired": True, "norm_method": "median"},
        "gene_col": "id",
    },
    "MLE_with_batch": {
        "run_func": mageck_mle,
        "params": {
            "design_matrix": "design_with_batch.tsv",  # Enthält Batch-Spalte
            "norm_method": "median"
        },
        "gene_col": "Gene",
    },
}
```

**Erwartung:**
- MLE sollte höhere Spearman Korrelation zeigen
- MLE sollte mehr konsistente Top-Hits haben

---

### Szenario 2: Normalisierungsmethoden vergleichen

**Symptome:**
- Control sgRNAs zeigen starke Varianz
- Viele Controls in Top-Hits

**Lösung:**
```python
methods = {
    "RRA_median_norm": {
        "run_func": mageck_test,
        "params": {"paired": True, "norm_method": "median"},
        "gene_col": "id",
    },
    "RRA_control_norm": {
        "run_func": mageck_test,
        "params": {
            "paired": True,
            "norm_method": "control",
            "control_sgrnas": "controls.txt"
        },
        "gene_col": "id",
    },
    "RRA_total_norm": {
        "run_func": mageck_test,
        "params": {"paired": True, "norm_method": "total"},
        "gene_col": "id",
    },
}
```

**Prüfe:**
- Control False-Positives pro Methode
- Welche Norm hat wenigste Controls in Top-N?

---

### Szenario 3: Paired vs. Unpaired Design

**Symptome:**
- Biologische Replikate sind matched (z.B. gleiche Zelllinie)
- Unsicher ob Pairing hilft

**Lösung:**
```python
methods = {
    "RRA_paired": {
        "run_func": mageck_test,
        "params": {"paired": True, "norm_method": "median"},
        "gene_col": "id",
    },
    "RRA_unpaired": {
        "run_func": mageck_test,
        "params": {"paired": False, "norm_method": "median"},
        "gene_col": "id",
    },
}
```

**Erwartung:**
- Wenn Replikate wirklich matched: Paired sollte besser sein
- Höhere Replicate Consistency
- Mehr Power (mehr signifikante Hits)

---

## Best Practices

### 1. Anzahl Permutationen

```python
# Für schnelle Tests: 3-5 Permutationen
n_permutations = 3

# Für robuste Schätzung: 10-20 Permutationen
n_permutations = 10

# Für Publication: 50-100 Permutationen
n_permutations = 50
```

### 2. Top-N Listen

```python
# Standard
top_n_list = [50, 100, 200]

# Großer Screen
top_n_list = [100, 250, 500, 1000]

# Kleiner Screen
top_n_list = [25, 50, 100]
```

### 3. Minimal Replicates

```python
# Minimum für Leave-One-Out: 3 Replikate
# Mit 2 Replikaten: Nur 1 Leave-Out möglich, keine Vergleiche

# Empfohlen: 3 Replikate
control_ids = ["Total_Rep1", "Total_Rep2", "Total_Rep3"]

# Gut: 4+ Replikate
control_ids = ["Total_Rep1", "Total_Rep2", "Total_Rep3", "Total_Rep4"]
```

---

## Troubleshooting

### Problem: Leave-One-Out schlägt fehl

**Mögliche Ursachen:**
1. Zu wenig Replikate (< 2)
2. MAGeCK Run schlägt fehl
3. Design Matrix fehlt (für MLE)

**Lösung:**
```python
# Prüfe Logs in:
# results/method_comparison/<method>/leave_one_out/leave_out_rep1/

# Teste MAGeCK manuell:
from crispr_screens.core.mageck import mageck_test

mageck_test(
    count_table="counts.txt",
    control_ids=["Total_Rep1", "Total_Rep2"],  # Rep3 weggelassen
    treatment_ids=["Sort_Rep1", "Sort_Rep2"],
    out_dir="test_run",
    prefix="test",
    paired=True,
    norm_method="median"
)
```

### Problem: Keine Gene in sgRNA Coherence

**Mögliche Ursachen:**
1. Gene- und sgRNA-Summary passen nicht zusammen
2. Falsche Column Names
3. Keine Top-Genes gefunden

**Lösung:**
```python
# Prüfe Spalten-Namen
import pandas as pd
gene_df = pd.read_csv("gene_summary.tsv", sep="\t")
print(gene_df.columns)

sgrna_df = pd.read_csv("sgrna_summary.txt", sep="\t")
print(sgrna_df.columns)

# Passe gene_col an:
sgrna_coherence_job(
    ...,
    gene_col="Gene",  # oder "id"
    sgrna_gene_col="Gene",
)
```

### Problem: Permutations-Test findet zu viele Hits

**Das ist normal wenn:**
- Du viele echte Effekte hast
- Einige technische Artefakte bestehen bleiben

**Prüfe:**
```python
# Vergleiche mit Original
original_hits = 500
perm_mean_hits = 150

fraction = perm_mean_hits / original_hits
print(f"Permutation finds {fraction:.1%} of original hits")

# < 10% ist gut
# 10-25% ist akzeptabel
# > 25% ist problematisch
```

---

## API Reference

Siehe auch: `examples/method_comparison_example.py`

### Hauptfunktionen

```python
from crispr_screens import (
    mageck_method_comparison_job,      # Comprehensive comparison
    leave_one_replicate_out_job,       # Only replicate consistency
    sgrna_coherence_job,                # Only sgRNA coherence
    control_false_positive_job,         # Only control FP check
    permutation_test_job,               # Only permutation test
)
```

### Core Functions (ohne PyPipeGraph)

```python
from crispr_screens.core.method_comparison import (
    leave_one_replicate_out_analysis,
    analyze_sgrna_coherence,
    analyze_control_false_positives,
    permutation_test_analysis,
    compare_mageck_methods,
)
```

---

## Zitation

Wenn du diese Methode in deiner Publikation verwendest, zitiere bitte:

```
Diese Analyse verwendet MAGeCK Method Comparison Tools aus dem
crispr-screens Package (Version X.Y.Z) für die systematische
Bewertung von MAGeCK RRA und MLE Methoden.

MAGeCK Original Paper:
Li et al. (2014) MAGeCK enables robust identification of essential 
genes from genome-scale CRISPR/Cas9 knockout screens. 
Genome Biology 15:554.

Li et al. (2015) MAGeCK-VISPR: an integrated analysis and 
visualization platform for CRISPR screens. Genome Biology 16:281.
```

---

## Zusammenfassung

**Use Cases:**
1. Du hast Batch-Effekte → Vergleiche RRA vs. MLE
2. Du bist unsicher über Normalisierung → Vergleiche median/control/total
3. Du willst Paired vs. Unpaired testen
4. Du willst die robusteste Methode für dein Paper

**Output:**
- Eine klare Empfehlung welche Methode besser ist
- Quantitative Metriken für Reproduzierbarkeit
- Qualitätschecks für Library und Design

**Interpretation:**
- Höhere Metriken = bessere Methode
- Vergleiche relativ zwischen Methoden
- Wähle Methode mit besten Overall-Performance

**Next Steps:**
- Führe finale Analyse mit gewählter Methode durch
- Verwende stringentere Filter wenn nötig
- Validiere Top-Hits experimentell
