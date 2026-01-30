# MAGeCK Method Comparison - Quick Start

## Was ist das?

Ein Tool zur systematischen Bewertung, welche MAGeCK-Methode (RRA vs MLE) für deinen CRISPR-Screen am besten geeignet ist.

## Wann brauchst du das?

- ✓ Du hast RRA und MLE laufen lassen und willst wissen, welche robuster ist
- ✓ Du hast Batch-Effekte und willst MLE mit/ohne Batch-Korrektur vergleichen
- ✓ Du willst verschiedene Normalisierungsmethoden vergleichen (median, control, total)
- ✓ Du willst Paired vs. Unpaired Design testen
- ✓ Du brauchst quantitative Metriken für dein Paper

## Die 4 Analysen

### 1. **Replikat-Konsistenz** (Leave-one-replicate-out)
→ Ist die Methode stabil wenn ein Replikat fehlt?

### 2. **sgRNA-Kohärenz** 
→ Zeigen mehrere sgRNAs pro Gen in die gleiche Richtung?

### 3. **Falsch-Positiv-Rate** (Control sgRNAs)
→ Wie viele Non-Targeting Controls sind in den Top-Hits?

### 4. **Permutations-Test**
→ Findet die Methode Hits in zufälligen Daten?

## Installation

```bash
# Bereits installiert im crispr_screens package!
# Keine zusätzliche Installation nötig
```

## Minimal Example (5 Zeilen)

```python
import pypipegraph2 as ppg
from crispr_screens import mageck_method_comparison_job
from crispr_screens.core.mageck import mageck_test, mageck_mle

ppg.new()

# Definiere Methoden
methods = {
    "RRA_paired": {
        "run_func": mageck_test,
        "params": {"paired": True, "norm_method": "median"},
        "gene_col": "id",
    },
    "MLE": {
        "run_func": mageck_mle,
        "params": {"design_matrix": "design.tsv"},
        "gene_col": "Gene",
    },
}

# Run comparison
job = mageck_method_comparison_job(
    count_table="results/mageck_count/counts.txt",
    control_ids=["Total_Rep1", "Total_Rep2", "Total_Rep3"],
    treatment_ids=["Sort_Rep1", "Sort_Rep2", "Sort_Rep3"],
    output_dir="results/method_comparison",
    control_sgrnas="cache/input/controls.txt",
    methods=methods,
)

ppg.run()
```

## Results ansehen

```python
import pandas as pd

# Haupt-Summary
summary = pd.read_csv(
    "results/method_comparison/method_comparison_summary.tsv",
    sep="\t"
)

print(summary)
```

**Beispiel Output:**

```
method        mean_spearman  mean_jaccard_top_100  controls_in_top_100
RRA_paired    0.78           0.62                  5
MLE           0.91           0.81                  1
```

**Interpretation:**
- MLE hat höhere Spearman Korrelation (0.91 vs 0.78) → stabilere Rankings
- MLE hat höheren Jaccard Index (0.81 vs 0.62) → konsistentere Top-Hits
- MLE hat weniger Controls in Top 100 (1 vs 5) → weniger False-Positives

**→ MLE ist besser für diesen Datensatz!**

## Komplettes Beispiel

Siehe: `examples/method_comparison_example.py`

## Ausführliche Dokumentation

Siehe: `docs/method_comparison_guide.md`

- Detaillierte Erklärung jeder Metrik
- Interpretation Guidelines
- Troubleshooting
- Best Practices

## Typische Use Cases

### Use Case 1: RRA vs MLE bei Batch-Effekten

```python
methods = {
    "RRA_paired": {
        "run_func": mageck_test,
        "params": {"paired": True, "norm_method": "median"},
        "gene_col": "id",
    },
    "MLE_with_batch": {
        "run_func": mageck_mle,
        "params": {
            "design_matrix": "design_with_batch.tsv",
            "norm_method": "median"
        },
        "gene_col": "Gene",
    },
}
```

### Use Case 2: Normalisierung vergleichen

```python
methods = {
    "RRA_median": {
        "run_func": mageck_test,
        "params": {"paired": True, "norm_method": "median"},
        "gene_col": "id",
    },
    "RRA_control": {
        "run_func": mageck_test,
        "params": {
            "paired": True,
            "norm_method": "control",
            "control_sgrnas": "controls.txt"
        },
        "gene_col": "id",
    },
    "RRA_total": {
        "run_func": mageck_test,
        "params": {"paired": True, "norm_method": "total"},
        "gene_col": "id",
    },
}
```

### Use Case 3: Paired vs Unpaired

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

## Output Struktur

```
results/method_comparison/
├── method_comparison_summary.tsv           # ← START HIER!
├── RRA_paired/
│   ├── leave_one_out/
│   │   └── RRA_paired_leave_one_out_comparison.tsv
│   ├── RRA_paired_sgrna_coherence.tsv
│   ├── RRA_paired_control_false_positives.tsv
│   └── permutations/
└── MLE/
    └── ... (gleiche Struktur)
```

## Interpretation Cheat Sheet

| Metrik                | Gut   | Akzeptabel | Schlecht |
| --------------------- | ----- | ---------- | -------- |
| Spearman Corr         | > 0.8 | 0.6 - 0.8  | < 0.6    |
| Jaccard Top-100       | > 0.7 | 0.4 - 0.7  | < 0.4    |
| Direction Consistency | > 0.8 | 0.6 - 0.8  | < 0.6    |
| Controls in Top-100   | 0-2   | 3-5        | > 5      |
| Perm Hits (% of orig) | < 10% | 10-25%     | > 25%    |

**Höhere Werte = Besser** (außer Controls und Perm Hits, da niedriger = besser)

## Einzelne Analysen laufen lassen

Falls du nicht die komplette Comparison brauchst:

```python
from crispr_screens import (
    leave_one_replicate_out_job,      # Nur Replicate Consistency
    sgrna_coherence_job,               # Nur sgRNA Coherence
    control_false_positive_job,        # Nur Control FP Check
    permutation_test_job,              # Nur Permutation Test
)

# Beispiel: Nur Control False-Positives checken
job = control_false_positive_job(
    gene_summary="results/rra/gene_summary.tsv",
    control_sgrnas="cache/input/controls.txt",
    output_file="results/comparison/control_fp.tsv",
)
```

## Tests

```bash
# Prüfe dass alles funktioniert
cd code/crispr_screens
python tests/test_method_comparison.py
```

## Häufige Fragen

**Q: Wie lange dauert die Analyse?**
A: Für 3 Replikate ca. 5-30 Minuten pro Methode (abhängig von MAGeCK Runtime)

**Q: Brauche ich mindestens 3 Replikate?**
A: Ja, für Leave-one-out. Mit 2 Replikaten kannst du nur Coherence/FP/Perm laufen lassen.

**Q: Kann ich mehr als 2 Methoden vergleichen?**
A: Ja! Füge einfach mehr Einträge zum `methods` Dictionary hinzu.

**Q: Was wenn eine Methode klar besser ist?**
A: Nimm diese für deine finale Analyse und erwähne die Comparison im Methods-Teil deines Papers.

**Q: Was wenn beide ähnlich sind?**
A: Wähle RRA wenn simpler Design, MLE wenn Batch-Effekte. Oder nimm Intersection der Top-Hits beider.

## Support

- Beispiele: `examples/method_comparison_example.py`
- Dokumentation: `docs/method_comparison_guide.md`
- Tests: `tests/test_method_comparison.py`

## Zitation

Wenn du dieses Tool verwendest, zitiere bitte:

```
Method comparison performed using crispr-screens package (v X.Y.Z).
MAGeCK: Li et al. (2014) Genome Biology 15:554.
```

## Nächste Schritte

1. ✓ Lade das komplette Beispiel: `examples/method_comparison_example.py`
2. ✓ Passe Pfade und Sample-IDs an
3. ✓ Run die Analyse
4. ✓ Check `method_comparison_summary.tsv`
5. ✓ Wähle beste Methode
6. ✓ Führe finale Analyse mit gewählter Methode durch
