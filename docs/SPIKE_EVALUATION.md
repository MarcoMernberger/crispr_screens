# Spike-In Evaluation für CRISPR Screen Qualitätskontrolle

Diese Module evaluieren die Performance verschiedener MAGeCK-Normalisierungs- und Analysevarianten anhand von Spike-In-Kontrollen.

## Übersicht

Die `evaluate_spike_in_performance_job` Funktion berechnet umfassende Qualitätsmetriken für CRISPR-Screen-Analysen mit Spike-In-Kontrollen. Die Spike-Ins werden in drei Kategorien eingefügt:

- **SPIKE_POS_G001..G020**: Positive Kontrollen (enriched)
- **SPIKE_NEG_G001..G020**: Negative Kontrollen (depleted)
- **SPIKE_NEUTRAL_G001..G020**: Neutrale Kontrollen (unverändert)

## Berechnete Metriken

### 1. Detection Performance

#### Precision, Recall, F1
- **Precision**: Anteil der korrekt identifizierten Hits unter allen als signifikant markierten Genen
- **Recall**: Anteil der detektierten Spike-Ins unter allen erwarteten Spike-Ins
- **F1 Score**: Harmonisches Mittel von Precision und Recall (wichtigste Metrik)

Formeln:
```
Precision = TP / (TP + FP)
Recall = TP / (TP + FN)
F1 = 2 × (Precision × Recall) / (Precision + Recall)
```

#### AUC-ROC und AUC-PR
- **AUC-ROC**: Area Under the ROC Curve - misst die Trennbarkeit zwischen True Positives und False Positives
- **AUC-PR**: Area Under the Precision-Recall Curve - besonders nützlich bei unbalancierten Datensätzen

Beide Werte liegen zwischen 0 und 1, wobei 1 = perfekt.

### 2. Separation Metrics

#### Cohen's d
Effektstärke, die die Trennung zwischen Spike-In-Gruppen misst:

```
Cohen's d = (μ₁ - μ₂) / σ_pooled
```

- |d| < 0.2: kleiner Effekt
- 0.2 ≤ |d| < 0.8: mittlerer Effekt
- |d| ≥ 0.8: großer Effekt

Wichtige Vergleiche:
- `neg_vs_neutral_cohens_d`: Trennung zwischen NEG und NEUTRAL Spike-Ins
- `pos_vs_neutral_cohens_d`: Trennung zwischen POS und NEUTRAL Spike-Ins
- `pos_vs_neg_cohens_d`: Trennung zwischen POS und NEG Spike-Ins

#### Mean/Median Differences
Absolute Unterschiede in den Log Fold Changes zwischen Spike-In-Gruppen.

#### Mann-Whitney U p-values
Statistische Signifikanz der Trennung zwischen Gruppen.

### 3. Ranking Power

#### AUCC (Area Under Cumulative Curve)
Misst, wie früh die wahren Hits im Ranking erscheinen:

```
AUCC = ∫ (cumulative fraction of hits) d(normalized rank)
```

- AUCC = 1.0: alle Hits stehen ganz oben im Ranking (perfekt)
- AUCC = 0.5: Hits sind zufällig verteilt
- AUCC < 0.5: Hits stehen eher unten im Ranking (schlecht)

#### Top-k Enrichment
Anteil der erwarteten Hits in den Top-k Genen:
- `frac_in_top_10`: Fraktion der Hits in Top 10
- `frac_in_top_50`: Fraktion der Hits in Top 50
- `frac_in_top_100`: Fraktion der Hits in Top 100
- `frac_in_top_500`: Fraktion der Hits in Top 500

#### Median/Mean Rank
Mittlerer Rang der erwarteten Hits. **Niedriger ist besser.**

### 4. Consistency

#### Coefficient of Variation (CV)
Relative Streuung innerhalb einer Spike-In-Gruppe:

```
CV = σ / |μ|
```

**Niedriger CV = höhere Konsistenz = besser**

Metriken:
- `neg_cv`: CV der NEG Spike-Ins
- `pos_cv`: CV der POS Spike-Ins
- `neutral_cv`: CV der NEUTRAL Spike-Ins

#### Interquartile Range (IQR)
Robustes Streuungsmaß:

```
IQR = Q₃ - Q₁
```

Weniger sensitiv gegenüber Ausreißern als Standardabweichung.

## Composite Score

Der `composite_score` kombiniert alle Metriken zu einer Gesamtbewertung:

```python
composite_score = Σ (weight_i × normalized_metric_i) / Σ weight_i
```

Standard-Gewichte:
- F1 Score: 3.0 (wichtigste Metrik)
- AUC-ROC: 2.0
- AUC-PR: 2.0
- AUCC: 2.0
- Cohen's d (separation): 1.5
- CV (consistency): -1.0 (niedriger ist besser)
- Median Rank: -1.0 (niedriger ist besser)

## Verwendung

```python
from crispr_screens import evaluate_spike_in_performance_job

# Dictionary mit MAGeCK-Ergebnissen
spike_test_methods = {
    "Paired_Median": Path("results/.../gene_summary.tsv"),
    "Unpaired_Median": Path("results/.../gene_summary.tsv"),
    "Paired_Control": Path("results/.../gene_summary.tsv"),
    "Unpaired_Control": Path("results/.../gene_summary.tsv"),
}

# Evaluations-Job erstellen
eval_job = evaluate_spike_in_performance_job(
    output_file=Path("results/spike_evaluation.tsv"),
    mageck_results=spike_test_methods,
    direction="neg",  # "neg" oder "pos"
    fdr_threshold=0.05,
    lfc_threshold=1.0,
    weights={
        "f1": 3.0,
        "auc_roc": 2.0,
        "aucc": 2.0,
        "neg_vs_neutral_cohens_d": 1.5,
        "neg_cv": -1.0,
    },
    dependencies=[rra1, rra2, rra3, rra4],
)
```

## Output

Die Evaluierung erzeugt:

1. **TSV-Datei** mit allen Metriken für jede MAGeCK-Variante
2. **Console Output** mit einer Zusammenfassung der Top-3-Methoden:

```
================================================================================
SPIKE-IN EVALUATION SUMMARY
================================================================================

Direction: neg selection
FDR threshold: 0.05
LFC threshold: 1.0

Number of methods evaluated: 4

--------------------------------------------------------------------------------
TOP 3 METHODS (by composite score):
--------------------------------------------------------------------------------

1. Paired_Median
   Composite Score: 0.847
   F1: 0.923 | Precision: 0.950 | Recall: 0.900
   AUC-ROC: 0.982 | AUCC: 0.891
   Median Rank: 23.5
   Detected: 18/20 expected hits
   Separation (Cohen's d): 3.45

2. Paired_Control
   Composite Score: 0.812
   ...
```

## Interpretation

**Beste Methode:** Höchster `composite_score` und `rank = 1`

**Wichtige Kriterien:**
- **F1 > 0.8**: Gute Balance zwischen Precision und Recall
- **AUC-ROC > 0.9**: Sehr gute Diskriminierung
- **AUCC > 0.7**: Gutes frühes Ranking
- **Cohen's d > 2.0**: Sehr gute Separation
- **CV < 0.3**: Gute Konsistenz

## Troubleshooting

### NaN-Werte in Metriken
- Zu wenige Spike-Ins detektiert
- Fehlende Spalten in MAGeCK-Output
- Alle Werte gleich (keine Varianz)

### Niedrige Scores
- **Niedriger F1**: Schlechte Detection → evtl. zu strenge Thresholds
- **Niedriger AUCC**: Schlechtes Ranking → Normalisierung überprüfen
- **Niedriger Cohen's d**: Schlechte Separation → mehr Replikate oder anderes log_effect
- **Hoher CV**: Inkonsistenz → technische Variabilität reduzieren

## Weiterführende Analysen

Nach der Evaluation können Sie:

1. **Volcano Plots** für die beste Methode erstellen
2. **Scatter Plots** zum Vergleich von zwei Methoden
3. **Heatmaps** der Separation-Metriken
4. **Box Plots** der LFCs pro Spike-In-Gruppe

## Literatur

- Precision/Recall: Standard ML metrics
- Cohen's d: Cohen, J. (1988). Statistical Power Analysis
- AUC-ROC: Fawcett, T. (2006). ROC analysis
- AUCC: Inspired by enrichment score metrics
