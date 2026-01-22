# Control sgRNA Qualitätskontrolle - Deutsche Kurzanleitung

## 🎯 Was wurde implementiert?

Vollständige Qualitätskontrolle für Control sgRNAs in CRISPR Screens mit mehreren Bedingungen.
Prüft ob die Control sgRNAs stabil und für Normalisierung geeignet sind.

## 🚀 Schnellstart

### 1. Anysnake2 Umgebung aktualisieren

```bash
cd /clara/ffs/e/20260105_AG_Stiewe_Andrea_Nist_sgRNA_Brunello_screen
anysnake2 shell
```

Dies installiert automatisch:
- crispr_screens Paket (editable mode)
- scipy, seaborn, scikit-learn, pandas

### 2. Quickstart Script ausführen

```bash
# Im anysnake2 shell:
cd code/crispr_screens
python quickstart_control_qc.py
```

Das Script wird:
1. Deine Daten überprüfen (Pfade, Spalten, Controls)
2. QC Analyse durchführen
3. Plots und Metriken generieren

### 3. Oder programmatisch verwenden

```python
from crispr_screens.core.qc import generate_control_qc_report

report = generate_control_qc_report(
    count_table="results/mageck_count/all/counts.count.tsv",
    control_sgrnas="incoming/control_sgRNAs.txt",
    baseline_condition="Total",  # Deine Baseline Bedingung
    output_dir="results/qc/control_sgrnas",
)
```

## 📊 Welche QC Checks werden durchgeführt?

### 1. **Verteilung pro Bedingung** (Univariat)
- Histogram der Δc = log2(CPM_Bedingung / CPM_Baseline)
- Median, IQR, Tail-Rate
- Wilcoxon Test (H0: Median = 0)

### 2. **Paarweise Vergleiche** (Matrix)
- Heatmap: Median-Shift zwischen allen Bedingungspaaren
- Sollte überall ~0 sein
- Zeigt globale Drifts

### 3. **Replikat-Konsistenz**
- Korrelationsmatrizen zwischen Replikaten
- Pro Bedingung separat
- Sollte hoch sein (>0.9)

### 4. **PCA Analyse** (Multivariat)
- PCA nur auf Controls
- 2 Versionen: nach Bedingung und nach Replikat gefärbt
- **Wichtig**: Controls sollten NICHT nach Bedingung trennen!

### 5. **Statistische Tests**
- Wilcoxon signed-rank Test pro Bedingung
- Effektgrößen (nicht nur p-Werte)

## 📈 Interpretation

### ✅ GUTE Controls (für Normalisierung geeignet)
- Median Δc ≈ 0 (± 0.1)
- Tail Rate < 5% (|Δc| > 1)
- Replikat-Korrelation > 0.9
- PCA: Proben clustern nach Replikat, NICHT nach Bedingung

### ⚠️ AKZEPTABLE Controls
- Median Δc: 0.1 - 0.3
- Tail Rate: 5-10%
- Korrelation: 0.8-0.9
- Mit Vorsicht verwendbar

### ❌ PROBLEMATISCHE Controls
- Median Δc > 0.5 → Systematischer Shift
- Tail Rate > 15% → Instabilität
- Korrelation < 0.8 → Batch Effekte
- PCA: Bedingungen trennen auf PC1 → Controls zeigen Behandlungseffekt

**→ In diesem Fall: Median- statt Control-Normalisierung verwenden!**

## 📁 Output Dateien

Nach dem Run findest du in `results/qc/control_sgrnas/`:

**Metriken:**
- `control_qc_metrics.tsv` - Alle Metriken pro Bedingung
- `control_qc_pairwise_shifts.tsv` - Paarweise Shifts Matrix

**Plots (PNG + PDF):**
- `control_qc_distribution` - Histogramme pro Bedingung
- `control_qc_pairwise_heatmap` - Paarweise Shifts Heatmap
- `control_qc_replicate_correlation` - Replikat-Korrelationen
- `control_qc_pca_condition` - PCA nach Bedingung gefärbt
- `control_qc_pca_replicate` - PCA nach Replikat gefärbt

## 🔧 Integration in dein run.py

```python
from crispr_screens import control_qc_job

# Nach deinem mageck_count_job:
qc_job = control_qc_job(
    count_table=mageck_count_dir / "counts.count.tsv",
    control_sgrnas=Path("incoming/control_sgRNAs.txt"),
    baseline_condition="Total",  # An deine Spalten anpassen
    output_dir=results_dir / "qc" / "control_sgrnas",
    prefix="control_qc",
    dependencies=[count_job],  # Hängt von count job ab
)
```

## 💡 Wichtige Hinweise

1. **Baseline Bedingung identifizieren:**
   - Muss mit Spalten-Präfix übereinstimmen
   - Typisch: "Total", "T0", "unsorted"
   - Bei dir wahrscheinlich: **"Total"** (für Total_Rep1, Total_Rep2, Total_Rep3)

2. **Spalten-Format:**
   - Erwartet: `Bedingung_Replikat` (z.B. "Sort1_Rep2")
   - Delimiter standardmäßig: `_`
   - Automatisches Gruppieren von Replikaten

3. **CPM Berechnung:**
   - Verwendet **raw counts**, nicht MAGeCK-normalisierte Werte
   - Vermeidet zirkuläre Logik
   - Manuell: CPM = (count / library_size) × 10^6

4. **Priorisierung der Plots:**
   - **Zuerst**: PCA Plots (zeigen globale Muster)
   - **Dann**: Verteilungsplots (Details pro Bedingung)
   - **Zuletzt**: Paarweise Heatmap (Übersicht)

## 🐛 Troubleshooting

### "No module named 'crispr_screens'"
→ `anysnake2 shell` ausführen (baut Umgebung neu)

### "Baseline condition not found"
→ Spaltennamen prüfen: `pd.read_csv(..., sep='\t').columns`
→ Baseline-Namen anpassen

### "No control sgRNAs found"
→ sgRNA IDs zwischen Count-Tabelle und Control-Datei vergleichen
→ Spaltenname korrekt? (default: "sgRNA")

### Speicherprobleme
→ Bei >100K sgRNAs: Count-Tabelle vorher auf Controls filtern

## 📚 Dokumentation

- **Vollständig:** `docs/control_qc_readme.md`
- **Beispiele:** `examples/control_qc_example.py`
- **Tests:** `tests/test_control_qc.py`
- **Zusammenfassung:** `IMPLEMENTATION_SUMMARY.md`

## ✨ Features

✅ Alle gewünschten QC Checks implementiert  
✅ Sauberer, dokumentierter, type-hinted Code  
✅ PyPipeGraph2 Integration  
✅ Schöne Plots (seaborn/matplotlib)  
✅ Statistische Tests (Wilcoxon, PCA, Korrelationen)  
✅ Umfangreiche Dokumentation  
✅ Interpretations-Richtlinien  
✅ Sofort einsatzbereit  

Viel Erfolg mit deiner CRISPR Screen Analyse! 🎉
