# Control sgRNA QC - File Index

## 📚 Documentation & Guides

### English
- **[IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)** - Complete implementation summary with all features
- **[docs/control_qc_readme.md](docs/control_qc_readme.md)** - Comprehensive documentation (installation, usage, interpretation)
- **[ARCHITECTURE.txt](ARCHITECTURE.txt)** - Visual architecture overview with data flow diagrams

### German
- **[KURZANLEITUNG_DE.md](KURZANLEITUNG_DE.md)** - Deutsche Kurzanleitung (Quick-Start auf Deutsch)

---

## 🎯 Quick Start Files

- **[quickstart_control_qc.py](quickstart_control_qc.py)** ⭐ 
  - Executable script to run QC on your data
  - Includes data validation and helpful output
  - Usage: `python quickstart_control_qc.py`

- **[examples/control_qc_example.py](examples/control_qc_example.py)**
  - Multiple usage examples
  - Standalone, programmatic, and PyPipeGraph integration
  - Interpretation guidelines

---

## 💻 Source Code

### Core Implementation
- **[src/crispr_screens/core/qc.py](src/crispr_screens/core/qc.py)** - Main QC module (700+ lines)
  - `control_sgrna_qc()` - Core analysis function
  - `generate_control_qc_report()` - High-level report generator
  - `plot_control_distribution_per_condition()` - Univariate plots
  - `plot_pairwise_control_shifts()` - Pairwise heatmap
  - `plot_control_replicate_correlation()` - Correlation matrices
  - `plot_control_pca()` - PCA visualization
  - Helper functions for data loading and processing

### Service Layer
- **[src/crispr_screens/services/io.py](src/crispr_screens/services/io.py)**
  - `control_qc_report()` - I/O wrapper for service layer

### PyPipeGraph Integration
- **[src/crispr_screens/jobs/qc_jobs.py](src/crispr_screens/jobs/qc_jobs.py)**
  - `control_qc_job()` - PyPipeGraph2 job wrapper
  - Handles dependencies, output files, and invariants

### Package Configuration
- **[src/crispr_screens/__init__.py](src/crispr_screens/__init__.py)**
  - Exports `control_qc_job` for easy import
  - Package-level API

---

## 🧪 Testing

- **[tests/test_control_qc.py](tests/test_control_qc.py)**
  - Functional tests with mock data
  - Tests imports, analysis, and report generation
  - Run with: `python tests/test_control_qc.py`

---

## ⚙️ Configuration Files

### Python Package
- **[pyproject.toml](pyproject.toml)**
  - Updated with new dependencies:
    - pandas, numpy, matplotlib, seaborn
    - scipy, scikit-learn

### Anysnake2 Environment
- **[../../anysnake2.toml](../../anysnake2.toml)**
  - Updated with:
    - `crispr_screens = {editable = true, path = "./code/crispr_screens"}`
    - Added: pandas, scipy, seaborn, scikit-learn

---

## 📊 Expected Output Structure

After running QC, you'll get:

```
results/qc/control_sgrnas/
├── control_qc_metrics.tsv                 # Metrics per condition
├── control_qc_pairwise_shifts.tsv         # Pairwise comparison matrix
├── control_qc_distribution.png            # Univariate histograms
├── control_qc_distribution.pdf
├── control_qc_pairwise_heatmap.png        # Pairwise heatmap
├── control_qc_pairwise_heatmap.pdf
├── control_qc_replicate_correlation.png   # Replicate consistency
├── control_qc_replicate_correlation.pdf
├── control_qc_pca_condition.png           # PCA by condition
├── control_qc_pca_condition.pdf
├── control_qc_pca_replicate.png           # PCA by replicate
└── control_qc_pca_replicate.pdf
```

---

## 🚀 Getting Started

### Step 1: Update Environment
```bash
cd /clara/ffs/e/20260105_AG_Stiewe_Andrea_Nist_sgRNA_Brunello_screen
anysnake2 shell
```

### Step 2: Choose Your Method

**Option A: Quick Start Script (Recommended for First Run)**
```bash
cd code/crispr_screens
python quickstart_control_qc.py
```

**Option B: Direct Python**
```python
from crispr_screens.core.qc import generate_control_qc_report

report = generate_control_qc_report(
    count_table="results/mageck_count/all/counts.count.tsv",
    control_sgrnas="incoming/control_sgRNAs.txt",
    baseline_condition="Total",
    output_dir="results/qc/control_sgrnas",
)
```

**Option C: Integrate into Workflow**
```python
from crispr_screens import control_qc_job

qc_job = control_qc_job(
    count_table="results/mageck_count/all/counts.count.tsv",
    control_sgrnas="incoming/control_sgRNAs.txt",
    baseline_condition="Total",
    output_dir="results/qc",
    dependencies=[count_job],
)
```

### Step 3: Interpret Results

1. **Check metrics file** - Look for median Δc close to 0
2. **View PCA plots** - Controls should NOT separate by condition
3. **Check distribution plots** - Should be centered at 0
4. **Verify replicate correlation** - Should be > 0.9

See [docs/control_qc_readme.md](docs/control_qc_readme.md) for detailed interpretation.

---

## 📖 Key Concepts

### What is Δc?
```
Δc = log2(CPM_condition / CPM_baseline)
```
- Measures log2 fold-change of controls vs baseline
- Should be ~0 if controls are neutral
- Large values indicate systematic shifts

### Why Check Controls?
Controls are used for normalization in MAGeCK. If they:
- Show systematic shifts → Normalization introduces bias
- Are unstable → Increases noise
- Respond to treatment → Can't be used as neutral reference

### Multi-Condition Specifics
With >2 conditions, need to check:
1. Per-condition stability (vs baseline)
2. Pairwise consistency (no global drifts)
3. Replicate reproducibility (within conditions)
4. Multivariate patterns (PCA for hidden batch effects)

---

## 🔍 Troubleshooting

### Common Issues

**"Module not found"**
- Solution: Run `anysnake2 shell` to rebuild environment

**"Baseline condition not found"**
- Check column names: `pd.read_csv(..., sep='\t').columns`
- Adjust `baseline_condition` parameter

**"No control sgRNAs found"**
- Verify sgRNA IDs match between files
- Check column name (default: "sgRNA")

**Plotting errors**
- Ensure matplotlib backend is set (use 'Agg' for servers)
- Check for sufficient memory with large datasets

---

## 📞 Support

For questions or issues:
1. Check [docs/control_qc_readme.md](docs/control_qc_readme.md) - Comprehensive documentation
2. Review [examples/control_qc_example.py](examples/control_qc_example.py) - Usage examples
3. Run [tests/test_control_qc.py](tests/test_control_qc.py) - Verify installation

---

## ✨ Features at a Glance

✅ **5 QC Checks**: Univariate, Pairwise, Replicate, PCA, Statistical  
✅ **Rich Visualizations**: Histograms, heatmaps, PCA biplots  
✅ **Statistical Tests**: Wilcoxon, correlations, effect sizes  
✅ **Multi-format Output**: TSV metrics + PNG/PDF plots  
✅ **Workflow Integration**: PyPipeGraph2 jobs  
✅ **Flexible Usage**: Standalone, programmatic, or workflow  
✅ **Comprehensive Documentation**: English + German  
✅ **Production Ready**: Type hints, error handling, testing  

---

## 📝 Quick Reference

| Need | File |
|------|------|
| Quick start | [quickstart_control_qc.py](quickstart_control_qc.py) |
| Full documentation | [docs/control_qc_readme.md](docs/control_qc_readme.md) |
| Examples | [examples/control_qc_example.py](examples/control_qc_example.py) |
| German guide | [KURZANLEITUNG_DE.md](KURZANLEITUNG_DE.md) |
| Architecture | [ARCHITECTURE.txt](ARCHITECTURE.txt) |
| Source code | [src/crispr_screens/core/qc.py](src/crispr_screens/core/qc.py) |
| Tests | [tests/test_control_qc.py](tests/test_control_qc.py) |

---

## 📅 Version Info

**Implementation Date**: January 8, 2026  
**Python Version**: 3.12+  
**Key Dependencies**: pandas, numpy, matplotlib, seaborn, scipy, scikit-learn  
**Integration**: PyPipeGraph2, Anysnake2  

---

Viel Erfolg! 🎉
