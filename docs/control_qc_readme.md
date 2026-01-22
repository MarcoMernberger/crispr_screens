# Control sgRNA Quality Control for CRISPR Screens

Comprehensive QC analysis to validate that control sgRNAs are stable and suitable for normalization across multiple experimental conditions.

## Overview

When analyzing CRISPR screens with MAGeCK, control sgRNAs (non-targeting controls) are used for normalization. This QC module validates that these controls are:

1. **Neutral** - No systematic shifts between conditions
2. **Stable** - Consistent across replicates
3. **Suitable for normalization** - No treatment-specific response

## Features

### Implemented QC Checks

1. **Univariate Analysis** (per condition vs baseline)
   - Histogram/density plots of Δc = log2(CPM_condition / CPM_baseline)
   - Metrics: median, IQR, tail-rate (|Δc| > 1)
   - Wilcoxon signed-rank test for systematic shift

2. **Pairwise Comparisons** (condition vs condition)
   - Heatmap of median control shifts between all condition pairs
   - Identifies global drifts or batch effects

3. **Replicate Consistency**
   - Correlation heatmaps between replicates per condition
   - Pearson correlation on log2(CPM+1) values

4. **Multivariate PCA**
   - PCA on control sgRNAs across all samples
   - Visualize whether controls cluster by condition (bad) or replicate (good)
   - Scree plot of variance explained

5. **Statistical Testing**
   - Wilcoxon signed-rank test per condition (H0: median Δc = 0)
   - Effect size metrics (not just p-values)

## Installation

The module requires the following dependencies (automatically added to `pyproject.toml`):

```toml
dependencies = [
  "pandas>=2.2.2,<3.0.0",
  "numpy>=1.26.0,<3.0.0",
  "matplotlib>=3.8.0,<4.0.0",
  "seaborn>=0.13.0,<0.14.0",
  "scipy>=1.11.0,<2.0.0",
  "scikit-learn>=1.5.0,<2.0.0",
]
```

## Usage

### Standalone Analysis

```python
from crispr_screens.core.qc import generate_control_qc_report

# Generate full QC report
report = generate_control_qc_report(
    count_table="results/mageck_count/all/counts.count.tsv",
    control_sgrnas="incoming/control_sgRNAs.txt",
    baseline_condition="Total",  # or "T0", "unsorted", etc.
    output_dir="results/qc/control_sgrnas",
    prefix="control_qc",
)

# Access results
metrics = report["qc_results"]["metrics"]
for condition, metric in metrics.items():
    print(f"{condition}: median Δc = {metric['median']:.3f}")
```

### PyPipeGraph Integration

```python
from crispr_screens.jobs.qc_jobs import control_qc_job
from crispr_screens.jobs.mageck_jobs import mageck_count_job

# Create MAGeCK count job
count_job = mageck_count_job(...)

# Create control QC job (depends on count job)
qc_job = control_qc_job(
    count_table="results/mageck_count/all/counts.count.tsv",
    control_sgrnas="incoming/control_sgRNAs.txt",
    baseline_condition="Total",
    output_dir="results/qc/control_sgrnas",
    dependencies=[count_job],
)
```

### Programmatic Access

```python
from crispr_screens.core.qc import control_sgrna_qc

# Run QC analysis
qc_results = control_sgrna_qc(
    count_table="counts.tsv",
    control_sgrnas="controls.txt",
    baseline_condition="Total",
)

# Access components
qc_results["metrics"]  # Per-condition metrics
qc_results["pairwise_median"]  # Pairwise shifts
qc_results["replicate_correlations"]  # Replicate consistency
qc_results["delta"]  # Raw Δc values
qc_results["cpm"]  # CPM values
```

## Input Files

### Count Table Format

MAGeCK count table (`.count.tsv`) with structure:

```
sgRNA       Gene    Total_Rep1  Total_Rep2  Sort1_Rep1  Sort1_Rep2  ...
s_8567      HNMT    1195        1211        38          511         ...
s_76462     NonTargetingControl  1408  1294  75        1765        ...
```

- **sgRNA column**: sgRNA IDs (must match control file)
- **Gene column**: Gene names
- **Sample columns**: Raw read counts (format: `Condition_Replicate`)

### Control sgRNA File Format

Text file with one control sgRNA ID per line:

```
s_76442
s_76443
s_76444
...
```

## Output Files

The QC report generates:

1. **Metrics TSV** (`control_qc_metrics.tsv`)
   - Per-condition: median, IQR, tail-rate, Wilcoxon p-value

2. **Pairwise Shifts TSV** (`control_qc_pairwise_shifts.tsv`)
   - Matrix of median shifts between all condition pairs

3. **Distribution Plots** (`control_qc_distribution.png/pdf`)
   - Histograms with density overlays per condition

4. **Pairwise Heatmap** (`control_qc_pairwise_heatmap.png/pdf`)
   - Heatmap of condition-condition shifts

5. **Replicate Correlation** (`control_qc_replicate_correlation.png/pdf`)
   - Correlation matrices per condition

6. **PCA Plots** (`control_qc_pca_condition.png/pdf`, `control_qc_pca_replicate.png/pdf`)
   - PCA colored by condition and by replicate

## Interpretation Guidelines

### ✓ GOOD Controls (suitable for normalization)

- **Median Δc ≈ 0** (± 0.1)
- **Tail rate < 5%** (|Δc| > 1)
- **Replicate correlation > 0.9**
- **PCA**: Samples cluster by replicate, NOT by condition
- **Pairwise shifts**: All ≈ 0

### ⚠ ACCEPTABLE Controls

- **Median Δc**: 0.1-0.3
- **Tail rate**: 5-10%
- **Replicate correlation**: 0.8-0.9
- Some technical variation but consistent

### ✗ PROBLEMATIC Controls (NOT suitable)

- **Median Δc > 0.5** → Systematic shift
- **Tail rate > 15%** → High instability
- **Replicate correlation < 0.8** → Batch effects
- **PCA**: Conditions separate on PC1 → Treatment response
- **Pairwise**: One condition shows consistent offset → Global drift

### Decision Tree

```
IF all checks pass:
   → Use control-based normalization (MAGeCK default)

ELSE IF one condition problematic:
   → Exclude that condition OR use median normalization

ELSE IF controls show treatment response:
   → Use different control set OR median/total-count normalization
```

## Examples

See `examples/control_qc_example.py` for:
- Standalone analysis
- Programmatic access
- PyPipeGraph integration
- Custom baseline detection
- Detailed interpretation guide

## Implementation Details

### Module Structure

```
crispr_screens/
├── core/
│   └── qc.py                    # Core QC functions and plotting
├── services/
│   └── io.py                    # I/O wrapper (control_qc_report)
└── jobs/
    └── qc_jobs.py               # PyPipeGraph job wrapper
```

### Key Functions

#### `control_sgrna_qc()`
Core analysis function that:
1. Loads count table and filters for control sgRNAs
2. Calculates CPM from raw counts (NOT using MAGeCK-normalized values)
3. Computes Δc = log2(CPM_condition + 1) - log2(CPM_baseline + 1)
4. Calculates metrics per condition
5. Performs pairwise comparisons
6. Computes replicate correlations
7. Returns comprehensive results dict

#### Plotting Functions
- `plot_control_distribution_per_condition()` - Univariate histograms
- `plot_pairwise_control_shifts()` - Pairwise heatmap
- `plot_control_replicate_correlation()` - Correlation matrices
- `plot_control_pca()` - PCA biplot with variance

#### `generate_control_qc_report()`
High-level function that runs all analyses and saves outputs.

### Normalization Approach

**Important**: Uses **raw counts** to calculate CPM independently:

```python
CPM = (count / library_size) × 10^6
```

Does NOT use MAGeCK-normalized counts (which may already be control-normalized), as this would be circular reasoning.

### Baseline Selection

The baseline condition (reference) must be specified:
- Common choices: `"Total"`, `"T0"`, `"unsorted"`, `"Plasmid"`
- All other conditions are compared to this baseline
- Must match condition name in column headers (before delimiter)

### Column Name Parsing

By default, splits column names on `"_"`:
- `Total_Rep1` → condition: `"Total"`, replicate: `"Rep1"`
- `Sort1_Rep2` → condition: `"Sort1"`, replicate: `"Rep2"`

Customize with `delimiter` parameter if using different naming scheme.

## Theory and Background

### Why This QC Matters

MAGeCK's default normalization uses control sgRNAs to estimate size factors. If controls:
- Show systematic shifts → Normalization introduces bias
- Are unstable → Increases noise
- Respond to treatment → Cannot be used as "neutral" reference

### Multi-Condition Specifics

For screens with >2 conditions, QC must check:
1. **Per-condition** stability (vs baseline)
2. **Pairwise** consistency (no global drifts)
3. **Replicate** reproducibility (within conditions)
4. **Multivariate** patterns (PCA reveals hidden batch effects)

### Statistical Considerations

- **Δc values**: Log2 fold-change is symmetric and biologically interpretable
- **Wilcoxon test**: Non-parametric test for median shift (robust to outliers)
- **Effect size**: Always check median/IQR, not just p-values
- **Multiple testing**: Not corrected (exploratory QC, not hypothesis testing)

## Troubleshooting

### Common Issues

1. **"No control sgRNAs found"**
   - Check sgRNA IDs match between count table and control file
   - Verify column name: default is `"sgRNA"`

2. **"Baseline condition not found"**
   - Check condition name matches column prefix
   - List available: `pd.read_csv(count_file, sep='\t').columns`

3. **PCA fails / empty plots**
   - Ensure >3 samples for meaningful PCA
   - Check for zero-variance sgRNAs (filtered automatically)

4. **High memory usage**
   - Large screens (>100K sgRNAs) may need substantial RAM
   - Filter count table to controls before loading if needed

## References

This implementation follows best practices from:
- Li et al. (2014) MAGeCK paper
- Hart et al. (2015) High-resolution CRISPR screens
- Doench et al. (2016) Optimized sgRNA design
- Community standards for CRISPR screen QC

## License

Same as parent package (see main LICENSE file).

## Authors

Implementation: GitHub Copilot & Marco Mernberger
Theory/Design: Based on user requirements and CRISPR screen best practices
