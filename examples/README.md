# Examples

This directory contains example scripts demonstrating how to use the crispr_screens package.

## Available Examples

### example_new_reports.py

Demonstrates the new report generation classes:
- `QCReport`: Comprehensive QC analysis with control validation and normalization comparison
- `ResultReport`: MAGeCK result reporting with PDF export

**Features shown:**
- Running comprehensive QC analysis
- Validating control sgRNAs
- Comparing multiple normalization methods
- Generating PDF reports
- Creating result reports from MAGeCK outputs
- Exporting ranked gene lists for GSEA

**Usage:**
```bash
python examples/example_new_reports.py
```

**Customization:**
```python
from crispr_screens.models import QCReport, QCConfig

config = QCConfig(
    project_name="My Screen",
    out_dir="my_qc",
    baseline_condition="Total",
    norm_methods=["total", "median", "stable_set"]
)

qc = QCReport(
    config=config,
    count_table_path="my_counts.tsv",
    control_sgrnas_path="my_controls.txt"
)

qc.build(generate_pdf=True)
```

## Migration Examples

### Old Approach (Deprecated)

```python
from crispr_screens.core.report import ResultReport, ReportConfig

# Old import still works but shows deprecation warning
```

### New Approach (Recommended)

```python
from crispr_screens.models import ResultReport, ReportConfig

# New import from models module
```

## PyPipeGraph2 Integration

### Using comprehensive_qc_job

```python
from crispr_screens import comprehensive_qc_job

job = comprehensive_qc_job(
    count_table="counts.tsv",
    baseline_condition="Total",
    output_dir="results/qc",
    control_sgrnas="controls.txt",
    run_control_qc=True,
    run_library_qc=True,
    generate_pdf=True
)
```

### Complete Pipeline

```python
from crispr_screens import (
    mageck_count_job,
    comprehensive_qc_job,
    mageck_mle_job,
    mageck_report_job
)

# Step 1: Count reads
count_job = mageck_count_job(
    fastq_files=fastq_files,
    library="brunello.csv",
    output_prefix="results/counts"
)

# Step 2: QC analysis
qc_job = comprehensive_qc_job(
    count_table="results/counts.count.tsv",
    baseline_condition="Total",
    output_dir="results/qc",
    control_sgrnas="controls.txt",
    generate_pdf=True,
    dependencies=[count_job]
)

# Step 3: MAGeCK MLE
mle_job = mageck_mle_job(
    count_file="results/counts.count.tsv",
    design_matrix="design.tsv",
    output_prefix="results/mle",
    dependencies=[count_job, qc_job]
)

# Step 4: Generate report
report_job = mageck_report_job(
    gene_summary_path="results/mle.gene_summary.txt",
    output_dir="results/report",
    qc_json_path="results/qc/qc_summary.json",
    dependencies=[mle_job]
)
```

## Output Files

### QC Report Outputs

```
qc_results/
├── qc_summary.md              # Markdown report
├── qc_summary.json            # Machine-readable summary
├── qc_summary.html            # HTML report
├── qc_report.pdf              # PDF report
└── qc_assets/
    ├── plots/
    │   ├── control_median_shift.png
    │   ├── control_iqr.png
    │   ├── replicate_correlations_comparison.png
    │   └── logfc_distribution_comparison.png
    └── tables/
        ├── library_stats.tsv
        ├── size_factors.tsv
        ├── size_factors_comparison.tsv
        ├── control_neutrality_qc.json
        ├── best_normalization.json
        └── analysis_recommendation.json
```

### Result Report Outputs

```
report_results/
├── report.md                   # Markdown report
├── report.html                 # HTML report
├── result_report.pdf           # PDF report
└── report_assets/
    ├── plots/
    │   ├── volcano_Time21.png
    │   ├── waterfall_Time21.png
    │   └── effect_vs_reproducibility_Time21.png
    └── tables/
        ├── top_hits_Time21.tsv
        └── ranklist_Time21.rnk
```

## Requirements

All examples require:
- Python >= 3.12
- crispr_screens package installed
- reportlab (for PDF generation)
- markdown (for HTML generation)

Install dependencies:
```bash
pip install crispr_screens[pdf]
```

Or if installing from source:
```bash
cd code/crispr_screens
pip install -e .
```

## See Also

- [API Documentation](../docs/api.md)
- [Migration Guide](../docs/MIGRATION.md)
- [Control QC README](../docs/control_qc_readme.md)
