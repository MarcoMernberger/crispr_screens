# API Reference

## Models

The models module contains report generators and configuration classes.

### ResultReport

Comprehensive report generator for MAGeCK CRISPR screen results (MLE/RRA).

```python
from crispr_screens.models import ResultReport, ReportConfig

config = ReportConfig(
    project_name="My Screen",
    out_dir="report_out",
    fdr_threshold=0.1,
    effect_threshold=1.0
)

report = ResultReport(
    config=config,
    gene_summary_path="gene_summary.tsv",
    sgrna_summary_path="sgrna_summary.tsv",  # optional
    count_path="counts.tsv",  # optional
    metadata_path="metadata.tsv",  # optional
    qc_json_path="qc_summary.json"  # optional
)

# Generate reports
report.build(generate_pdf=True)  # Generates MD, HTML, and PDF

# Export ranked gene list for GSEA
report.export_ranklist(eff="Time21")
```

**Output Files:**
- `report.md`: Markdown report
- `report.html`: HTML report (if markdown available)
- `result_report.pdf`: PDF report (if generate_pdf=True)
- `report_assets/plots/*.png`: Volcano, waterfall, and other plots
- `report_assets/tables/*.tsv`: Top hits tables

### QCReport

Comprehensive quality control report generator.

```python
from crispr_screens.models import QCReport, QCConfig

config = QCConfig(
    project_name="My Screen",
    out_dir="qc_out",
    baseline_condition="Total",
    norm_methods=["total", "median", "stable_set"]
)

qc = QCReport(
    config=config,
    count_table_path="counts.tsv",
    control_sgrnas_path="controls.txt",  # optional
    metadata_path="metadata.tsv"  # optional
)

# Run QC analysis
qc.build(
    run_control_qc=True,  # Validate control sgRNAs
    run_library_qc=True,  # Compare normalization methods
    generate_pdf=True  # Generate PDF report
)
```

**Output Files:**
- `qc_summary.md`: Markdown QC report
- `qc_summary.json`: Machine-readable summary
- `qc_summary.html`: HTML report
- `qc_report.pdf`: PDF report (if generate_pdf=True)
- `qc_assets/plots/*.png`: QC plots
- `qc_assets/tables/*.tsv`: QC metrics and size factors

## Jobs

PyPipeGraph2 job wrappers for MAGeCK workflows.

### comprehensive_qc_job

Unified QC job that replaces control_qc_job and standard_qc_job.

```python
from crispr_screens import comprehensive_qc_job

job = comprehensive_qc_job(
    count_table="counts.tsv",
    baseline_condition="Total",
    output_dir="results/qc",
    control_sgrnas="controls.txt",
    project_name="My Screen",
    run_control_qc=True,
    run_library_qc=True,
    generate_pdf=True
)
```

### mageck_mle_job

Run MAGeCK MLE analysis.

```python
from crispr_screens import mageck_mle_job

job = mageck_mle_job(
    count_file="counts.tsv",
    design_matrix="design.tsv",
    output_prefix="results/mle"
)
```

### mageck_rra_job

Run MAGeCK RRA analysis.

```python
from crispr_screens import mageck_rra_job

job = mageck_rra_job(
    count_file="counts.tsv",
    control_samples=["Total_Rep1", "Total_Rep2"],
    treatment_samples=["Sorted_Rep1", "Sorted_Rep2"],
    output_prefix="results/rra"
)
```

## Core Functions

Low-level functions for custom workflows.

### QC Functions

```python
from crispr_screens.core.qc import (
    load_control_sgrnas,
    read_counts,
    calculate_cpm,
    calculate_delta_logfc,
    qc_controls_neutrality,
    choose_best_normalization,
    recommend_analysis_method
)
```

### Plotting Functions

```python
from crispr_screens.core.plots import (
    volcano_plot,
    ma_plot,
    plot_replicate_correlation,
    plot_pca
)
```

## Migration from core.report

The `ResultReport` class has been moved from `core.report` to `models`:

**Old (deprecated):**
```python
from crispr_screens.core.report import ResultReport, ReportConfig
```

**New:**
```python
from crispr_screens.models import ResultReport, ReportConfig
```

See [MIGRATION.md](MIGRATION.md) for detailed migration guide.
