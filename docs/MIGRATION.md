# Migration Guide: core.report → models

## Overview

The report classes have been moved from `crispr_screens.core.report` to `crispr_screens.models` as part of a refactoring to improve code organization and add new functionality.

## What Changed

### 1. Module Location
- **Old:** `from crispr_screens.core.report import ResultReport, ReportConfig`
- **New:** `from crispr_screens.models import ResultReport, ReportConfig`

### 2. New Classes Added
- `QCReport`: Unified QC reporting class that combines control and library QC
- `QCConfig`: Configuration for QC reports
- `comprehensive_qc_job`: New PyPipeGraph2 job for unified QC

### 3. Enhanced Functionality
- **PDF Export**: Both `ResultReport` and `QCReport` now support PDF generation via reportlab
- **Unified QC**: `QCReport` combines control_qc and standard_qc functionality
- **Better Organization**: Report classes are now in `models/` where they conceptually belong

## Migration Steps

### Step 1: Update Imports

**Before:**
```python
from crispr_screens.core.report import ResultReport, ReportConfig
```

**After:**
```python
from crispr_screens.models import ResultReport, ReportConfig
```

### Step 2: Optional - Use New PDF Export

**Generate PDF with ResultReport:**
```python
from crispr_screens.models import ResultReport, ReportConfig

config = ReportConfig(project_name="My Screen", out_dir="results")
report = ResultReport(config, gene_summary_path="gene_summary.tsv")
report.build(generate_pdf=True)  # NEW: generate_pdf parameter
```

**Generate PDF with QCReport:**
```python
from crispr_screens.models import QCReport, QCConfig

config = QCConfig(project_name="My Screen", out_dir="qc_results")
qc = QCReport(config, count_table_path="counts.tsv", control_sgrnas_path="controls.txt")
qc.build(generate_pdf=True)
```

### Step 3: Optional - Use New Unified QC

**Replace separate control_qc_job and standard_qc_job:**

**Before:**
```python
from crispr_screens import control_qc_job, standard_qc_job

control_job = control_qc_job(
    count_table="counts.tsv",
    control_sgrnas="controls.txt",
    baseline_condition="Total",
    output_dir="results/control_qc"
)

qc_job = standard_qc_job(
    count_table="counts.tsv",
    baseline_condition="Total",
    output_dir="results/qc"
)
```

**After:**
```python
from crispr_screens import comprehensive_qc_job

# Single job that does both control and library QC
qc_job = comprehensive_qc_job(
    count_table="counts.tsv",
    baseline_condition="Total",
    output_dir="results/qc",
    control_sgrnas="controls.txt",
    run_control_qc=True,
    run_library_qc=True,
    generate_pdf=True
)
```

## Backward Compatibility

The old import path still works but will show a deprecation warning:

```python
# Still works, but shows warning
from crispr_screens.core.report import ResultReport

# DeprecationWarning: The 'crispr_screens.core.report' module is deprecated.
# Please use 'crispr_screens.models' instead
```

**Recommendation:** Update your imports to avoid warnings and benefit from improved organization.

## New Dependencies

The following packages are now required for PDF generation:
- `reportlab>=4.0.0` (for PDF generation)
- `markdown>=3.5.0` (for HTML report generation)

These are automatically installed as dependencies.

## Benefits of Migration

1. **Better Organization**: Report classes are in `models/` where they conceptually belong
2. **PDF Export**: Professional PDF reports with embedded plots and tables
3. **Unified QC**: Single job for all QC analysis instead of multiple separate jobs
4. **Cleaner API**: Clear separation between core functionality and report generation
5. **Future-Proof**: New development will focus on the `models` module

## Support

If you encounter issues during migration, please:
1. Check that all imports are updated
2. Ensure `reportlab` and `markdown` are installed
3. Verify that `build()` method calls include correct parameters

For questions, contact the development team.
