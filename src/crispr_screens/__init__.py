"""
CRISPR screens – demultiplexing and read processing tools.

This package provides:
- core
- models
- services
- jobs
- api
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("crispr-screens")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.1.0"


from .jobs.mageck_jobs import (
    combine_comparison_output_job,
    write_filtered_mageck_comparison_job,
    run_mageck_scatterview_job,
    mageck_mle_job,
    mageck_count_job,
    mageck_rra_job,
    create_query_control_sgrna_frames_job,
    create_combine_gene_info_with_mageck_output_job,
    create_spiked_count_table_job,
    evaluate_spike_in_performance_job,
    spike_evaluation_report_job,
    write_significant_genes_rra_job,
)
from .jobs.plot_jobs import write_venn_job, write_volcano_plot_job
from .jobs.qc_jobs import (
    control_qc_job,
    standard_qc_job,
    mageck_report_job,
    pairing_qc_job,
    pairing_qc_plots_job,
    # comprehensive_qc_job,
)

# from .models import ReportConfig, ResultReport, QCConfig, QCReport

__all__ = [
    "combine_comparison_output_job",
    "write_filtered_mageck_comparison_job",
    "run_mageck_scatterview_job",
    "write_venn_job",
    "mageck_mle_job",
    "mageck_count_job",
    "mageck_rra_job",
    "create_query_control_sgrna_frames_job",
    "create_combine_gene_info_with_mageck_output_job",
    "write_volcano_plot_job",
    "control_qc_job",
    "standard_qc_job",
    "comprehensive_qc_job",
    "mageck_report_job",
    "ReportConfig",
    "ResultReport",
    "QCConfig",
    "QCReport",
    "pairing_qc_job",
    "pairing_qc_plots_job",
    "create_spiked_count_table_job",
    "evaluate_spike_in_performance_job",
    "spike_evaluation_report_job",
    "write_significant_genes_rra_job",
]
