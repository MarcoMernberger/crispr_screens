"""
PyPipeGraph2 job wrappers for CRISPR screen quality control.
"""

from importlib import util as _import_util
from pypipegraph2 import (
    Job,
    MultiFileGeneratingJob,
    FunctionInvariant,
    ParameterInvariant,
)
from pathlib import Path
from typing import List, Union
from crispr_screens.services.io import (
    generate_control_qc_report,
    standard_qc_report,
    mageck_report,
)
from crispr_screens.core.qc import (
    export_control_counts_and_cpm,
    generate_standard_qc_report,
    compute_size_factors_total,
    compute_size_factors_median_ratio,
    compute_size_factors_stable_set,
    qc_logfc_distribution,
    qc_replicate_consistency,
    choose_best_normalization,
    recommend_analysis_method,
)


def control_qc_job(
    count_table: Union[Path, str],
    control_sgrnas: Union[Path, str],
    baseline_condition: str,
    output_dir: Union[Path, str],
    prefix: str = "control_qc",
    sgrna_col: str = "sgRNA",
    gene_col: str = "Gene",
    delimiter: str = "_",
    save_formats: List[str] = ["png", "pdf"],
    dependencies: List[Job] = [],
) -> MultiFileGeneratingJob:
    """
    Create pypipegraph job for control sgRNA QC analysis.

    Generates comprehensive QC report validating that control sgRNAs are
    stable and suitable for normalization across experimental conditions.

    Parameters
    ----------
    count_table : Path or str
        Path to MAGeCK count table (.count.tsv).
    control_sgrnas : Path or str
        Path to control sgRNA file (one ID per line).
    baseline_condition : str
        Baseline condition name (e.g., "Total", "T0", "unsorted").
    output_dir : Path or str
        Output directory for QC reports.
    prefix : str
        Filename prefix for output files.
    sgrna_col : str
        sgRNA column name in count table.
    gene_col : str
        Gene column name in count table.
    delimiter : str
        Delimiter for parsing condition_replicate from column names.
    save_formats : list
        List of formats to save ("png", "pdf", "svg").
    dependencies : list
        List of pypipegraph Jobs to depend on.

    Returns
    -------
    MultiFileGeneratingJob
        Job that generates QC metrics and plots.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Define output files
    outfiles = [
        output_dir / f"{prefix}_metrics.tsv",
        output_dir / f"{prefix}_pairwise_shifts.tsv",
        output_dir / f"{prefix}_raw_counts.tsv",
        output_dir / f"{prefix}_cpm.tsv",
    ]

    # Add plot files
    plot_names = [
        "distribution",
        "pairwise_heatmap",
        "replicate_correlation",
        "pca_condition",
        "pca_replicate",
    ]

    for plot_name in plot_names:
        for fmt in save_formats:
            outfiles.append(output_dir / f"{prefix}_{plot_name}.{fmt}")

    def __dump(
        outfiles,
        count_table=count_table,
        control_sgrnas=control_sgrnas,
        baseline_condition=baseline_condition,
        output_dir=output_dir,
        prefix=prefix,
        sgrna_col=sgrna_col,
        gene_col=gene_col,
        delimiter=delimiter,
        save_formats=save_formats,
    ):
        generate_control_qc_report(
            count_table=count_table,
            control_sgrnas=control_sgrnas,
            baseline_condition=baseline_condition,
            output_dir=output_dir,
            prefix=prefix,
            sgrna_col=sgrna_col,
            gene_col=gene_col,
            delimiter=delimiter,
            save_formats=save_formats,
        )

    job = MultiFileGeneratingJob(outfiles, __dump).depends_on(dependencies)

    # Add invariants
    job.depends_on(
        FunctionInvariant(
            f"{prefix}_export_control_counts_and_cpm_func",
            export_control_counts_and_cpm,
        )
    )
    job.depends_on(
        FunctionInvariant(
            f"{prefix}_generate_control_qc_report_func",
            generate_control_qc_report,
            allowed_non_locals=[
                "control_sgrna_qc",
                "export_control_counts_and_cpm",
                "plot_control_distribution_per_condition",
                "plot_control_pca",
                "plot_control_pca",
                "plot_control_replicate_correlation",
                "plot_pairwise_control_shifts",
                "plt",
            ],
        )
    )
    job.depends_on(
        ParameterInvariant(
            f"{prefix}_control_qc_params",
            (
                str(count_table),
                str(control_sgrnas),
                baseline_condition,
                sgrna_col,
                gene_col,
                delimiter,
                tuple(save_formats),
            ),
        )
    )

    return job


# def export_control_data_job(
#     count_table: Union[Path, str],
#     control_sgrnas: Union[Path, str],
#     output_dir: Union[Path, str],
#     prefix: str = "control_data",
#     sgrna_col: str = "sgRNA",
#     gene_col: str = "Gene",
#     dependencies: List[Job] = [],
# ) -> MultiFileGeneratingJob:
#     """
#     Create pypipegraph job to export control sgRNA counts and CPM.

#     Exports raw counts and CPM values (calculated on full library) for
#     control sgRNAs only to TSV files.

#     Parameters
#     ----------
#     count_table : Path or str
#         Path to MAGeCK count table (.count.tsv).
#     control_sgrnas : Path or str
#         Path to control sgRNA file (one ID per line).
#     output_dir : Path or str
#         Output directory.
#     prefix : str
#         Filename prefix for output files.
#     sgrna_col : str
#         sgRNA column name in count table.
#     gene_col : str
#         Gene column name in count table.
#     dependencies : list
#         List of pypipegraph Jobs to depend on.

#     Returns
#     -------
#     MultiFileGeneratingJob
#         Job that generates control data TSV files.
#     """
#     output_dir = Path(output_dir)
#     output_dir.mkdir(exist_ok=True, parents=True)

#     # Define output files
#     outfiles = [
#         output_dir / f"{prefix}_raw_counts.tsv",
#         output_dir / f"{prefix}_cpm.tsv",
#     ]

#     def __dump(
#         outfiles,
#         count_table=count_table,
#         control_sgrnas=control_sgrnas,
#         output_dir=output_dir,
#         prefix=prefix,
#         sgrna_col=sgrna_col,
#         gene_col=gene_col,
#     ):
#         export_control_data(
#             count_table=count_table,
#             control_sgrnas=control_sgrnas,
#             output_dir=output_dir,
#             prefix=prefix,
#             sgrna_col=sgrna_col,
#             gene_col=gene_col,
#         )

#     job = MultiFileGeneratingJob(outfiles, __dump).depends_on(dependencies)

#     # Add invariants
#     job.depends_on(
#         FunctionInvariant(
#             f"{prefix}_export_control_data_func", export_control_data
#         )
#     )
#     job.depends_on(
#         ParameterInvariant(
#             f"{prefix}_export_control_data_params",
#             (
#                 str(count_table),
#                 str(control_sgrnas),
#                 sgrna_col,
#                 gene_col,
#             ),
#         )
#     )

#     return job


def standard_qc_job(
    count_tsv: Union[Path, str],
    output_dir: Union[Path, str],
    metadata_tsv: Union[Path, str, None] = None,
    control_sgrna_txt: Union[Path, str, None] = None,
    baseline_condition: str = "total",
    sgrna_col: str = "sgRNA",
    gene_col: str = "Gene",
    delimiter: str = "_",
    norm_methods: Union[List[str], None] = None,
    pseudocount: float = 1.0,
    save_formats: Union[List[str], None] = None,
    dependencies: List[Job] = [],
) -> MultiFileGeneratingJob:
    """
    Create pypipegraph job for standard CRISPR screen QC analysis.

    Generates comprehensive QC report analyzing full library with multiple
    normalization methods and recommending analysis strategy (RRA vs MLE).

    Parameters
    ----------
    count_tsv : Path or str
        Path to MAGeCK count table.
    output_dir : Path or str
        Output directory for QC reports.
    metadata_tsv : Path or str, optional
        Path to metadata file with sample/condition/replicate columns.
        If None, will parse from column names.
    control_sgrna_txt : Path or str, optional
        Path to control sgRNA file. If provided, includes control QC.
    baseline_condition : str
        Baseline condition name (e.g., "total", "T0").
    sgrna_col : str
        sgRNA column name.
    gene_col : str
        Gene column name.
    delimiter : str
        Delimiter for parsing condition_replicate.
    norm_methods : list of str, optional
        Normalization methods to test. Default: ["median", "total",
        "stable_set"].
    pseudocount : float
        Pseudocount for log transformation.
    save_formats : list of str, optional
        Plot formats to save. Default: ["png"].
    dependencies : list
        List of pypipegraph Jobs to depend on.

    Returns
    -------
    MultiFileGeneratingJob
        Job that generates comprehensive QC analysis.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    if save_formats is None:
        save_formats = ["png"]

    if norm_methods is None:
        norm_methods = ["median", "total", "stable_set"]

    # Define output files
    outfiles = [
        output_dir / "library_stats.tsv",
        output_dir / "size_factors.tsv",
        output_dir / "size_factors_comparison.tsv",
        output_dir / "normalization_recommendation.json",
        output_dir / "analysis_recommendation.json",
        output_dir / "qc_summary.md",
    ]

    # Add control QC file if controls provided
    if control_sgrna_txt is not None:
        outfiles.append(output_dir / "control_neutrality_qc.json")

    # Add per-normalization method output directories
    for method in norm_methods:
        method_dir = output_dir / f"norm_{method}"
        outfiles.extend(
            [
                method_dir / "logfc_distribution.tsv",
                method_dir / "replicate_consistency.tsv",
            ]
        )
        for fmt in save_formats:
            outfiles.append(method_dir / f"ma_plots.{fmt}")

    # Add summary plots
    for fmt in save_formats:
        outfiles.extend(
            [
                output_dir / f"pca_full_library.{fmt}",
                output_dir / f"sample_correlations.{fmt}",
            ]
        )

    def __dump(
        outfiles,
        count_tsv=count_tsv,
        output_dir=output_dir,
        metadata_tsv=metadata_tsv,
        control_sgrna_txt=control_sgrna_txt,
        baseline_condition=baseline_condition,
        sgrna_col=sgrna_col,
        gene_col=gene_col,
        delimiter=delimiter,
        norm_methods=norm_methods,
        pseudocount=pseudocount,
        save_formats=save_formats,
    ):
        standard_qc_report(
            count_tsv=count_tsv,
            output_dir=output_dir,
            metadata_tsv=metadata_tsv,
            control_sgrna_txt=control_sgrna_txt,
            baseline_condition=baseline_condition,
            sgrna_col=sgrna_col,
            gene_col=gene_col,
            delimiter=delimiter,
            norm_methods=norm_methods,
            pseudocount=pseudocount,
            save_formats=save_formats,
        )

    job = MultiFileGeneratingJob(outfiles, __dump).depends_on(dependencies)

    # Add invariants for core functions
    job.depends_on(
        FunctionInvariant(
            "standard_qc_generate_standard_qc_report",
            generate_standard_qc_report,
            allowed_non_locals=[
                "calculate_cpm",
                "calculate_delta_logfc",
                "compute_size_factors_total",
                "compute_size_factors_median_ratio",
                "compute_size_factors_stable_set",
                "compute_size_factors_control",
                "qc_logfc_distribution",
                "qc_replicate_consistency",
                "qc_controls_neutrality",
                "choose_best_normalization",
                "recommend_analysis_method",
                "plot_ma_grid",
                "plot_library_pca",
                "plot_sample_correlations",
                "plt",
                "Path",
                "json",
            ],
        )
    )

    job.depends_on(
        FunctionInvariant(
            "standard_qc_compute_size_factors_total",
            compute_size_factors_total,
        )
    )

    job.depends_on(
        FunctionInvariant(
            "standard_qc_compute_size_factors_median_ratio",
            compute_size_factors_median_ratio,
        )
    )

    job.depends_on(
        FunctionInvariant(
            "standard_qc_compute_size_factors_stable_set",
            compute_size_factors_stable_set,
        )
    )

    job.depends_on(
        FunctionInvariant(
            "standard_qc_qc_logfc_distribution",
            qc_logfc_distribution,
        )
    )

    job.depends_on(
        FunctionInvariant(
            "standard_qc_qc_replicate_consistency",
            qc_replicate_consistency,
        )
    )

    job.depends_on(
        FunctionInvariant(
            "standard_qc_choose_best_normalization",
            choose_best_normalization,
        )
    )

    job.depends_on(
        FunctionInvariant(
            "standard_qc_recommend_analysis_method",
            recommend_analysis_method,
        )
    )

    # Add parameter invariant
    job.depends_on(
        ParameterInvariant(
            "standard_qc_params",
            (
                str(count_tsv),
                str(metadata_tsv) if metadata_tsv else None,
                str(control_sgrna_txt) if control_sgrna_txt else None,
                baseline_condition,
                sgrna_col,
                gene_col,
                delimiter,
                tuple(norm_methods),
                pseudocount,
                tuple(save_formats),
            ),
        )
    )

    return job


def mageck_report_job(
    job_id,
    output_dir,
    gene_summary_path=None,
    sgrna_summary_path=None,
    rra_summary_path=None,
    mle_summary_path=None,
    readout="full",
    effect_cols=None,
    gene_col="Gene",
    fdr_threshold=0.25,
    top_n=50,
    gene_sets=None,
    pathway_enrichment_path=None,
    depend_on=None,
    **kwargs,
):
    """
    PyPipeGraph2 job for comprehensive MAGeCK result reporting.

    Creates MultiFileGeneratingJob that produces publication-ready plots
    and analysis metrics for MAGeCK RRA and/or MLE results.

    Parameters
    ----------
    job_id : str
        Unique identifier for the job
    output_dir : str or Path
        Directory to save report outputs
    gene_summary_path : str or Path, optional
        Path to MAGeCK gene summary (for MLE or single RRA)
    sgrna_summary_path : str or Path, optional
        Path to MAGeCK sgRNA summary
    rra_summary_path : str or Path, optional
        Path to MAGeCK RRA gene summary (for comparison with MLE)
    mle_summary_path : str or Path, optional
        Path to MAGeCK MLE gene summary (for comparison with RRA)
    readout : str
        Report mode: "minimal" or "full" (default "full")
    effect_cols : list, optional
        For MLE: list of beta column names for effect decomposition
    gene_col : str
        Column name for gene identifiers
    fdr_threshold : float
        FDR threshold for significance
    top_n : int
        Number of top genes to highlight
    gene_sets : dict, optional
        Dictionary mapping gene set names to lists of genes
    pathway_enrichment_path : str or Path, optional
        Path to pathway enrichment results
    depend_on : Job or list of Jobs, optional
        Jobs to depend on
    **kwargs : dict
        Additional parameters

    Returns
    -------
    MultiFileGeneratingJob
        PyPipeGraph2 job that generates report files
    """
    output_dir = Path(output_dir)

    # Define output files based on readout mode
    outfiles = [output_dir / "mageck_report_summary.json"]

    # Minimal mode plots
    minimal_plots = [
        "01_volcano_plot.png",
        "02_effect_vs_reproducibility.png",
        "04_direction_consistency.png",
    ]

    # Add rank stability if both RRA and MLE provided
    if rra_summary_path and mle_summary_path:
        minimal_plots.append("03_rank_stability.png")

    # Add pathway enrichment if provided
    if pathway_enrichment_path:
        minimal_plots.append("05_pathway_enrichment.png")

    outfiles.extend([output_dir / plot for plot in minimal_plots])

    # Full mode additional plots
    if readout == "full":
        full_plots = [
            "06_replicate_heatmap.png",
            "07_effect_decomposition.png",
            "09_gene_set_distributions.png",
            "10_beta_vs_se.png",
            "11_wald_z_distribution.png",
            "12_qq_plot.png",
        ]

        # Contrast plot for 2-condition comparisons
        if effect_cols and len(effect_cols) == 2:
            full_plots.append("08_contrast_plot.png")

        outfiles.extend([output_dir / plot for plot in full_plots])

    def __dump(
        outfiles,
        output_dir=output_dir,
        gene_summary_path=gene_summary_path,
        sgrna_summary_path=sgrna_summary_path,
        rra_summary_path=rra_summary_path,
        mle_summary_path=mle_summary_path,
        readout=readout,
        effect_cols=effect_cols,
        gene_col=gene_col,
        fdr_threshold=fdr_threshold,
        top_n=top_n,
        gene_sets=gene_sets,
        pathway_enrichment_path=pathway_enrichment_path,
        **kwargs,
    ):
        mageck_report(
            output_dir=output_dir,
            gene_summary_path=gene_summary_path,
            sgrna_summary_path=sgrna_summary_path,
            rra_summary_path=rra_summary_path,
            mle_summary_path=mle_summary_path,
            readout=readout,
            effect_cols=effect_cols,
            gene_col=gene_col,
            fdr_threshold=fdr_threshold,
            top_n=top_n,
            gene_sets=gene_sets,
            pathway_enrichment_path=pathway_enrichment_path,
            **kwargs,
        )

    job = MultiFileGeneratingJob(outfiles, __dump).depends_on(depend_on)
    job.job_id = job_id

    # Add function invariants
    from crispr_screens.core.mageck_report import generate_mageck_report
    from crispr_screens.core.plots import (
        volcano_plot,
        plot_effect_size_vs_reproducibility,
        plot_rank_stability,
        plot_direction_consistency,
        plot_replicate_effect_heatmap,
        plot_effect_decomposition,
        plot_contrast,
        plot_pathway_enrichment_summary,
        plot_gene_set_score_distribution,
        plot_beta_vs_standard_error,
        plot_wald_z_distribution,
        plot_qq,
    )

    job.depends_on(
        FunctionInvariant(
            "mageck_report.generate_mageck_report",
            generate_mageck_report,
        )
    )

    job.depends_on(
        FunctionInvariant(
            "plots.volcano_plot",
            volcano_plot,
        )
    )

    job.depends_on(
        FunctionInvariant(
            "plots.plot_effect_size_vs_reproducibility",
            plot_effect_size_vs_reproducibility,
        )
    )

    job.depends_on(
        FunctionInvariant(
            "plots.plot_rank_stability",
            plot_rank_stability,
        )
    )

    job.depends_on(
        FunctionInvariant(
            "plots.plot_direction_consistency",
            plot_direction_consistency,
        )
    )

    if readout == "full":
        job.depends_on(
            FunctionInvariant(
                "plots.plot_replicate_effect_heatmap",
                plot_replicate_effect_heatmap,
            )
        )

        job.depends_on(
            FunctionInvariant(
                "plots.plot_effect_decomposition",
                plot_effect_decomposition,
            )
        )

        job.depends_on(
            FunctionInvariant(
                "plots.plot_contrast",
                plot_contrast,
            )
        )

        job.depends_on(
            FunctionInvariant(
                "plots.plot_pathway_enrichment_summary",
                plot_pathway_enrichment_summary,
            )
        )

        job.depends_on(
            FunctionInvariant(
                "plots.plot_gene_set_score_distribution",
                plot_gene_set_score_distribution,
            )
        )

        job.depends_on(
            FunctionInvariant(
                "plots.plot_beta_vs_standard_error",
                plot_beta_vs_standard_error,
            )
        )

        job.depends_on(
            FunctionInvariant(
                "plots.plot_wald_z_distribution",
                plot_wald_z_distribution,
            )
        )

        job.depends_on(
            FunctionInvariant(
                "plots.plot_qq",
                plot_qq,
            )
        )

    # Add parameter invariant
    job.depends_on(
        ParameterInvariant(
            "mageck_report_params",
            (
                str(gene_summary_path) if gene_summary_path else None,
                str(sgrna_summary_path) if sgrna_summary_path else None,
                str(rra_summary_path) if rra_summary_path else None,
                str(mle_summary_path) if mle_summary_path else None,
                readout,
                tuple(effect_cols) if effect_cols else None,
                gene_col,
                fdr_threshold,
                top_n,
                (
                    str(pathway_enrichment_path)
                    if pathway_enrichment_path
                    else None
                ),
            ),
        )
    )

    return job


# def comprehensive_qc_job(
#     count_table: Union[Path, str],
#     baseline_condition: str,
#     output_dir: Union[Path, str],
#     control_sgrnas: Union[Path, str] = None,
#     metadata_tsv: Union[Path, str] = None,
#     prefix: str = "qc",
#     project_name: str = "CRISPR Screen QC",
#     sgrna_col: str = "sgRNA",
#     gene_col: str = "Gene",
#     delimiter: str = "_",
#     norm_methods: List[str] = None,
#     run_control_qc: bool = True,
#     run_library_qc: bool = True,
#     generate_pdf: bool = False,
#     dependencies: List[Job] = [],
# ) -> MultiFileGeneratingJob:
#     """
#     Create pypipegraph job for comprehensive QC analysis.

#     Unified QC job that runs both control sgRNA QC and library-wide QC
#     using the new QCReport class from models module.

#     Parameters
#     ----------
#     count_table : Path or str
#         Path to MAGeCK count table.
#     baseline_condition : str
#         Name of baseline condition (e.g., "Total", "T0").
#     output_dir : Path or str
#         Directory for output files.
#     control_sgrnas : Path or str, optional
#         Path to control sgRNA file. If None, control QC is skipped.
#     metadata_tsv : Path or str, optional
#         Path to metadata TSV. If None, metadata is inferred from column names.
#     prefix : str
#         Prefix for output files.
#     project_name : str
#         Project name for report title.
#     sgrna_col : str
#         Name of sgRNA ID column.
#     gene_col : str
#         Name of gene column.
#     delimiter : str
#         Delimiter for parsing condition_replicate.
#     norm_methods : List[str], optional
#         Normalization methods to test. Default: ["total", "median",
#         "stable_set"]
#     run_control_qc : bool
#         Whether to run control sgRNA QC (requires control_sgrnas).
#     run_library_qc : bool
#         Whether to run library-wide QC.
#     generate_pdf : bool
#         Whether to generate PDF report (requires reportlab).
#     dependencies : List[Job]
#         Pypipegraph jobs that this job depends on.

#     Returns
#     -------
#     MultiFileGeneratingJob
#         Job that generates QC reports and plots.

#     Notes
#     -----
#     This function replaces control_qc_job and standard_qc_job with a unified
#     approach using the QCReport class. For backward compatibility, the old
#     functions remain available but are deprecated.

#     Output Files
#     ------------
#     - qc_summary.md: Markdown QC report
#     - qc_summary.json: Machine-readable QC summary
#     - qc_summary.html: HTML QC report (if markdown available)
#     - qc_report.pdf: PDF QC report (if generate_pdf=True and reportlab present)
#     - qc_assets/plots/*.png: All QC plots
#     - qc_assets/tables/*.tsv: QC data tables
#     - qc_assets/tables/*.json: QC metrics in JSON format

#     Examples
#     --------
#     >>> job = comprehensive_qc_job(
#     ...     count_table="brunello.count.tsv",
#     ...     baseline_condition="Total",
#     ...     output_dir="results/qc",
#     ...     control_sgrnas="control_sgRNAs.txt",
#     ...     project_name="My Screen",
#     ...     generate_pdf=True,
#     ... )
#     """
#     from crispr_screens.models import QCReport, QCConfig

#     output_dir = Path(output_dir)
#     count_table = Path(count_table)

#     # Determine output files
#     output_files = [
#         output_dir / "qc_summary.md",
#         output_dir / "qc_summary.json",
#     ]

#     # Add HTML if markdown available
#     if _import_util.find_spec("markdown") is not None:
#         output_files.append(output_dir / "qc_summary.html")

#     # Add PDF if requested and reportlab available
#     if generate_pdf:
#         if _import_util.find_spec("reportlab") is not None:
#             output_files.append(output_dir / "qc_report.pdf")

#     def run_qc():
#         """Execute comprehensive QC analysis."""
#         # Create QC configuration
#         config = QCConfig(
#             project_name=project_name,
#             out_dir=output_dir,
#             sgrna_col=sgrna_col,
#             gene_col=gene_col,
#             delimiter=delimiter,
#             baseline_condition=baseline_condition,
#             norm_methods=norm_methods or ["total", "median", "stable_set"],
#         )

#         # Create QC report
#         qc_report = QCReport(
#             config=config,
#             count_table_path=count_table,
#             control_sgrnas_path=control_sgrnas,
#             metadata_path=metadata_tsv,
#         )

#         # Run analysis
#         qc_report.build(
#             run_control_qc=run_control_qc and control_sgrnas is not None,
#             run_library_qc=run_library_qc,
#             generate_pdf=generate_pdf,
#         )

#     # Create job
#     job = MultiFileGeneratingJob(
#         output_files=output_files,
#         calc_function=run_qc,
#     ).depends_on(*dependencies)

#     # Add function invariants for QC functions
#     from crispr_screens.core import qc

#     job.depends_on(
#         FunctionInvariant("qc.load_control_sgrnas", qc.load_control_sgrnas)
#     )
#     job.depends_on(FunctionInvariant("qc.read_counts", qc.read_counts))
#     job.depends_on(FunctionInvariant("qc.calculate_cpm", qc.calculate_cpm))
#     job.depends_on(
#         FunctionInvariant("qc.calculate_delta_logfc", qc.calculate_delta_logfc)
#     )
#     job.depends_on(
#         FunctionInvariant(
#             "qc.qc_controls_neutrality", qc.qc_controls_neutrality
#         )
#     )
#     job.depends_on(
#         FunctionInvariant(
#             "qc.choose_best_normalization", qc.choose_best_normalization
#         )
#     )
#     job.depends_on(
#         FunctionInvariant(
#             "qc.recommend_analysis_method", qc.recommend_analysis_method
#         )
#     )

#     # Add parameter invariant
#     job.depends_on(
#         ParameterInvariant(
#             "comprehensive_qc_params",
#             (
#                 str(count_table),
#                 baseline_condition,
#                 str(output_dir),
#                 str(control_sgrnas) if control_sgrnas else None,
#                 str(metadata_tsv) if metadata_tsv else None,
#                 project_name,
#                 tuple(norm_methods) if norm_methods else None,
#                 run_control_qc,
#                 run_library_qc,
#                 generate_pdf,
#             ),
#         )
#     )

#     return job
