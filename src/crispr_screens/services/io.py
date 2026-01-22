import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Union, List, Dict, Optional
from pandas import DataFrame

from crispr_screens.core.qc import (
    control_sgrna_qc,
    export_control_counts_and_cpm,
    generate_standard_qc_report as _generate_standard_qc_report,
)
from crispr_screens.core.mageck_report import (
    generate_mageck_report as _generate_mageck_report,
)
from crispr_screens.core.plots import (
    plot_control_distribution_per_condition,
    plot_control_pca,
    plot_control_replicate_correlation,
    plot_pairwise_control_shifts,
)
from crispr_screens.core.mageck import (
    mageck_pathway as _mageck_pathway,
    mageck_plot as _mageck_plot,
)


def save_figure(f, folder, name, bbox_inches="tight"):
    folder.mkdir(exist_ok=True, parents=True)
    for suffix in [".png", ".svg", ".pdf"]:
        f.savefig(folder / (name + suffix), bbox_inches=bbox_inches)


# def control_qc_report(
#     count_table: Union[Path, str],
#     control_sgrnas: Union[Path, str],
#     baseline_condition: str,
#     output_dir: Union[Path, str],
#     prefix: str = "control_qc",
#     sgrna_col: str = "sgRNA",
#     gene_col: str = "Gene",
#     delimiter: str = "_",
#     save_formats: List[str] = ["png", "pdf"],
# ):
#     """
#     Generate comprehensive control sgRNA QC report.

#     Wrapper around core.qc.generate_control_qc_report for service layer.

#     Parameters
#     ----------
#     count_table : Path or str
#         Path to MAGeCK count table.
#     control_sgrnas : Path or str
#         Path to control sgRNA file.
#     baseline_condition : str
#         Baseline condition name (e.g., "Total", "T0").
#     output_dir : Path or str
#         Output directory for QC reports.
#     prefix : str
#         Filename prefix for output files.
#     sgrna_col : str
#         sgRNA column name in count table.
#     gene_col : str
#         Gene column name in count table.
#     delimiter : str
#         Delimiter for parsing condition_replicate.
#     save_formats : list
#         List of formats to save ("png", "pdf", "svg").

#     Returns
#     -------
#     dict
#         QC results with file paths.
#     """
#     from crispr_screens.core.qc import generate_control_qc_report

#     return generate_control_qc_report(
#         count_table=count_table,
#         control_sgrnas=control_sgrnas,
#         baseline_condition=baseline_condition,
#         output_dir=output_dir,
#         prefix=prefix,
#         sgrna_col=sgrna_col,
#         gene_col=gene_col,
#         delimiter=delimiter,
#         save_formats=save_formats,
#     )


def generate_control_qc_report(
    count_table: Union[Path, str, DataFrame],
    control_sgrnas: Union[Path, str, set],
    baseline_condition: str,
    output_dir: Union[Path, str],
    prefix: str = "control_qc",
    sgrna_col: str = "sgRNA",
    gene_col: str = "Gene",
    delimiter: str = "_",
    save_formats: List[str] = ["png", "pdf"],
) -> Dict:
    """
    Generate comprehensive control sgRNA QC report.

    Creates all QC plots and saves metrics to file.

    Parameters
    ----------
    count_table : Path, str, or DataFrame
        MAGeCK count table.
    control_sgrnas : Path, str, or set
        Control sgRNA IDs.
    baseline_condition : str
        Baseline condition name.
    output_dir : Path or str
        Output directory.
    prefix : str
        Filename prefix.
    sgrna_col : str
        sgRNA column name.
    gene_col : str
        Gene column name.
    delimiter : str
        Condition/replicate delimiter.
    save_formats : list
        List of formats to save ("png", "pdf", "svg").

    Returns
    -------
    dict
        QC results with file paths.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Run QC analysis
    print("Running control sgRNA QC analysis...")
    qc_results = control_sgrna_qc(
        count_table=count_table,
        control_sgrnas=control_sgrnas,
        baseline_condition=baseline_condition,
        sgrna_col=sgrna_col,
        gene_col=gene_col,
        delimiter=delimiter,
    )

    # Save metrics
    metrics_file = output_dir / f"{prefix}_metrics.tsv"
    metrics_df = pd.DataFrame(qc_results["metrics"]).T
    metrics_df.to_csv(metrics_file, sep="\t")
    print(f"Saved metrics to {metrics_file}")

    # Save pairwise median shifts
    pairwise_file = output_dir / f"{prefix}_pairwise_shifts.tsv"
    qc_results["pairwise_median"].to_csv(pairwise_file, sep="\t")
    print(f"Saved pairwise shifts to {pairwise_file}")
    cpm_files = export_control_counts_and_cpm(
        count_table=qc_results["raw_counts"],
        control_sgrnas=control_sgrnas,
        output_dir=output_dir,
        prefix=f"{prefix}",
        sgrna_col=sgrna_col,
        gene_col=gene_col,
    )
    saved_files = {
        "metrics": metrics_file,
        "pairwise_shifts": pairwise_file,
    }
    saved_files.update(cpm_files)

    # Generate and save plots
    plots_to_generate = [
        ("distribution", plot_control_distribution_per_condition),
        ("pairwise_heatmap", plot_pairwise_control_shifts),
        ("replicate_correlation", plot_control_replicate_correlation),
        (
            "pca_condition",
            lambda qc: plot_control_pca(qc, color_by="condition"),
        ),
        (
            "pca_replicate",
            lambda qc: plot_control_pca(qc, color_by="replicate"),
        ),
    ]

    for plot_name, plot_func in plots_to_generate:
        print(f"Generating {plot_name} plot...")
        try:
            fig_result = plot_func(qc_results)
            if isinstance(fig_result, tuple):
                fig = fig_result[0]
            else:
                fig = fig_result

            for fmt in save_formats:
                outfile = output_dir / f"{prefix}_{plot_name}.{fmt}"
                fig.savefig(outfile, dpi=300, bbox_inches="tight")
                saved_files[f"{plot_name}_{fmt}"] = outfile

            plt.close(fig)
            print(f"  Saved {plot_name}")
        except Exception as e:
            print(f"  Warning: Failed to generate {plot_name}: {e}")

    print(f"\nControl QC report complete. Files saved to {output_dir}")

    return {
        "qc_results": qc_results,
        "files": saved_files,
    }


def mageck_pathway(
    gene_ranking: Union[Path, str],
    gmt_file: Union[Path, str],
    output_dir: Union[Path, str],
    prefix: str = "pathway",
    method: str = "gsea",
    single_ranking: bool = False,
    sort_criteria: str = "neg",
    keep_tmp: bool = False,
    ranking_column: Optional[Union[str, int]] = None,
    ranking_column_2: Optional[Union[str, int]] = None,
    pathway_alpha: Optional[float] = None,
    permutation: Optional[int] = None,
    other_parameter: Dict[str, str] = [],
) -> Dict:
    """
    Service wrapper for `mageck pathway`.

    Returns a dict from the underlying core function containing stdout/stderr
    and output files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    return _mageck_pathway(
        gene_ranking=gene_ranking,
        gmt_file=gmt_file,
        out_dir=output_dir,
        prefix=prefix,
        method=method,
        single_ranking=single_ranking,
        sort_criteria=sort_criteria,
        keep_tmp=keep_tmp,
        ranking_column=ranking_column,
        ranking_column_2=ranking_column_2,
        pathway_alpha=pathway_alpha,
        permutation=permutation,
        other_parameter=other_parameter,
    )


def mageck_plot(
    gene_summary: Optional[Union[Path, str]] = None,
    sgrna_summary: Optional[Union[Path, str]] = None,
    output_dir: Union[Path, str] = ".",
    prefix: str = "plot",
    other_parameter: Dict[str, str] = [],
) -> Dict:
    """
    Service wrapper for `mageck plot`.

    Returns a dict from the underlying core function containing stdout/stderr
    and output files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    return _mageck_plot(
        gene_summary=gene_summary,
        sgrna_summary=sgrna_summary,
        out_dir=output_dir,
        prefix=prefix,
        other_parameter=other_parameter,
    )


def read_dataframe(path: Union[str, Path], **kwargs) -> DataFrame:
    """
    Read a tabular file into a pandas DataFrame based on file extension.

    Rules:
    - .csv        -> read as CSV
    - .tsv        -> read as TSV
    - .txt        -> treated as TSV
    - .xls/.xlsx  -> read as Excel
    - other       -> try TSV, raise error if that fails

    Additional keyword arguments are forwarded to the pandas reader.
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"File does not exist: {path}")

    suffix = path.suffix.lower()

    try:
        if suffix == ".csv":
            return pd.read_csv(path, **kwargs)

        if suffix in {".tsv", ".txt"}:
            return pd.read_csv(path, sep="\t", **kwargs)

        if suffix in {".xls", ".xlsx"}:
            return pd.read_excel(path, **kwargs)

        # Fallback: try TSV for unknown extensions
        try:
            return pd.read_csv(path, sep="\t", **kwargs)
        except Exception as exc:
            raise ValueError(
                f"Unsupported file extension '{suffix}'. "
                "Tried to read as TSV but failed."
            ) from exc

    except Exception as exc:
        raise RuntimeError(f"Failed to read file '{path}': {exc}") from exc


def standard_qc_report(
    count_tsv: Union[Path, str],
    output_dir: Union[Path, str],
    metadata_tsv: Union[Path, str, None] = None,
    control_sgrna_txt: Union[Path, str, None] = None,
    baseline_condition: str = "total",
    **kwargs,
) -> Dict:
    """
    Generate standard QC report for CRISPR screen.

    Service wrapper for core.qc.generate_standard_qc_report.

    Parameters
    ----------
    count_tsv : Path or str
        Path to MAGeCK count table.
    output_dir : Path or str
        Output directory.
    metadata_tsv : Path or str, optional
        Path to metadata file.
    control_sgrna_txt : Path or str, optional
        Path to control sgRNA file.
    baseline_condition : str
        Baseline condition name.
    **kwargs
        Additional arguments passed to generate_standard_qc_report.

    Returns
    -------
    dict
        QC results.
    """
    return _generate_standard_qc_report(
        count_tsv=count_tsv,
        output_dir=output_dir,
        metadata_tsv=metadata_tsv,
        control_sgrna_txt=control_sgrna_txt,
        baseline_condition=baseline_condition,
        **kwargs,
    )


def mageck_report(
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
    **kwargs,
):
    """
    Generate comprehensive MAGeCK reporting with specialized plots.

    Service layer wrapper around core.mageck_report.generate_mageck_report.

    Creates publication-ready visualizations addressing key questions:
    - "Is this a real hit?" → Effect-size vs reproducibility, rank stability
    - "Is this experimentally stable?" → Replicate heatmaps, direction
       consistency
    - "What is functionally connected?" → Pathway enrichment, gene-set
      distributions
    - "How good is the model?" → Beta vs SE, Wald-Z, QQ-plots

    Parameters
    ----------
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
        Report mode: "minimal" (5 key plots) or "full" (all 11+ plots)
        - minimal: volcano, effect-vs-reproducibility, rank stability,
                   direction consistency, pathway enrichment
        - full: adds replicate heatmap, effect decomposition, contrast,
                gene-set distributions, beta-vs-SE, Wald-Z, QQ-plot,
                plus executive summary metrics
    effect_cols : list, optional
        For MLE: list of beta column names for effect decomposition
        E.g., ['beta|Time', 'beta|Met', 'beta|High']
    gene_col : str
        Column name for gene identifiers (default "Gene")
    fdr_threshold : float
        FDR threshold for significance (default 0.25)
    top_n : int
        Number of top genes to highlight in plots (default 50)
    gene_sets : dict, optional
        Dictionary mapping gene set names to lists of genes
        For gene-set score distribution analysis
        Example: {'Kinases': ['BRAF', 'RAF1', ...], 'TFs': ['TP53', ...]}
    pathway_enrichment_path : str or Path, optional
        Path to pathway enrichment results (CSV/TSV)
        Expected columns: pathway, pvalue, gene_ratio
    **kwargs : dict
        Additional parameters passed to plotting functions

    Returns
    -------
    dict
        Summary metrics and paths to generated plots containing:
        - 'output_dir': Path to output directory
        - 'readout': Mode used ('minimal' or 'full')
        - 'plots': Dictionary mapping plot types to file paths
        - 'metrics': Dictionary of computed metrics per plot
        - 'executive_summary': Overall summary metrics (full mode only)

    Examples
    --------
    Minimal mode (quick check):
    >>> results = mageck_report(
    ...     output_dir="results/mageck_report",
    ...     gene_summary_path="results/mageck_mle/gene_summary.txt",
    ...     sgrna_summary_path="results/mageck_mle/sgrna_summary.txt",
    ...     readout="minimal"
    ... )

    Full mode with method comparison:
    >>> results = mageck_report(
    ...     output_dir="results/mageck_report_full",
    ...     rra_summary_path="results/mageck_test/gene_summary.txt",
    ...     mle_summary_path="results/mageck_mle/gene_summary.txt",
    ...     sgrna_summary_path="results/mageck_mle/sgrna_summary.txt",
    ...     readout="full",
    ...     effect_cols=['beta|Treatment', 'beta|Time'],
    ...     gene_sets={'Kinases': kinase_genes, 'TFs': tf_genes},
    ...     fdr_threshold=0.25
    ... )
    """
    return _generate_mageck_report(
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
