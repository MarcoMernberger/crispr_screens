from pypipegraph2 import Job, FileGeneratingJob, MultiFileGeneratingJob
from pathlib import Path
from typing import Dict, Optional, Union, List, Tuple, Literal
from crispr_screens.services.mageck_io import (
    write_filtered_mageck_comparison,
    combine_comparison_output,
    write_venn,
)
from crispr_screens.r_integration.mageck_wrapper import run_mageck_scatterview
from crispr_screens.core.mageck import mageck_count, mageck_test, mageck_mle


def combine_comparison_output_job(
    mageck_results: Dict[str, Path],
    output_file: Path,
    combine_on: Union[str, Dict[str, str]] = "id",
    how: str = "inner",
    dependencies: List[Job] = [],
):
    def __dump(
        output_file,
        mageck_results=mageck_results,
        combine_on=combine_on,
        how=how,
    ):
        combine_comparison_output(
            output_file=output_file,
            mageck_results=mageck_results,
            combine_on=combine_on,
            how=how,
        )

    return FileGeneratingJob(output_file, __dump).depends_on(dependencies)


def write_filtered_mageck_comparison_job(
    output_file: Path,
    combined_frame_input_file: Path,
    comparisons_to_filter: List[str],
    fdr_threshold: Union[float, Dict[str, float]] = 0.05,
    change_threshold: Union[float, Dict[str, float]] = 1.0,
    z_thresholds: Optional[Union[float, Dict[str, float]]] = None,
    direction: str = "both",  # "both", "pos", "neg"
    require_all: bool = True,  # AND (True) vs OR (False)
    dependencies: List[Job] = [],
):

    def __dump(
        output_file,
        combined_frame_input_file=combined_frame_input_file,
        comparisons_to_filter=comparisons_to_filter,
        fdr_threshold=fdr_threshold,
        change_threshold=change_threshold,
        z_thresholds=z_thresholds,
        direction=direction,
        require_all=require_all,
    ):
        write_filtered_mageck_comparison(
            output_file=output_file,
            combined_frame_input_file=combined_frame_input_file,
            comparisons_to_filter=comparisons_to_filter,
            fdr_threshold=fdr_threshold,
            change_threshold=change_threshold,
            z_thresholds=z_thresholds,
            direction=direction,
            require_all=require_all,
        )

    return FileGeneratingJob(output_file, __dump).depends_on(dependencies)


def run_mageck_scatterview_job(
    input_file: Union[str, Path],
    x_col: str,
    y_col: str,
    output_dir: Union[str, Path] = ".",
    filebase_name: str = "scatterview",
    datatype: Literal["rra", "mle"] = "rra",
    gene_col: str = "id",
    sep: str = "\t",
    normalize: bool = True,
    top: int = 10,
    groups: Tuple[str] = (
        "bottomleft",
        "topcenter",
        "midright",
    ),  # selection only those groups
    select: Optional[
        Literal["positive", "negative", "both", "none"]
    ] = "negative",  # selection by diagonal None, "positive", "negative", "both", "none"
    neg_effect_cutoff: float = -0.4,
    pos_effect_cutoff: float = 0.4,
    delta_cutoff_k: float = 2,
    filter_fdr_x: bool = False,  # selection filter by FDR column
    filter_fdr_y: bool = False,  # selection filter by FDR column
    filter_groups: bool = True,  # if false groups will be displayed but not filtered
    fdr_cutoff: float = 0.05,  # FDR cutoff
    # plot parameters
    toplabels: bool = True,
    label_selected_only: bool = False,
    xlab: Optional[str] = None,
    ylab: Optional[str] = None,
    jpeg_width: int = 20,
    jpeg_height: int = 15,
    # additional plot parameter for MLE
    auto_cut_diag: float = 2,
    auto_cut_x: float = 2,
    auto_cut_y: float = 2,
    # additional plot parameter for RRA
    x_cut: Optional[Tuple[float, float]] = None,  # e.g. (-0.5, 0.5),
    y_cut: Optional[Tuple[float, float]] = None,  # e.g. (-0.5, 0.5),
    dependencies: List[Job] = [],
):
    outfiles = [
        Path(output_dir) / f"{filebase_name}_data.tsv",
        Path(output_dir) / f"{filebase_name}_hits_selected.tsv",
        Path(output_dir) / f"{filebase_name}.jpeg",
        Path(output_dir) / f"{filebase_name}_lmfit.jpeg",
    ]

    def __dump(
        outfiles,
        input_file=input_file,
        x_col=x_col,
        y_col=y_col,
        output_dir=output_dir,
        filebase_name=filebase_name,
        datatype=datatype,
        gene_col=gene_col,
        sep=sep,
        normalize=normalize,
        top=top,
        groups=groups,
        select=select,
        neg_effect_cutoff=neg_effect_cutoff,
        pos_effect_cutoff=pos_effect_cutoff,
        delta_cutoff_k=delta_cutoff_k,
        filter_fdr_x=filter_fdr_x,
        filter_fdr_y=filter_fdr_y,
        filter_groups=filter_groups,
        fdr_cutoff=fdr_cutoff,
        toplabels=toplabels,
        label_selected_only=label_selected_only,
        xlab=xlab,
        ylab=ylab,
        jpeg_width=jpeg_width,
        jpeg_height=jpeg_height,
        auto_cut_diag=auto_cut_diag,
        auto_cut_x=auto_cut_x,
        auto_cut_y=auto_cut_y,
        x_cut=x_cut,
        y_cut=y_cut,
    ):
        run_mageck_scatterview(
            input_file=input_file,
            x_col=x_col,
            y_col=y_col,
            output_dir=output_dir,
            filebase_name=filebase_name,
            datatype=datatype,
            gene_col=gene_col,
            sep=sep,
            normalize=normalize,
            top=top,
            groups=groups,
            select=select,
            neg_effect_cutoff=neg_effect_cutoff,
            pos_effect_cutoff=pos_effect_cutoff,
            delta_cutoff_k=delta_cutoff_k,
            filter_fdr_x=filter_fdr_x,
            filter_fdr_y=filter_fdr_y,
            filter_groups=filter_groups,
            fdr_cutoff=fdr_cutoff,
            toplabels=toplabels,
            label_selected_only=label_selected_only,
            xlab=xlab,
            ylab=ylab,
            jpeg_width=jpeg_width,
            jpeg_height=jpeg_height,
            auto_cut_diag=auto_cut_diag,
            auto_cut_x=auto_cut_x,
            auto_cut_y=auto_cut_y,
            x_cut=x_cut,
            y_cut=y_cut,
        )

    return MultiFileGeneratingJob(outfiles, __dump).depends_on(dependencies)


def write_venn_job(
    outdir: Union[Path, str],
    filebasename: str,
    label_to_file: Dict[str, str],
    id_cols: Union[List[str], str] = "id",
    sep: str = "\t",
    figsize: Tuple[float, float] = (6, 6),
    title: str | None = None,
    dependencies: List[Job] = [],
):
    out_venn = outdir / f"{filebasename}_venn.png"
    out_dataframe = outdir / f"{filebasename}.tsv"

    def __dump(
        outfiles,
        outdir=outdir,
        filebasename=filebasename,
        label_to_file=label_to_file,
        id_cols=id_cols,
        sep=sep,
        figsize=figsize,
        title=title,
    ):
        write_venn(
            outdir=outdir,
            filebasename=filebasename,
            label_to_file=label_to_file,
            id_cols=id_cols,
            sep=sep,
            figsize=figsize,
            title=title,
        )

    return MultiFileGeneratingJob([out_venn, out_dataframe], __dump).depends_on(
        dependencies
    )


def mageck_count_job(
    sgrna_list: Union[Path, str],
    samples: dict,
    out_dir: Union[Path, str],
    prefix: str,
    control_sgrnas=Optional[: Union[Path, str]],
    norm_method: str = None,
    pdf_report: bool = False,
    other_parameter: Dict[str, str] = [],
    dependencies: List[Job] = [],
):
    outfile = Path(out_dir) / f"{prefix}.count.txt"

    def __dump(
        outfile,
        sgrna_list=sgrna_list,
        samples=samples,
        out_dir=out_dir,
        prefix=prefix,
        control_sgrnas=control_sgrnas,
        norm_method=norm_method,
        pdf_report=pdf_report,
        other_parameter=other_parameter,
    ):
        mageck_count(
            sgrna_list=sgrna_list,
            samples=samples,
            out_dir=out_dir,
            prefix=prefix,
            control_sgrnas=control_sgrnas,
            norm_method=norm_method,
            pdf_report=pdf_report,
            other_parameter=other_parameter,
        )

    return FileGeneratingJob(outfile, __dump).depends_on(dependencies)


def mageck_rra_job(
    count_table: Union[Path, str],
    treatment_ids: List[str],
    control_ids: List[str],
    out_dir: Union[Path, str],
    prefix: str,
    control_sgrnas: Optional[Union[Path, str]],
    norm_method: str = None,
    paired: bool = False,
    pdf_report: bool = False,
    other_parameter: Dict[str, str] = [],
    dependencies: List[Job] = [],
):
    outfile = Path(out_dir) / f"{prefix}.gene_summary.tsv"

    def __dump(
        outfile,
        count_table=count_table,
        treatment_ids=treatment_ids,
        control_ids=control_ids,
        out_dir=out_dir,
        prefix=prefix,
        control_sgrnas=control_sgrnas,
        norm_method=norm_method,
        paired=paired,
        pdf_report=pdf_report,
        other_parameter=other_parameter,
    ):
        mageck_test(
            count_table=count_table,
            treatment_ids=treatment_ids,
            control_ids=control_ids,
            out_dir=out_dir,
            prefix=prefix,
            control_sgrnas=control_sgrnas,
            norm_method=norm_method,
            paired=paired,
            pdf_report=pdf_report,
            other_parameter=other_parameter,
        )

    return FileGeneratingJob(outfile, __dump).depends_on(dependencies)


def mageck_mle_job(
    count_table: Union[Path, str],
    design_matrix: str,
    out_dir: Union[Path, str],
    prefix: str,
    control_sgrnas: Optional[Union[Path, str]] = None,
    norm_method: str = None,
    other_parameter: Dict[str, str] = [],
    dependencies: List[Job] = [],
):
    outfile = Path(out_dir) / f"{prefix}.gene_summary.tsv"

    def __dump(
        outfile,
        count_table=count_table,
        design_matrix=design_matrix,
        out_dir=out_dir,
        prefix=prefix,
        control_sgrnas=control_sgrnas,
        norm_method=norm_method,
        other_parameter=other_parameter,
    ):
        mageck_mle(
            count_table=count_table,
            design_matrix=design_matrix,
            out_dir=out_dir,
            prefix=prefix,
            control_sgrnas=control_sgrnas,
            norm_method=norm_method,
            other_parameter=other_parameter,
        )

    return FileGeneratingJob(outfile, __dump).depends_on(dependencies)
