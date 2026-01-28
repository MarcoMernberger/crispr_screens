import pandas as pd  # type: ignore
from typing import Dict, Optional, Union, List, Tuple, Callable
from pathlib import Path
from crispr_screens.core.mageck import (
    filter_multiple_mageck_comparisons,
    combine_comparisons,
    split_frame_to_control_and_query,
    combine_gene_info_with_mageck_output,
)
from crispr_screens.core.mageck_spikein import create_spiked_count_table


def write_filtered_mageck_comparison(
    output_file: Path,
    combined_frame_input_file: Path,
    comparisons_to_filter: List[str],
    fdr_threshold: Union[float, Dict[str, float]] = 0.05,
    change_threshold: Union[float, Dict[str, float]] = 1.0,
    z_thresholds: Optional[Union[float, Dict[str, float]]] = None,
    direction: str = "both",  # "both", "pos", "neg"
    require_all: bool = True,  # AND (True) vs OR (False)
):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined_frame = pd.read_csv(combined_frame_input_file, sep="\t")
    filtered_frame = filter_multiple_mageck_comparisons(
        combined_frame=combined_frame,
        comparisons_to_filter=comparisons_to_filter,
        fdr_threshold=fdr_threshold,
        change_threshold=change_threshold,
        z_thresholds=z_thresholds,
        direction=direction,
        require_all=require_all,
    )
    filtered_frame.to_csv(output_file, sep="\t", index=False)


def combine_comparison_output(
    output_file: Path,
    mageck_results: Dict[str, Path],
    combine_on: Union[str, Dict[str, str]] = "id",
    how: str = "inner",
):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    mageck_frames = {}
    for name, mageck_file in mageck_results.items():
        mageck_frames[name] = pd.read_csv(mageck_file, sep="\t")
    merged_frame = combine_comparisons(
        mageck_frames, combine_on=combine_on, how=how
    )
    merged_frame.to_csv(output_file, sep="\t", index=False)


def create_query_control_sgrna_frames(
    infile: Path,
    outfiles: Tuple[Path],
    control_prefix: str,
    id_col: Optional[str] = None,
    name_column: str = "name",
    sgRNA_column: str = "sgRNA",
    infer_genes: Optional[Callable] = None,
    read_csv_kwargs: Optional[Dict] = None,
):
    outfiles[0].parent.mkdir(parents=True, exist_ok=True)
    read_csv_kwargs = read_csv_kwargs or {"sep": "\t"}
    df = pd.read_csv(infile, **read_csv_kwargs)
    split_dfs = split_frame_to_control_and_query(
        df,
        control_prefix=control_prefix,
        id_col=id_col,
        name_column=name_column,
        sgRNA_column=sgRNA_column,
        infer_genes=infer_genes,
    )
    split_dfs["control"].to_csv(
        outfiles[1], sep="\t", index=False, header=False
    )
    split_dfs["query"].to_csv(outfiles[0], sep="\t", index=False, header=False)


def create_combine_gene_info_with_mageck_output(
    mageck_file: Path,
    gene_info_file: Path,
    output_file: Path,
    name_column_mageck: str = "id",
    name_column_genes: str = "name_given",
    how: str = "inner",
    columns_to_add: List[str] = [
        "gene_stable_id",
        "name",
        "chr",
        "start",
        "stop",
        "strand",
        "tss",
        "tes",
        "biotype",
    ],
    read_csv_kwargs: Optional[Dict] = None,
):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    read_csv_kwargs = read_csv_kwargs or {"sep": "\t"}
    mageck_df = pd.read_csv(mageck_file, **read_csv_kwargs)
    gene_info_df = pd.read_csv(gene_info_file, **read_csv_kwargs)
    combined_df = combine_gene_info_with_mageck_output(
        mageck_df,
        gene_info_df,
        name_column_mageck=name_column_mageck,
        name_column_genes=name_column_genes,
        how=how,
        columns_to_add=columns_to_add,
    )
    combined_df.to_csv(output_file, sep="\t", index=False)


def write_spiked_count_table(
    outfile: Union[Path, str],
    count_file: Union[str, Path],
    replicate_of: Dict[str, str],
    sample_to_group: Dict[str, str],
    group_contrast: Tuple[str] = ("sorted", "total"),
    n_genes: int = 20,
    log_effect: float = 2.0,
    baseline_mean: float = 300.0,
    dispersion: float = 0.08,
):
    outfile.parent.mkdir(parents=True, exist_ok=True)
    count_df = pd.read_csv(count_file, sep="\t")
    sample_cols = [c for c in count_df.columns if (c not in ("sgRNA", "Gene"))]
    count_spike = create_spiked_count_table(
        count_df,
        replicate_of=replicate_of,
        sample_to_group=sample_to_group,
        sample_cols=sample_cols,
        group_contrast=group_contrast,
        n_genes=n_genes,
        log_effect=log_effect,
        baseline_mean=baseline_mean,
        dispersion=dispersion,
    )
    count_spike.to_csv(outfile, sep="\t", index=False)
