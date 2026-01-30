import pandas as pd  # type: ignore
import re
import shutil
import subprocess
import rpy2.robjects as ro
from typing import Dict, Optional, Union, List, Callable
from pandas import DataFrame  # type: ignore
from pathlib import Path


def combine_comparisons(
    mageck_results: Dict[str, DataFrame],
    combine_on: Union[str, Dict[str, str]] = "id",
    how: str = "inner",
) -> DataFrame:
    keys = list(mageck_results.keys())
    combine_on_columns = combine_on
    if isinstance(combine_on, str):
        combine_on_columns = dict([(key, combine_on) for key in keys])
    for key in mageck_results:
        mageck_results[key] = mageck_results[key].rename(
            columns={
                c: f"{key}|{c}"
                for c in mageck_results[key].columns
                if c != combine_on_columns[key]
            }
        )
    combine_on_columns = combine_on
    if isinstance(combine_on, str):
        combine_on_columns = dict([(key, combine_on) for key in keys])
    merged_frames = mageck_results[keys[0]].merge(
        mageck_results[keys[1]],
        left_on=combine_on_columns[keys[0]],
        right_on=combine_on_columns[keys[1]],
    )
    if len(keys) > 2:
        for key2 in keys[2:]:
            print(key2)
            merged_frames = merged_frames.merge(
                mageck_results[key2],
                left_on=combine_on_columns[keys[0]],
                right_on=combine_on_columns[key2],
            )
    return merged_frames


def filter_multiple_mageck_comparisons(
    combined_frame: DataFrame,
    comparisons_to_filter: List[str],
    fdr_threshold: Union[float, Dict[str, float]] = 0.05,
    change_threshold: Union[float, Dict[str, float]] = 1.0,
    z_thresholds: Optional[Union[float, Dict[str, float]]] = None,
    direction: str = "both",  # "both", "pos", "neg"
    require_all: bool = True,  # AND (True) vs OR (False)
) -> DataFrame:
    """
    Filter combined MAGeCK RRA and/or MLE output frames.

    This function automatically detects whether a comparison column comes
    from MAGeCK-MLE (beta, z, wald-fdr) or MAGeCK-RRA (lfc, fdr).

    Expected patterns (before suffixing by combine_comparisons):
        MLE:
            <condition>|beta
            <condition>|z
            <condition>|fdr
            <condition>|wald-fdr

        RRA:
            neg|lfc, pos|lfc
            neg|fdr, pos|fdr

    Suffixing:
        When merging, columns become for example:
            neg|lfc_mycomparison
            high_T21|beta_othercomparison

    The function automatically detects whether a comparison is MLE-like
    (beta, z, wald-fdr) or RRA-like (lfc, fdr) based on column names.
    """

    df = combined_frame

    def _to_dict(val, keys):
        """convert scalar thresholds to dicts keyed by comparison names"""
        if isinstance(val, dict):
            return val
        return {k: val for k in keys}

    def find_matching(pattern_list):
        """find columns that match any pattern + the comparison suffix"""
        regexes = [re.compile(key + p + r"$") for p in pattern_list]
        matches = []
        for col in df.columns:
            if any(r.search(col) for r in regexes):
                matches.append(col)
        return matches

    fdr_thr = _to_dict(fdr_threshold, comparisons_to_filter)
    change_thr = _to_dict(change_threshold, comparisons_to_filter)
    z_thr = (
        _to_dict(z_thresholds, comparisons_to_filter)
        if z_thresholds is not None
        else None
    )

    # Initialize global mask:
    if require_all:
        # require_all=True → logical AND
        global_mask = pd.Series(True, index=df.index)
    else:
        # require_all=False → logical OR
        global_mask = pd.Series(False, index=df.index)

    # regex patterns help us identify relevant columns, even after suffixing
    fdr_patterns = [
        r"\|fdr",  # RRA and MLE fdr after merging
        r"\|wald-fdr",  # MLE wald-fdr
    ]  # FDR columns (MLE and RRA)
    eff_patterns = [
        r"\|beta",  # MLE effect size
        r"\|lfc",  # RRA effect size
    ]  # effect size columns
    z_patterns = [r"\|z_"]  # Z score columns (MLE only)
    for key in comparisons_to_filter:
        # loop over comparisons
        eff_cols = find_matching(eff_patterns)
        fdr_cols = find_matching(fdr_patterns)
        z_cols = find_matching(z_patterns)
        if not fdr_cols or not eff_cols:
            # no valid columns for this comparison
            raise ValueError(f"No fdr or effect column found for key={key}.")

        if len(fdr_cols) > 1:
            fdr_cols = [x for x in fdr_cols if "wald-fdr" in x]
        fdr_col = fdr_cols[0]
        eff_col = eff_cols[0]
        z_col = z_cols[0] if z_cols else None

        # build filter mask for this comparison
        mask = pd.Series(True, index=df.index)
        if fdr_thr[key] is not None:
            mask = df[fdr_col] <= fdr_thr[key]  # FDR filter
        if direction == "both":
            mask &= df[eff_col].abs() >= change_thr[key]  # effect size filter
        elif direction == "pos":
            mask &= df[eff_col] >= change_thr[key]
        elif direction == "neg":
            mask &= df[eff_col] <= -change_thr[key]

        else:
            raise ValueError(
                "direction must be one of: 'both', 'pos', or 'neg'"
            )
        df["test"] = df[eff_col] >= change_thr[key]
        if z_col is not None and z_thr is not None and key in z_thr:
            mask &= df[z_col].abs() >= z_thr[key]  # Z-filter (MLE only)

        # Combine mask depending on AND/OR logic
        if require_all:
            global_mask &= mask
        else:
            global_mask |= mask
    return df[global_mask].copy()


def mageck_count(
    sgrna_list: Union[Path, str],
    samples: dict,
    out_dir: Union[Path, str],
    prefix: str,
    control_sgrnas=Optional[: Union[Path, str]],
    norm_method: str = "median",
    pdf_report: bool = False,
    other_parameter: Dict[str, str] = [],
):

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    if isinstance(list(samples.values())[1], list):
        fastqs = " ".join(f"{",".join(fastq)}" for fastq in samples.values())
    else:
        fastqs = " ".join(fastq for fastq in samples.values())

    sample_labels = ",".join(sample_name for sample_name in samples.keys())

    command_parameters = [
        "-l",
        sgrna_list,
        "--fastq",
        fastqs,
        "--sample-label",
        sample_labels,
        "-n",
        f"{out_dir}/{prefix}",
    ]

    if control_sgrnas is not None:
        command_parameters.extend(["--control-sgrna", str(control_sgrnas)])

    if norm_method is not None:
        command_parameters.extend(["--norm-method", str(norm_method)])
    else:
        command_parameters.extend(["--norm-method", "none"])
    if pdf_report:
        command_parameters.append("--pdf-report")

    command_parameters.extend(other_parameter)

    command = ["mageck count"]
    command.extend(command_parameters)
    cmd = " ".join(command)

    print(cmd)

    mageck_count = subprocess.run(
        cmd, capture_output=True, text=True, shell=True
    )
    print(mageck_count.stdout)
    print(mageck_count.stderr)

    # if pdf_report:
    ro.r("rmarkdown::render")(f"{out_dir}/{prefix}.count_report.Rmd")

    count_txt = Path(f"{out_dir}/{prefix}.count.txt")
    count_tsv = Path(f"{out_dir}/{prefix}.count.tsv")
    if count_txt.is_file():
        shutil.copy(count_txt, count_tsv)


def mageck_test(
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
):

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    command_parameters = [
        "-k",
        str(count_table),
        "-t",
        ",".join(treatment_ids),
        "-c",
        ",".join(control_ids),
        "-n",
        f"{out_dir}/{prefix}",
    ]

    if control_sgrnas is not None:
        command_parameters.extend(["--control-sgrna", str(control_sgrnas)])

    if norm_method is not None:
        command_parameters.extend(["--norm-method", str(norm_method)])

    if paired:
        command_parameters.append("--paired")

    if pdf_report:
        command_parameters.append("--pdf-report")
    for k in other_parameter:
        command_parameters.extend([k, other_parameter[k]])

    command = ["mageck test"]
    command.extend(command_parameters)
    cmd = " ".join(command)
    print(cmd)

    mageck_count = subprocess.run(
        cmd, capture_output=True, text=True, shell=True
    )
    print(mageck_count.stdout)
    print(mageck_count.stderr)

    # if pdf_report:
    ro.r("rmarkdown::render")(f"{out_dir}/{prefix}.report.Rmd")

    gene_summary_txt = Path(f"{out_dir}/{prefix}.gene_summary.txt")
    gene_summary_tsv = Path(f"{out_dir}/{prefix}.gene_summary.tsv")

    sgrna_summary_txt = Path(f"{out_dir}/{prefix}.sgrna_summary.txt")
    sgrna_summary_tsv = Path(f"{out_dir}/{prefix}.sgrna_summary.tsv")

    if gene_summary_txt.is_file():
        shutil.copy(gene_summary_txt, gene_summary_tsv)
    if sgrna_summary_txt.is_file():
        shutil.copy(sgrna_summary_txt, sgrna_summary_tsv)


def mageck_mle(
    count_table: Union[Path, str],
    design_matrix: str,
    out_dir: Union[Path, str],
    prefix: str,
    control_sgrnas: Optional[Union[Path, str]] = None,
    norm_method: str = None,
    other_parameter: Dict[str, str] = [],
):

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    command_parameters = [
        "-k",
        count_table,
        "-d",
        design_matrix,
        "-n",
        f"{out_dir}/{prefix}",
    ]

    if control_sgrnas is not None:
        command_parameters.extend(["--control-sgrna", str(control_sgrnas)])

    if norm_method is not None:
        command_parameters.extend(["--norm-method", str(norm_method)])

    command_parameters.extend(other_parameter)

    command = ["mageck mle"]
    command.extend(command_parameters)
    cmd = " ".join(command)
    print(cmd)

    mageck_count = subprocess.run(
        cmd, capture_output=True, text=True, shell=True
    )
    print(mageck_count.stdout)
    print(mageck_count.stderr)

    gene_summary_txt = Path(f"{out_dir}/{prefix}.gene_summary.txt")
    gene_summary_tsv = Path(f"{out_dir}/{prefix}.gene_summary.tsv")

    sgrna_summary_txt = Path(f"{out_dir}/{prefix}.sgrna_summary.txt")
    sgrna_summary_tsv = Path(f"{out_dir}/{prefix}.sgrna_summary.tsv")

    if gene_summary_txt.is_file():
        shutil.copy(gene_summary_txt, gene_summary_tsv)
    if sgrna_summary_txt.is_file():
        shutil.copy(sgrna_summary_txt, sgrna_summary_tsv)


def split_frame_to_control_and_query(
    mageck_frame: DataFrame,
    control_prefix: str,
    id_col: Optional[str] = None,
    name_column: str = "Name",
    sgRNA_column: str = "sgRNA",
    infer_genes: Optional[Callable] = None,
) -> Dict[str, DataFrame]:

    if id_col is None:
        mageck_frame["id"] = "s_" + mageck_frame.index.astype(str)
    else:
        mageck_frame["id"] = mageck_frame[id_col].str.replace(
            r"\s+", "_", regex=True
        )

    if infer_genes is not None:
        mageck_frame[name_column] = infer_genes(mageck_frame)
    else:
        mageck_frame[name_column] = mageck_frame[name_column].str.replace(
            r"\s+", "_", regex=True
        )
    control_rows = mageck_frame[
        mageck_frame[name_column].str.startswith(control_prefix)
    ].index
    df_control = mageck_frame.loc[control_rows][["id"]].copy()
    df_query = mageck_frame[["id", sgRNA_column, name_column]].copy()
    return {"control": df_control, "query": df_query}


def mageck_pathway(
    gene_ranking: Union[Path, str],
    gmt_file: Union[Path, str],
    out_dir: Union[Path, str],
    prefix: str = "pathway",
    method: str = "gsea",
    single_ranking: bool = False,
    output_prefix: Optional[str] = None,
    sort_criteria: str = "neg",
    keep_tmp: bool = False,
    ranking_column: Optional[Union[str, int]] = None,
    ranking_column_2: Optional[Union[str, int]] = None,
    pathway_alpha: Optional[float] = None,
    permutation: Optional[int] = None,
    other_parameter: Dict[str, str] = [],
):
    """
    Wrapper for `mageck pathway` subcommand.

    Required:
      - gene_ranking: gene ranking file
      - gmt_file: GMT pathways file

    Optional parameters mirror MAGeCK CLI options.

    Returns a dict with stdout/stderr and list of output files in out_dir
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    command_parameters = [
        "--gene-ranking",
        str(gene_ranking),
        "--gmt-file",
        str(gmt_file),
    ]

    # output prefix
    out_pref = output_prefix if output_prefix is not None else prefix
    command_parameters.extend(["-n", f"{out_dir}/{out_pref}"])

    # method
    if method is not None:
        command_parameters.extend(["--method", str(method)])

    if single_ranking:
        command_parameters.append("--single-ranking")

    if sort_criteria is not None:
        command_parameters.extend(["--sort-criteria", str(sort_criteria)])

    if keep_tmp:
        command_parameters.append("--keep-tmp")

    if ranking_column is not None:
        command_parameters.extend(["--ranking-column", str(ranking_column)])
    if ranking_column_2 is not None:
        command_parameters.extend(["--ranking-column-2", str(ranking_column_2)])
    if pathway_alpha is not None:
        command_parameters.extend(["--pathway-alpha", str(pathway_alpha)])
    if permutation is not None:
        command_parameters.extend(["--permutation", str(permutation)])

    for k in other_parameter:
        command_parameters.extend([k, other_parameter[k]])

    command = ["mageck pathway"]
    command.extend(command_parameters)
    cmd = " ".join(command)
    print(cmd)

    proc = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    print(proc.stdout)
    print(proc.stderr)

    # Collect output files with provided prefix
    outputs = list(Path(out_dir).glob(f"{out_pref}*"))

    return {
        "cmd": cmd,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "outputs": [str(p) for p in outputs],
    }


def mageck_plot(
    gene_summary: Optional[Union[Path, str]] = None,
    sgrna_summary: Optional[Union[Path, str]] = None,
    out_dir: Union[Path, str] = ".",
    prefix: str = "plot",
    other_parameter: Dict[str, str] = [],
):
    """
    Generic wrapper for `mageck plot` subcommand.

    It will pass any provided files to the CLI and collect generated plots.
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    command_parameters = []
    if gene_summary is not None:
        command_parameters.extend(["-k", str(gene_summary)])
    if sgrna_summary is not None:
        command_parameters.extend(["-s", str(sgrna_summary)])

    command_parameters.extend(["-n", f"{out_dir}/{prefix}"])

    for k in other_parameter:
        command_parameters.extend([k, other_parameter[k]])

    command = ["mageck plot"]
    command.extend(command_parameters)
    cmd = " ".join(command)
    print(cmd)

    proc = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    print(proc.stdout)
    print(proc.stderr)

    outputs = list(Path(out_dir).glob(f"{prefix}*"))
    return {
        "cmd": cmd,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "outputs": [str(p) for p in outputs],
    }


def combine_gene_info_with_mageck_output(
    df_mageck: DataFrame,
    df_genes: DataFrame,
    name_column_mageck: str = "id",
    name_column_genes: str = "Gene",
    how: str = "left",
    columns_to_add: List[str] = [
        "gene_stable_id",
        "name",
        "chr",
        "start",
        "stop",
        "strand",
        "biotype",
    ],
) -> DataFrame:
    """
    Combine gene information dataframe with MAGeCK output dataframe.

    Parameters:
    - df_mageck: DataFrame containing MAGeCK results.
    - df_genes: DataFrame containing gene information.
    - name_column_mageck: Column name in df_mageck to match with df_genes.
    - name_column_genes: Column name in df_genes to match with df_mageck.

    Returns:
    - Merged DataFrame with gene information added to MAGeCK results.
    """
    df_genes = df_genes.drop_duplicates(keep="first", subset=name_column_genes)
    merged_df = df_mageck.merge(
        df_genes[columns_to_add + [name_column_genes]],
        left_on=name_column_mageck,
        right_on=name_column_genes,
        how=how,
    )
    merged_df.drop(columns=[name_column_genes], inplace=True)
    return merged_df


def get_significant_genes(
    df_mageck: DataFrame,
    fdr_column: str = "pos|fdr",
    fdr_threshold: float = 0.05,
    logfc_column: str = "pos|lfc",
    logfc_threshold: float = 1.0,
    direction: str = "both",  # "both", "pos", "neg"
) -> DataFrame:
    """
    get_significant_genes filters the gene summary table for significant hits.

        Parameters
    ----------
    df_mageck : DataFrame
        Mageck output dataframe, gene summary table.
    fdr_column : str, optional
        False Discovery Rate column, by default "pos|fdr"
    fdr_threshold : float, optional
        maximum FDR, by default 0.05
    logfc_column : str, optional
        logFC column, by default "pos|lfc"
    logfc_threshold : float, optional
        minimum fold change, by default 1.0
    direction : str, optional
        positive or nagative selection, by default "both"

    Returns
    -------
    DataFrame
        Filtered dataframe with significant genes only.
    """
    if direction == "both":
        sig_genes = df_mageck[
            (df_mageck[fdr_column] <= fdr_threshold)
            & (df_mageck[logfc_column].abs() >= logfc_threshold)
        ]
    elif direction == "pos":
        sig_genes = df_mageck[
            (df_mageck[fdr_column] <= fdr_threshold)
            & (df_mageck[logfc_column] >= logfc_threshold)
        ]
    elif direction == "neg":
        sig_genes = df_mageck[
            (df_mageck[fdr_column] <= fdr_threshold)
            & (df_mageck[logfc_column] <= -logfc_threshold)
        ]
    else:
        raise ValueError("direction must be one of: 'both', 'pos', or 'neg'")
    return sig_genes
