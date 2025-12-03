import pandas as pd
import re
import shutil
import subprocess
import rpy2.robjects as ro
from typing import Dict, Optional, Union, List
from pandas import DataFrame
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

        fdr_col = [x for x in fdr_cols if "wald-fdr" in x][0]
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
    norm_method: str = None,
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
        count_table,
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

    command_parameters.extend(other_parameter)

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
