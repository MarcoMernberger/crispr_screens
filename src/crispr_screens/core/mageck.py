import pandas as pd
import re
from typing import Dict, Optional, Union, List
from pandas import DataFrame


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


# def get_top_genes_from_selected_groups(
#     scatter_frame: DataFrame, group: str = "top", gene_id_column: str = "id"
# ) -> DataFrame:
#     scatter_frame_selected = scatter_frame[scatter_frame["group"] == "top"]
