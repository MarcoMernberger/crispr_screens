import pandas as pd
from typing import Dict, Optional, Union, List, Tuple
from pathlib import Path
from crispr_screens.core.mageck import (
    filter_multiple_mageck_comparisons,
    combine_comparisons,
)
from crispr_screens.core.plots import plot_selected_venn


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


def write_venn(
    outdir: Union[Path, str],
    filebasename: str,
    label_to_file: Dict[str, str],
    id_cols: Union[List[str], str] = "id",
    sep: str = "\t",
    figsize: Tuple[float, float] = (6, 6),
    title: str | None = None,
):
    outdir = Path(outdir)
    outdir.parent.mkdir(exist_ok=True, parents=True)
    out_venn = f"{filebasename}_venn"
    out_dataframe = outdir / f"{filebasename}.tsv"
    results = plot_selected_venn(
        label_to_file=label_to_file,
        id_cols=id_cols,
        sep=sep,
        figsize=figsize,
        title=title,
    )
    save_figure(results["figure"], outdir, out_venn)
    results["memberships"].to_csv(out_dataframe, sep="t", index=False)


def save_figure(f, folder, name, bbox_inches="tight"):
    folder.mkdir(exist_ok=True, parents=True)
    for suffix in [".png", ".svg", ".pdf"]:
        f.savefig(folder / (name + suffix), bbox_inches=bbox_inches)
