import itertools
from pathlib import Path
from typing import Dict, Iterable, Tuple, Union, List

import matplotlib.pyplot as plt
import pandas as pd

try:
    from matplotlib_venn import venn2, venn3

    _HAS_MPL_VENN = True
except ImportError:
    _HAS_MPL_VENN = False


def plot_selected_venn(
    label_to_file: Dict[str, str],
    id_cols: Union[List[str], str] = "id",
    sep: str = "\t",
    figsize: Tuple[float, float] = (6, 6),
    title: str | None = None,
):
    """
    Plot Venn/UpSet-style diagram from 2–4 selected.tsv files.

    Parameters
    ----------
    label_to_file : dict
        Mapping label -> path to selected.tsv
        (values are read with pandas.read_csv(..., sep=sep)).
    id_col : str
        Column name that contains gene IDs.
    sep : str
        Separator for the TSV files.
    figsize : (float, float)
        Figure size for matplotlib.
    title : str or None
        Optional plot title.

    Returns
    -------
    result : dict
        {
          "sets": dict[label -> set(ids)],
          "memberships": pd.DataFrame with boolean membership per label,
          "figure": matplotlib.figure.Figure
        }
    """
    if not (2 <= len(label_to_file) <= 4):
        raise ValueError("Please provide between 2 and 4 files (labels).")

    sets: Dict[str, set] = {}
    ii = 0
    for label, path in label_to_file.items():
        if isinstance(id_cols, str):
            id_col = id_cols
        else:
            id_col = id_cols[ii]
        df = pd.read_csv(path, sep=sep)
        if id_col not in df.columns:
            raise ValueError(
                f"Column '{id_col}' not found in file '{path}'. "
                f"Available columns: {list(df.columns)}"
            )
        ids = set(df[id_col].dropna().astype(str))
        sets[label] = ids
        ii += 1

    labels = list(sets.keys())
    n = len(labels)

    all_ids = sorted(set().union(*sets.values()))
    membership_data = {"id": all_ids}
    for label in labels:
        membership_data[label] = [gene in sets[label] for gene in all_ids]

    memberships = pd.DataFrame(membership_data).set_index("id")

    fig = plt.figure(figsize=figsize)

    if n in (2, 3) and _HAS_MPL_VENN:
        ax = fig.add_subplot(111)
        set_list = [sets[l] for l in labels]

        if n == 2:
            venn2(set_list, set_labels=labels, ax=ax)
        else:
            venn3(set_list, set_labels=labels, ax=ax)

        if title is not None:
            ax.set_title(title)

    else:
        ax = fig.add_subplot(111)

        combos: Iterable[Tuple[str, ...]] = []
        for r in range(1, n + 1):
            combos += itertools.combinations(labels, r)

        combo_counts = []
        combo_labels = []
        for combo in combos:
            mask = memberships[list(combo)].all(axis=1)
            count = mask.sum()
            if count > 0:
                combo_counts.append(count)
                combo_labels.append("&".join(combo))

        order = sorted(
            range(len(combo_counts)),
            key=lambda i: combo_counts[i],
            reverse=True,
        )

        combo_counts = [combo_counts[i] for i in order]
        combo_labels = [combo_labels[i] for i in order]

        ax.bar(range(len(combo_counts)), combo_counts)
        ax.set_xticks(range(len(combo_labels)))
        ax.set_xticklabels(combo_labels, rotation=45, ha="right")
        ax.set_ylabel("Number of genes")
        if title is None:
            ax.set_title("Set intersections (UpSet-style)")
        else:
            ax.set_title(title)

        fig.tight_layout()

    result = {
        "sets": sets,
        "memberships": memberships,
        "figure": fig,
    }
    return result
