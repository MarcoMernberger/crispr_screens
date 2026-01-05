import itertools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from typing import Dict, Iterable, Tuple, Union, List, Optional
from pandas import DataFrame

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


def volcano_plot(
    df: pd.DataFrame,
    log_fc_column: str,
    fdr_column: Union[str, Tuple[str, str]],
    *,
    name_column: Optional[str] = None,
    top_n_labels: int = 0,
    transform_y: bool = True,  # True -> plot -log10(FDR)
    log_threshold: float = 1.0,
    fdr_threshold: float = 0.05,
    point_size: float = 12,
    alpha: float = 0.75,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    figsize: Tuple[float, float] = (8, 6),
    y_clip_min: float = 1e-300,
    label_fontsize: int = 9,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Volcano plot supporting either:
      - a single FDR column (MLE-style), or
      - separate FDR columns for positive/negative effects (RRA-style).

    If fdr_column is a tuple:
        (fdr_pos, fdr_neg)
    then:
        logFC >= 0 -> fdr_pos
        logFC <  0 -> fdr_neg
    """

    if log_fc_column not in df.columns:
        raise KeyError(f"Missing column: {log_fc_column}")

    if isinstance(fdr_column, tuple):
        if len(fdr_column) != 2:
            raise ValueError("fdr_column tuple must be (pos_fdr, neg_fdr)")
        fdr_pos_col, fdr_neg_col = fdr_column
        for col in (fdr_pos_col, fdr_neg_col):
            if col not in df.columns:
                raise KeyError(f"Missing column: {col}")
    else:
        if fdr_column not in df.columns:
            raise KeyError(f"Missing column: {fdr_column}")
        fdr_pos_col = fdr_column
        fdr_neg_col = None

    cols = [log_fc_column, fdr_pos_col]
    if fdr_neg_col is not None:
        cols.append(fdr_neg_col)
    if name_column:
        cols.append(name_column)

    data = df[cols].copy()

    # Ensure numeric
    data[log_fc_column] = pd.to_numeric(data[log_fc_column], errors="coerce")
    for col in [fdr_pos_col, fdr_neg_col]:
        if col is not None:
            data[col] = pd.to_numeric(data[col], errors="coerce")

    data = data.dropna(subset=[log_fc_column])
    x = data[log_fc_column].to_numpy()

    # Select correct FDR depending on direction
    if fdr_neg_col is not None:
        fdr_raw = np.where(
            x >= 0,
            data[fdr_pos_col].to_numpy(),
            data[fdr_neg_col].to_numpy(),
        )
    else:
        fdr_raw = data[fdr_pos_col].to_numpy()

    fdr_raw = np.clip(fdr_raw, y_clip_min, 1.0)

    if transform_y:
        y = -np.log10(fdr_raw)
        y_thr_plot = -np.log10(max(fdr_threshold, y_clip_min))
        y_label = ylabel or r"$-\log_{10}(\mathrm{FDR})$"
    else:
        y = fdr_raw
        y_thr_plot = fdr_threshold
        y_label = ylabel or "FDR"

    # Significance logic (RAW FDR scale)
    sig = (fdr_raw <= fdr_threshold) & (np.abs(x) >= log_threshold)
    pos = sig & (x >= log_threshold)
    neg = sig & (x <= -log_threshold)
    nonsig = ~sig

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    # Non-significant first
    ax.scatter(
        x[nonsig],
        y[nonsig],
        s=point_size,
        alpha=alpha,
        c="grey",
        edgecolors="none",
    )

    # Significant
    ax.scatter(
        x[neg], y[neg], s=point_size, alpha=alpha, c="blue", edgecolors="none"
    )
    ax.scatter(
        x[pos], y[pos], s=point_size, alpha=alpha, c="red", edgecolors="none"
    )

    # Threshold lines
    ax.axvline(+log_threshold, linestyle="--", linewidth=1, color="grey")
    ax.axvline(-log_threshold, linestyle="--", linewidth=1, color="grey")
    ax.axhline(y_thr_plot, linestyle="--", linewidth=1, color="grey")

    ax.set_xlabel(xlabel or log_fc_column)
    ax.set_ylabel(y_label)
    if title:
        ax.set_title(title)

    # Optional labeling (top by smallest relevant FDR)
    if top_n_labels > 0:
        if not name_column:
            raise ValueError("top_n_labels requires name_column")

        label_df = data.loc[sig].copy()
        label_df["_fdr_used"] = fdr_raw[sig]
        label_df = pd.concat(
            [
                label_df.sort_values("_fdr_used").head(top_n_labels),
                label_df.sort_values("_fdr_used").tail(top_n_labels),
            ]
        )

        for _, row in label_df.iterrows():
            xi = row[log_fc_column]
            fdr_i = max(min(row["_fdr_used"], 1.0), y_clip_min)
            yi = -np.log10(fdr_i) if transform_y else fdr_i
            ax.text(xi, yi, str(row[name_column]), fontsize=label_fontsize)

    ax.margins(x=0.05, y=0.05)
    return fig, ax
