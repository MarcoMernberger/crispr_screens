import pandas as pd
from pathlib import Path
from typing import Union
from pandas import DataFrame


def save_figure(f, folder, name, bbox_inches="tight"):
    folder.mkdir(exist_ok=True, parents=True)
    for suffix in [".png", ".svg", ".pdf"]:
        f.savefig(folder / (name + suffix), bbox_inches=bbox_inches)


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
