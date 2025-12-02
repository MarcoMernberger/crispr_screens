"""
CRISPR screens – demultiplexing and read processing tools.

This package provides:
- core
- models
- services
- jobs
- api
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("crispr-screens")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.1.0"


from .jobs.mageck_jobs import (
    combine_comparison_output_job,
    write_filtered_mageck_comparison_job,
    run_mageck_scatterview_job,
    write_venn_job,
)

__all__ = [
    "combine_comparison_output_job",
    "write_filtered_mageck_comparison_job",
    "run_mageck_scatterview_job",
    "write_venn_job",
]
