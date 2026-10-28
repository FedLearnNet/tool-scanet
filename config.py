from dataclasses import dataclass
from pathlib import Path
from typing import Any

from pyfedappwrap.learning.run_runfig import AppConfig, AppInputConfig, AppOutputConfig


@dataclass
class MyAppConfig(AppConfig):
    cell_annotation: str = "leiden"
    num_rep_cells: int = 50
    n_pcs: int = 10
    n_neighbors: int = 5
    method: str = "UMAP"
    lowerbound: int = 500
    upperbound: int = 5000
    nbins: int = 10
    min_cell: int = 3
    threshold: float = 0.1
    top_n_genes: int = 500

    def __post_init__(self):
        if self.method not in ["PCA", "UMAP"]:
            raise ValueError("method must be either 'PCA' or 'UMAP'")
        if self.num_rep_cells <= 0:
            raise ValueError("num_rep_cells must be positive")
        if self.n_pcs <= 0:
            raise ValueError("n_pcs must be positive")
        if self.n_neighbors <= 0:
            raise ValueError("n_neighbors must be positive")
        if not (0 < self.threshold <= 1):
            raise ValueError("threshold must be between 0 and 1")
        if self.min_cell <= 0:
            raise ValueError("min_cell must be positive")
        if self.top_n_genes <= 0:
            raise ValueError("top_n_genes must be positive")


@dataclass
class MyAppInputConfig(AppInputConfig):
    data: Any = None


@dataclass
class MyAppOutputConfig(AppOutputConfig):
    cell_annotation: Path | None = None
    cluster: Path | None = None
    highly_variable_genes: Path | None = None
    violin_plot: Path | None = None
    processed_data: Path | None = None
    reduced_data: Path | None = None
