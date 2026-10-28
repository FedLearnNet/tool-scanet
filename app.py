import pathlib
from pathlib import Path
from typing import Optional

import numpy as np
import scanet as sn
import scanpy as sc
from anndata import AnnData
from pyfedappwrap.learning.app_types import BaseSelfLearnedApp
from scipy import sparse
from typing_extensions import override

from config import MyAppConfig, MyAppInputConfig, MyAppOutputConfig


class MyAPP(BaseSelfLearnedApp[MyAppConfig, MyAppInputConfig, MyAppOutputConfig]):

    def __init__(self):
        super().__init__()

    def _rename_figures(self, out_dir: Path, renames: dict[str, str]) -> None:
        for old_name, new_name in renames.items():
            old_path = out_dir / old_name
            new_path = out_dir / new_name
            if old_path.exists():
                self.logger.info(f"Renaming figure: {old_path.name} -> {new_path.name}")
                old_path.rename(new_path)
            else:
                self.logger.info(f"Figure not found, skipping rename: {old_path.name}")

    def _ensure_expected_figure(
            self,
            out_dir: Path,
            expected_name: str,
            candidates: list[str],
    ) -> Path:
        expected_path = out_dir / expected_name
        if expected_path.exists():
            self.logger.info(f"Expected figure already exists: {expected_path.name}")
            return expected_path

        for candidate in candidates:
            candidate_path = out_dir / candidate
            if candidate_path.exists():
                self.logger.info(
                    f"Normalizing figure filename: {candidate_path.name} -> {expected_path.name}"
                )
                candidate_path.rename(expected_path)
                return expected_path

        matching_files = sorted(out_dir.glob("*.png"))
        for file_path in matching_files:
            lower_name = file_path.name.lower()
            if "highly" in lower_name and "variable" in lower_name:
                self.logger.info(
                    f"Using discovered highly-variable-genes figure: {file_path.name} -> {expected_path.name}"
                )
                file_path.rename(expected_path)
                return expected_path

        self.logger.warning(
            f"Could not find figure for expected output '{expected_name}'. Checked candidates: {candidates}"
        )
        return expected_path

    def _preprocess_and_cluster(
            self,
            adata: AnnData,
            n_pcs: int,
            n_neighbors: int,
            method: str,
            cell_annotation: str,
    ) -> str:
        self.logger.info(
            f"Starting dimensionality reduction and clustering "
            f"(n_pcs={n_pcs}, n_neighbors={n_neighbors}, method={method}, annotation={cell_annotation})"
        )

        n_comp = min(n_pcs, min(adata.n_obs, adata.n_vars) - 1)
        if n_comp < 1:
            n_comp = 1

        self.logger.info(f"Using {n_comp} principal components")
        sc.pp.pca(adata, svd_solver="arpack", n_comps=n_comp)

        self.logger.info("Computing neighborhood graph")
        sc.pp.neighbors(adata, n_pcs=n_comp, n_neighbors=n_neighbors)

        if method == "UMAP":
            self.logger.info("Running UMAP embedding")
            sc.tl.umap(adata)
        else:
            self.logger.info("Using the PCA embedding computed above")

        self.logger.info("Running Leiden clustering")
        sc.tl.leiden(adata, resolution=0.5)

        return "leiden"

    def _prepare_for_representative_cells(self, adata: AnnData, cluster_col: str) -> AnnData:
        self.logger.info("Preparing AnnData for SCANet representative_cells")

        if cluster_col not in adata.obs.columns:
            raise ValueError(f"Cluster column '{cluster_col}' not found in adata.obs")

        pb_adata = adata.copy()
        self.logger.info(f"Copied AnnData for pseudobulk step. Matrix type: {type(pb_adata.X).__name__}")

        if sparse.issparse(pb_adata.X):
            self.logger.info(
                "Converting sparse expression matrix to dense float64 array because SCANet pseudobulk "
                "internally slices adata.X and otherwise produces SparseCSRMatrixView objects"
            )
            pb_adata.X = pb_adata.X.toarray().astype(np.float64, copy=False)
        else:
            self.logger.info("Converting dense expression matrix to contiguous float64 numpy array for compatibility")
            pb_adata.X = np.asarray(pb_adata.X, dtype=np.float64)

        pb_adata.obs = pb_adata.obs.copy()
        pb_adata.var = pb_adata.var.copy()

        self.logger.info(
            f"Prepared AnnData for representative_cells. Matrix type is now: {type(pb_adata.X).__name__}, "
            f"dtype: {pb_adata.X.dtype}"
        )
        return pb_adata

    def run_algorithm(self, data: MyAppInputConfig) -> MyAppOutputConfig:

        self.logger.info("Starting SCANet app run")

        out_dir = Path("./output")
        out_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Using output directory: {out_dir.resolve()}")

        sc.settings.figdir = str(out_dir) + "/"
        sc.settings.verbosity = 1

        self.logger.info(
            "Loaded hyperparameters: "
            f"lowerbound={self.config.lowerbound}, "
            f"upperbound={self.config.upperbound}, "
            f"min_cell={self.config.min_cell}, "
            f"threshold={self.config.threshold}, "
            f"top_n_genes={self.config.top_n_genes}, "
            f"num_rep_cells={self.config.num_rep_cells}, "
            f"n_pcs={self.config.n_pcs}, "
            f"n_neighbors={self.config.n_neighbors}, "
            f"method={self.config.method}, "
            f"cell_annotation={self.config.cell_annotation}"
        )

        self.logger.info("Reading input h5ad file with SCANet preprocessing")
        input_path = Path(data.data)
        if input_path.suffix != ".h5ad":
            renamed_input_path = input_path.parent / f"{input_path.name}.h5ad"
            self.logger.info(
                f"Input path has suffix '{input_path.suffix}', renaming file to '{renamed_input_path.name}'"
            )
            input_path.rename(renamed_input_path)
            input_path = renamed_input_path
            data.data = input_path
            self.logger.info(f"Updated input path to: {str(data.data)}")

        adata = sn.pp.read_h5ad(data.data, pr_process="Yes")
        self.logger.info(f"Loaded AnnData with {adata.n_obs} cells and {adata.n_vars} genes")

        self.logger.info(
            f"Filtering cells by number of genes: lowerbound={self.config.lowerbound}, upperbound={self.config.upperbound}"
        )
        adata = sn.pp.filter_cells_by_n_genes(adata, self.config.lowerbound, self.config.upperbound)
        self.logger.info(f"After cell filtering: {adata.n_obs} cells, {adata.n_vars} genes")

        mt_pct = int(self.config.threshold * 100)
        self.logger.info(f"Filtering mitochondrial genes with threshold={mt_pct}%")
        adata = sn.pp.filter_mitochondrial(adata, mt_pct)
        self.logger.info(
            f"After mitochondrial filtering: {adata.n_obs} cells, {adata.n_vars} genes")

        self.logger.info(f"Filtering genes with min_cell={self.config.min_cell}")
        adata = sn.pp.filter_genes(adata, min_cells=self.config.min_cell)
        self.logger.info(f"After gene filtering: {adata.n_obs} cells, {adata.n_vars} genes")

        if adata.n_obs == 0:
            raise ValueError("No cells remaining after gene filtering")
        if adata.n_vars == 0:
            raise ValueError("No genes remaining after filtering")

        self.logger.info("Applying log normalization")
        adata = sn.pp.log_normalize(adata, target_sum=1e6, use_log=True)

        self.logger.info(f"Selecting top {self.config.top_n_genes} highly variable genes")
        adata = sn.pp.extract_highly_variable_genes(
            adata,
            n_top_genes=self.config.top_n_genes,
            plot=False,
        )

        self.logger.info("Generating highly variable genes plot")
        sc.pl.highly_variable_genes(adata, show=False, save="_highly_variable_genes.png")

        self.logger.info("Generating QC violin plot")
        for candidate_set in (
                ["n_genes", "n_counts", "pct_counts_mt"],
                ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                ["n_genes_by_counts", "total_counts", "percent_mito"],
        ):
            qc_keys = [k for k in candidate_set if k in adata.obs.columns]
            if qc_keys:
                self.logger.info(f"Using QC keys for violin plot: {qc_keys}")
                sc.pl.violin(adata, keys=qc_keys, jitter=0.4, save="_qc_violin.png")
                break

        cluster_col = self._preprocess_and_cluster(
            adata=adata,
            n_pcs=self.config.n_pcs,
            n_neighbors=self.config.n_neighbors,
            method=self.config.method,
            cell_annotation=self.config.cell_annotation,
        )
        self.logger.info(f"Clustering finished. Using cluster column: {cluster_col}")

        cluster_sizes = adata.obs[cluster_col].value_counts()
        self.logger.info(f"Cluster sizes before representative cell selection: {cluster_sizes.to_dict()}")

        processed_path = out_dir / "processed.h5ad"
        self.logger.info(f"Writing processed AnnData to {processed_path}")
        adata.write_h5ad(processed_path)

        if self.config.method == "UMAP" and "X_umap" in adata.obsm:
            self.logger.info("Generating UMAP cell annotation plot")
            sc.pl.umap(adata, color=[cluster_col], save="_cell_annotation.png")

            self.logger.info("Generating UMAP cluster plot")
            sc.pl.umap(adata, color=[cluster_col], save="_SCANclusters.png")
        elif "X_pca" in adata.obsm:
            self.logger.info("Generating PCA plots")
            sc.pl.pca(adata, color=[cluster_col], save="_cell_annotation.png")
            sc.pl.pca(adata, color=[cluster_col], save="_SCANclusters.png")

        cluster_sizes = adata.obs[cluster_col].value_counts()
        too_small_clusters = cluster_sizes[cluster_sizes < self.config.num_rep_cells]
        if not too_small_clusters.empty:
            self.logger.info(
                "The following clusters are smaller than num_rep_cells and will be removed by SCANet: "
                f"{too_small_clusters.to_dict()}"
            )

        self.logger.info(
            f"Computing representative cells with num_rep_cells={self.config.num_rep_cells}"
        )
        pb_adata = self._prepare_for_representative_cells(adata, cluster_col)
        self.logger.info(
            f"Calling SCANet representative_cells with matrix type={type(pb_adata.X).__name__}, "
            f"dtype={pb_adata.X.dtype}, shape={pb_adata.X.shape}"
        )
        reduced_data = sn.pb.representative_cells(pb_adata, self.config.num_rep_cells, cluster_col)

        reduced_path = out_dir / "reduced.h5ad"
        self.logger.info(f"Writing reduced data to {reduced_path}")
        reduced_data.write_h5ad(reduced_path)

        self._rename_figures(
            out_dir,
            {
                "highly_variable_genes_highly_variable_genes.png": "highly_variable_genes.png",
                "violin_qc_violin.png": "violin_plot.png",
                "umap_cell_annotation.png": "cell_annotation.png",
                "umap_SCANclusters.png": "cluster.png",
                "pca_cell_annotation.png": "cell_annotation.png",
                "pca_SCANclusters.png": "cluster.png",
            },
        )
        self._ensure_expected_figure(
            out_dir,
            "highly_variable_genes.png",
            [
                "highly_variable_genes_highly_variable_genes.png",
                "filter_genes_dispersion_highly_variable_genes.png",
                "highly_variable_genes.png",
            ],
        )

        self.logger.info("Preparing app output payload")
        result = MyAppOutputConfig(
            reduced_data=reduced_path,
            cell_annotation=out_dir / "cell_annotation.png",
            cluster=out_dir / "cluster.png",
            highly_variable_genes=out_dir / "highly_variable_genes.png",
            violin_plot=out_dir / "violin_plot.png",
            processed_data=processed_path,
        )

        self.logger.info("SCANet app run completed successfully")
        return result

    @override
    def get_test_data(self) -> Optional[dict[str, Path]]:
        return {"data": "data/pbmc3k.h5ad"}
