# tool-scanet

Preprocesses single-cell expression data with SCANet and exports cell clusters, embeddings, and QC plots.

## Runtime inputs

Provide an AnnData `.h5ad` file through the `data` input. SCANet preprocessing
expects a sparse expression matrix in `adata.X` (CSR format is supported).

## Inputs and outputs

### Inputs

| Name | Type | Description |
| --- | --- | --- |
| `data` | PATH | Path to h5ad file (e.g. ./data/test_data.h5ad) |

### Outputs

| Name | Type | Description |
| --- | --- | --- |
| `reduced_data` | PATH | Reduced dimensionality data file |
| `cell_annotation` | IMAGE | Cell annotation visualization |
| `cluster` | IMAGE | Cluster visualization |
| `highly_variable_genes` | IMAGE | Highly variable genes plot |
| `violin_plot` | IMAGE | Violin plot of QC metrics (n_genes, total_counts, pct_mito) |
| `processed_data` | PATH | Processed data file |

## Configuration

| Name | Default | Description |
| --- | --- | --- |
| `min_cell` | `3` | Minimum number of cells for gene filtering |
| `threshold` | `0.1` | Threshold for mitochondrial gene filtering |
| `top_n_genes` | `500` | Number of top genes to select |
| `cell_annotation` | `leiden` | Cell annotation method for clustering |
| `num_rep_cells` | `50` | Number of representative cells per cluster |
| `n_pcs` | `10` | Number of principal components |
| `n_neighbors` | `5` | Number of neighbors for UMAP |
| `method` | `UMAP` | Dimensionality reduction method |
| `lowerbound` | `500` | Lower bound for cell filtering |
| `upperbound` | `5000` | Upper bound for cell filtering (max genes per cell) |
| `nbins` | `10` | Number of bins for histogram |

## Development

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-dev.txt
cp .env.example .env
# Set APP_ID and platform URLs before connecting.
python main.py
```

The Docker base image provides the FL-Net API. `requirements.txt` lists additional
runtime dependencies; `requirements-dev.txt` includes the API for local development.
Runtime datasets and trained models must be supplied separately; see `data/README.md`.

Build with `docker build -t tool-scanet .`.
