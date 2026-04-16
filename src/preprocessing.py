# src/preprocessing.py

import scanpy as sc
import numpy as np


def calculate_qc_metrics(adata, mito_prefix: str = "MT-"):
    """
    Annotate mitochondrial genes and compute QC metrics.

    Parameters:
        adata (AnnData): Annotated data matrix
        mito_prefix (str): Prefix for mitochondrial genes

    Returns:
        AnnData: Updated object with QC metrics
    """

    # Identify mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)

    # Compute QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    return adata


def filter_cells_genes(
    adata,
    min_genes: int = 200,
    min_cells: int = 3,
):
    """
    Filter low-quality cells and genes.

    Parameters:
        adata (AnnData)
        min_genes (int): Minimum genes per cell
        min_cells (int): Minimum cells per gene

    Returns:
        AnnData
    """

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    return adata


def filter_by_qc(
    adata,
    max_mt_pct: float = 20.0,
    min_counts: int = None,
    max_counts: int = None,
):
    """
    Filter cells based on QC thresholds.

    Parameters:
        adata (AnnData)
        max_mt_pct (float): Max mitochondrial percentage
        min_counts (int): Minimum counts per cell
        max_counts (int): Maximum counts per cell

    Returns:
        AnnData
    """

    # Mitochondrial filtering
    adata = adata[adata.obs.pct_counts_mt < max_mt_pct].copy()

    # Optional count filtering
    if min_counts is not None:
        adata = adata[adata.obs.total_counts > min_counts].copy()

    if max_counts is not None:
        adata = adata[adata.obs.total_counts < max_counts].copy()

    return adata


def normalize_log1p(adata, target_sum: float = 1e4):
    """
    Normalize counts per cell and log-transform.

    Parameters:
        adata (AnnData)
        target_sum (float): Scaling factor

    Returns:
        AnnData
    """

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    return adata


def highly_variable_genes(
    adata,
    n_top_genes: int = 2000,
    flavor: str = "seurat",
):
    """
    Identify highly variable genes (HVGs).

    Parameters:
        adata (AnnData)
        n_top_genes (int)
        flavor (str): 'seurat', 'cell_ranger', etc.

    Returns:
        AnnData
    """

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
    )

    # Subset to HVGs
    adata = adata[:, adata.var.highly_variable].copy()

    return adata


def scale_data(adata, max_value: float = 10):
    """
    Scale data to unit variance.

    Parameters:
        adata (AnnData)
        max_value (float): Clip values

    Returns:
        AnnData
    """

    sc.pp.scale(adata, max_value=max_value)
    return adata


def run_qc_pipeline(
    adata,
    mito_prefix: str = "MT-",
    min_genes: int = 200,
    min_cells: int = 3,
    max_mt_pct: float = 20.0,
    n_top_genes: int = 2000,
):
    """
    Full preprocessing pipeline:
    QC → filtering → normalization → HVG → scaling

    Parameters:
        adata (AnnData)

    Returns:
        AnnData
    """

    # Step 1: QC metrics
    adata = calculate_qc_metrics(adata, mito_prefix=mito_prefix)

    # Step 2: Basic filtering
    adata = filter_cells_genes(adata, min_genes=min_genes, min_cells=min_cells)

    # Step 3: QC-based filtering
    adata = filter_by_qc(adata, max_mt_pct=max_mt_pct)

    # Step 4: Normalize + log
    adata = normalize_log1p(adata)
    adata.raw = adata.copy()

    # Step 5: Highly variable genes
    adata = highly_variable_genes(adata, n_top_genes=n_top_genes)

    # Step 6: Scaling
    adata = scale_data(adata)

    return adata


def plot_qc(adata, save: bool = False):
    """
    Plot QC metrics.

    Parameters:
        adata (AnnData)
        save (bool)

    Returns:
        None
    """

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save="_qc_violin.png" if save else None,
        show=not save,
    )

    sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_mt",
        save="_qc_scatter.png" if save else None,
        show=not save,
    )