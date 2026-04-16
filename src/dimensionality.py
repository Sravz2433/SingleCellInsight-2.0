# src/dimensionality.py

import scanpy as sc


def run_pca(
    adata,
    n_comps: int = 50,
    svd_solver: str = "arpack",
):
    """
    Perform Principal Component Analysis (PCA).

    Parameters:
        adata (AnnData)
        n_comps (int): Number of principal components
        svd_solver (str): Solver type

    Returns:
        AnnData
    """

    sc.tl.pca(adata, n_comps=n_comps, svd_solver=svd_solver)
    return adata


def compute_neighbors(
    adata,
    n_neighbors: int = 15,
    n_pcs: int = 50,
    metric: str = "euclidean",
):
    """
    Compute neighborhood graph (required for UMAP & clustering).

    Parameters:
        adata (AnnData)
        n_neighbors (int)
        n_pcs (int)
        metric (str)

    Returns:
        AnnData
    """

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, metric=metric)
    return adata


def run_umap(
    adata,
    min_dist: float = 0.5,
    spread: float = 1.0,
):
    """
    Compute UMAP embedding.

    Parameters:
        adata (AnnData)
        min_dist (float)
        spread (float)

    Returns:
        AnnData
    """

    sc.tl.umap(adata, min_dist=min_dist, spread=spread)
    return adata


def run_tsne(
    adata,
    n_pcs: int = 50,
):
    """
    Compute t-SNE embedding (optional alternative to UMAP).

    Parameters:
        adata (AnnData)
        n_pcs (int)

    Returns:
        AnnData
    """

    sc.tl.tsne(adata, n_pcs=n_pcs)
    return adata


def run_dimensionality_reduction(
    adata,
    n_comps: int = 50,
    n_neighbors: int = 15,
    method: str = "umap",
):
    """
    Full dimensionality pipeline:
    PCA → Neighbors → UMAP/t-SNE

    Parameters:
        adata (AnnData)
        n_comps (int)
        n_neighbors (int)
        method (str): 'umap' or 'tsne'

    Returns:
        AnnData
    """

    # Step 1: PCA
    adata = run_pca(adata, n_comps=n_comps)

    # Step 2: Neighbors
    adata = compute_neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps)

    # Step 3: Embedding
    if method == "umap":
        adata = run_umap(adata)
    elif method == "tsne":
        adata = run_tsne(adata)
    else:
        raise ValueError("Invalid method. Choose 'umap' or 'tsne'.")

    return adata


def plot_pca_variance(adata, n_pcs: int = 50, save: bool = False):
    """
    Plot PCA variance ratio.

    Helps decide number of PCs.

    Parameters:
        adata (AnnData)
        n_pcs (int)
        save (bool)

    Returns:
        None
    """

    sc.pl.pca_variance_ratio(
        adata,
        n_pcs=n_pcs,
        log=True,
        save="_pca_variance.png" if save else None,
        show=not save,
    )


def plot_embedding(
    adata,
    basis: str = "umap",
    color=None,
    save: bool = False,
):
    """
    Plot embedding (UMAP/t-SNE).

    Parameters:
        adata (AnnData)
        basis (str): 'umap', 'tsne', or 'pca'
        color (list/str): Annotation keys (e.g., clusters)
        save (bool)

    Returns:
        None
    """

    sc.pl.embedding(
        adata,
        basis=basis,
        color=color,
        save=f"_{basis}.png" if save else None,
        show=not save,
    )