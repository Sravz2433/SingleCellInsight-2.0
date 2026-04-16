import scanpy as sc


def run_neighbors(
    adata,
    n_neighbors: int = 15,
    n_pcs: int = 50,
    metric: str = "euclidean",
):
    """
    Compute neighborhood graph.

    Parameters:
        adata (AnnData): Annotated data matrix
        n_neighbors (int): Number of neighbors
        n_pcs (int): Number of principal components
        metric (str): Distance metric

    Returns:
        AnnData: Updated object with neighbors computed
    """
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, metric=metric)
    return adata


def run_leiden(
    adata,
    resolution: float = 1.0,
    key_added: str = "leiden",
):
    """
    Perform Leiden clustering.

    Parameters:
        adata (AnnData): Annotated data matrix
        resolution (float): Clustering resolution
        key_added (str): Key to store results in adata.obs

    Returns:
        AnnData: Updated object with Leiden clusters
    """
    sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
    return adata


def run_louvain(
    adata,
    resolution: float = 1.0,
    key_added: str = "louvain",
):
    """
    Perform Louvain clustering.

    Parameters:
        adata (AnnData): Annotated data matrix
        resolution (float): Clustering resolution
        key_added (str): Key to store results

    Returns:
        AnnData: Updated object with Louvain clusters
    """
    sc.tl.louvain(adata, resolution=resolution, key_added=key_added)
    return adata


def run_clustering(
    adata,
    method: str = "leiden",
    resolution: float = 1.0,
    n_neighbors: int = 15,
    n_pcs: int = 50,
):
    """
    Full clustering pipeline:
    neighbors → clustering

    Parameters:
        adata (AnnData): Annotated data matrix
        method (str): 'leiden' or 'louvain'
        resolution (float): Clustering resolution
        n_neighbors (int): Number of neighbors
        n_pcs (int): Number of PCs

    Returns:
        AnnData: Updated object with clustering results
    """

    # Step 1: Compute neighbors
    adata = run_neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    # Step 2: Run clustering
    if method == "leiden":
        adata = run_leiden(adata, resolution=resolution)
    elif method == "louvain":
        adata = run_louvain(adata, resolution=resolution)
    else:
        raise ValueError("Invalid method. Choose 'leiden' or 'louvain'.")

    return adata


def plot_clusters(
    adata,
    method: str = "leiden",
    basis: str = "umap",
    save: bool = False,
):
    """
    Plot clustering results.

    Parameters:
        adata (AnnData): Annotated data matrix
        method (str): Clustering label key
        basis (str): Embedding ('umap' or 'pca')
        save (bool): Save figure

    Returns:
        None
    """

    sc.pl.embedding(
        adata,
        basis=basis,
        color=[method],
        save="_clusters.png" if save else None,
        show=not save,
    )