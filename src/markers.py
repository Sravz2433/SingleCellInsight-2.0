# src/markers.py

import scanpy as sc
import pandas as pd


def run_marker_detection(
    adata,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    n_genes: int = 100,
    use_raw: bool = True,
):
    """
    Identify marker genes using differential expression.

    Parameters:
        adata (AnnData)
        groupby (str): Cluster key in adata.obs
        method (str): 'wilcoxon', 't-test', 'logreg'
        n_genes (int): Number of genes to rank
        use_raw (bool): Whether to use raw counts

    Returns:
        AnnData
    """

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        n_genes=n_genes,
        use_raw=use_raw,
    )

    return adata


def get_marker_table(
    adata,
    group: str = None,
    key: str = "rank_genes_groups",
):
    """
    Extract marker genes into a DataFrame.

    Parameters:
        adata (AnnData)
        group (str): Specific cluster (optional)
        key (str): Key from rank_genes_groups

    Returns:
        pd.DataFrame
    """

    result = adata.uns[key]

    groups = result["names"].dtype.names
    dfs = []

    for g in groups:
        df = pd.DataFrame({
            "gene": result["names"][g],
            "logfoldchange": result["logfoldchanges"][g],
            "pvals": result["pvals"][g],
            "pvals_adj": result["pvals_adj"][g],
            "scores": result["scores"][g],
            "cluster": g
        })
        dfs.append(df)

    markers_df = pd.concat(dfs)

    if group is not None:
        markers_df = markers_df[markers_df["cluster"] == group]

    return markers_df


def filter_markers(
    markers_df,
    logfc_min: float = 0.25,
    pval_adj_max: float = 0.05,
):
    """
    Filter marker genes based on thresholds.

    Parameters:
        markers_df (pd.DataFrame)
        logfc_min (float)
        pval_adj_max (float)

    Returns:
        pd.DataFrame
    """

    filtered = markers_df[
        (markers_df["logfoldchange"] > logfc_min) &
        (markers_df["pvals_adj"] < pval_adj_max)
    ]

    return filtered


def save_markers(
    markers_df,
    filepath: str = "results/markers.csv",
):
    """
    Save marker genes to CSV.

    Parameters:
        markers_df (pd.DataFrame)
        filepath (str)

    Returns:
        None
    """

    markers_df.to_csv(filepath, index=False)


def plot_top_markers(
    adata,
    groupby: str = "leiden",
    n_genes: int = 20,
    save: bool = False,
):
    """
    Plot top marker genes.

    Parameters:
        adata (AnnData)
        groupby (str)
        n_genes (int)
        save (bool)

    Returns:
        None
    """

    sc.pl.rank_genes_groups(
        adata,
        n_genes=n_genes,
        sharey=False,
        save="_top_markers.png" if save else None,
        show=not save,
    )


def plot_heatmap(
    adata,
    groupby: str = "leiden",
    n_genes: int = 20,
    save: bool = False,
):
    """
    Plot heatmap of marker genes.

    Parameters:
        adata (AnnData)
        groupby (str)
        n_genes (int)
        save (bool)

    Returns:
        None
    """

    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=n_genes,
        groupby=groupby,
        show=not save,
        save="_heatmap.png" if save else None,
    )


def run_marker_pipeline(
    adata,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    logfc_min: float = 0.25,
    pval_adj_max: float = 0.05,
):
    """
    Full marker workflow:
    DE → extract → filter

    Returns:
        pd.DataFrame (filtered markers)
    """

    # Step 1: Differential expression
    adata = run_marker_detection(adata, groupby=groupby, method=method)

    # Step 2: Extract markers
    markers_df = get_marker_table(adata)

    # Step 3: Filter markers
    markers_df = filter_markers(
        markers_df,
        logfc_min=logfc_min,
        pval_adj_max=pval_adj_max,
    )

    return markers_df