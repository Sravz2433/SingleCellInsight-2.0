# src/annotation.py

import scanpy as sc
import celltypist


def run_celltypist_annotation(
    adata,
    model: str = "Immune_All_Low.pkl",
    majority_voting: bool = True,
):
    """
    Annotate cell types using CellTypist.

    Parameters:
        adata (AnnData)
        model (str): Pretrained model name
        majority_voting (bool)

    Returns:
        AnnData
    """

    # Ensure data is normalized/log-transformed
    if adata.raw is not None:
        input_data = adata.raw.to_adata()
    else:
        input_data = adata

    print("🔎 Running CellTypist annotation...")

    celltypist.models.download_models(force_update=False)
    
    predictions = celltypist.annotate(
        input_data,
        model=model,
        majority_voting=majority_voting,
    )

    # Add predictions to AnnData
    if hasattr(predictions.predicted_labels, "columns"):
        adata.obs["cell_type"] = predictions.predicted_labels.iloc[:, 0]
    else:
        adata.obs["cell_type"] = predictions.predicted_labels

    return adata


def plot_cell_types(adata, save: bool = False):
    """
    Plot UMAP with cell type annotations.
    """

    sc.pl.umap(
        adata,
        color="cell_type",
        save="_cell_types.png" if save else None,
        show=not save,
    )