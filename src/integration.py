# src/integration.py

import scanpy as sc
import anndata as ad
from typing import List


def load_datasets(file_paths: List[str]):
    """
    Load multiple datasets into a list of AnnData objects.
    Supports .h5ad and 10X formats.
    """
    adatas = []

    for path in file_paths:
        if path.endswith(".h5ad"):
            adata = sc.read_h5ad(path)
        else:
            # assume 10x format
            adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)

        adata.var_names_make_unique()
        adatas.append(adata)

    return adatas


def concatenate_datasets(adatas: List[ad.AnnData], batch_key="batch"):
    """
    Concatenate multiple AnnData objects and assign batch labels.
    """
    combined = ad.concat(adatas, label=batch_key, keys=[f"batch_{i}" for i in range(len(adatas))])
    return combined


def preprocess_for_integration(adata):
    """
    Basic preprocessing before integration.
    """
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")

    return adata


def run_harmony(adata, batch_key="batch"):
    """
    Run Harmony batch correction.
    Requires: harmonypy installed
    """
    try:
        import scanpy.external as sce
    except ImportError:
        raise ImportError("scanpy.external not found. Install harmonypy.")

    sce.pp.harmony_integrate(adata, key=batch_key)

    # Use corrected embeddings
    adata.obsm["X_pca"] = adata.obsm["X_pca_harmony"]

    return adata


def run_bbknn(adata, batch_key="batch"):
    """
    Run BBKNN batch correction.
    Requires: bbknn installed
    """
    try:
        import bbknn
    except ImportError:
        raise ImportError("bbknn not installed. pip install bbknn")

    bbknn.bbknn(adata, batch_key=batch_key)
    return adata


def run_integration_pipeline(file_paths: List[str],
                             method="harmony",
                             batch_key="batch"):
    """
    Full integration pipeline:
    Load → Concatenate → Preprocess → Integrate
    """

    print("📥 Loading datasets...")
    adatas = load_datasets(file_paths)

    print("🔗 Concatenating datasets...")
    adata = concatenate_datasets(adatas, batch_key=batch_key)

    print("⚙️ Preprocessing...")
    adata = preprocess_for_integration(adata)

    print(f"🧬 Running {method} integration...")

    if method.lower() == "harmony":
        adata = run_harmony(adata, batch_key=batch_key)
        sc.pp.neighbors(adata, use_rep="X_pca")

    elif method.lower() == "bbknn":
        adata = run_bbknn(adata, batch_key=batch_key)

    else:
        raise ValueError("Unsupported method. Choose 'harmony' or 'bbknn'")

    print("📉 Computing UMAP...")
    sc.tl.umap(adata)

    return adata