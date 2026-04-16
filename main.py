# main.py

import argparse
import yaml
import scanpy as sc

from src.preprocessing import run_qc_pipeline, plot_qc
from src.dimensionality import run_dimensionality_reduction, plot_pca_variance, plot_embedding
from src.clustering import run_clustering, plot_clusters
from src.markers import run_marker_pipeline, save_markers, plot_top_markers
from src.annotation import run_celltypist_annotation, plot_cell_types
from src.integration import run_integration_pipeline
from src.pathway import run_pathway_analysis
from src.report import generate_report

def load_config(config_path):
    """Load YAML config file"""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def load_data(input_path):
    """
    Load input dataset (10x or h5ad)
    """

    import os

    # Case 1: Preprocessed AnnData
    if input_path.endswith(".h5ad"):
        return sc.read_h5ad(input_path)

    # Case 2: 10x Genomics folder
    elif os.path.isdir(input_path):
        print("📂 Detected 10x Genomics format...")
        return sc.read_10x_mtx(
            input_path,
            var_names="gene_symbols",  # use gene symbols
            cache=True
        )

    else:
        raise ValueError(
            "Unsupported format. Provide a folder (10x) or .h5ad file"
        )

def main(config):
    # ========================
    # Load Data
    # ========================
    print("📥 Loading data...")
    if config.get("integration", {}).get("enabled", False):
        print("🔗 Running multi-sample integration...")
        adata = run_integration_pipeline(
            file_paths=config["integration"]["files"],
            method=config["integration"]["method"],
            batch_key=config["integration"].get("batch_key", "batch")
        )
    else:
        print("📥 Loading single dataset...")
        adata = load_data(config["input"]["path"])

    # ========================
    # Preprocessing
    # ========================
    print("🧪 Running preprocessing...")
    adata = run_qc_pipeline(
        adata,
        mito_prefix=config["preprocessing"]["mito_prefix"],
        min_genes=config["preprocessing"]["min_genes"],
        min_cells=config["preprocessing"]["min_cells"],
        max_mt_pct=config["preprocessing"]["max_mt_pct"],
        n_top_genes=config["preprocessing"]["n_top_genes"],
    )

    if config["plots"]["qc"]:
        plot_qc(adata, save=True)

    # ========================
    # Dimensionality Reduction
    # ========================
    print("📉 Running dimensionality reduction...")
    adata = run_dimensionality_reduction(
        adata,
        n_comps=config["dimensionality"]["n_comps"],
        n_neighbors=config["dimensionality"]["n_neighbors"],
        method=config["dimensionality"]["method"],
    )

    if config["plots"]["pca"]:
        plot_pca_variance(adata, save=True)

    if config["plots"]["embedding"]:
        plot_embedding(adata, basis=config["dimensionality"]["method"], save=True)

    # ========================
    # Clustering
    # ========================
    print("🔍 Running clustering...")
    adata = run_clustering(
        adata,
        method=config["clustering"]["method"],
        resolution=config["clustering"]["resolution"],
        n_neighbors=config["clustering"]["n_neighbors"],
        n_pcs=config["clustering"]["n_pcs"],
    )

    if config["plots"]["clusters"]:
        plot_clusters(adata, method=config["clustering"]["method"], save=True)
    
    # ========================
    # Cell Type Annotation
    # ========================
    if config["annotation"]["enabled"]:
        print("🧠 Annotating cell types...")

        adata = run_celltypist_annotation(
            adata,
            model=config["annotation"]["model"]
        )

        if config["plots"]["cell_types"]:
            plot_cell_types(adata, save=True)
    
    # ========================
    # Marker Genes
    # ========================
    print("🧬 Detecting marker genes...")
    markers_df = run_marker_pipeline(
        adata,
        groupby=config["clustering"]["method"],
        method=config["markers"]["method"],
        logfc_min=config["markers"]["logfc_min"],
        pval_adj_max=config["markers"]["pval_adj_max"],
    )

    save_markers(markers_df, filepath=config["output"]["markers_file"])

    if config["plots"]["markers"]:
        plot_top_markers(adata, save=True)
    # ========================
    # Pathway Analysis
    # ========================
    if config.get("pathway", {}).get("enabled", False):
        print("🧪 Running pathway enrichment...")

        pathway_results = run_pathway_analysis(
            markers_df,
            output_dir=config["pathway"]["output_dir"],
            top_n=config["pathway"]["top_n"]
        )
    
    # ========================
    # Report Generation
    # ========================
    if config.get("report", {}).get("enabled", False):
        print("📄 Generating report...")

        generate_report(
            output_dir=config["report"]["output_dir"],
            figures_dir="results/figures",
            markers_file=config["output"]["markers_file"],
            pathway_file=f"{config['pathway']['output_dir']}/combined_pathway_results.csv"
            if config.get("pathway", {}).get("enabled", False) else None
        )
        
    # ========================
    # Save Results
    # ========================
    print("💾 Saving processed data...")
    adata.write(config["output"]["adata_file"])

    print("✅ Pipeline completed successfully!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SingleCellInsight Pipeline")
    parser.add_argument("--config", type=str, required=True, help="Path to config.yaml")

    args = parser.parse_args()

    config = load_config(args.config)
    main(config)

