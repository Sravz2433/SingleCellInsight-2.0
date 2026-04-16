# src/pathway.py

import os
import pandas as pd
import gseapy as gp


def run_enrichment_for_cluster(gene_list, cluster_name, outdir,
                               organism="Human",
                               gene_sets=None):
    """
    Run enrichment analysis for a single cluster.
    """

    if gene_sets is None:
        gene_sets = [
            "GO_Biological_Process_2021",
            "KEGG_2021_Human"
        ]

    cluster_outdir = os.path.join(outdir, f"cluster_{cluster_name}")
    os.makedirs(cluster_outdir, exist_ok=True)

    print(f"🧬 Running enrichment for cluster {cluster_name}...")

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism,
            outdir=cluster_outdir,
            cutoff=0.05  # adjusted p-value cutoff
        )

        return enr

    except Exception as e:
        print(f"❌ Enrichment failed for cluster {cluster_name}: {e}")
        return None


def extract_top_genes(marker_df, cluster, top_n=50):
    """
    Extract top marker genes for a given cluster.
    Assumes Scanpy rank_genes_groups output converted to DataFrame.
    """

    df = marker_df[marker_df["cluster"] == cluster]
    # Determine sorting column dynamically
    if "logfoldchanges" in df.columns:
        sort_col = "logfoldchanges"
    elif "scores" in df.columns:
        sort_col = "scores"
    elif "log2fc" in df.columns:
        sort_col = "log2fc"
    else:
        raise ValueError("No valid ranking column found in marker dataframe")

    df = df.sort_values(by=sort_col, ascending=False)

    top_genes = df["gene"].head(top_n).tolist()
    return top_genes


def run_pathway_analysis(marker_df,
                         output_dir="results/pathways",
                         top_n=50):
    """
    Run pathway enrichment for all clusters.
    """

    os.makedirs(output_dir, exist_ok=True)

    clusters = marker_df["cluster"].unique()

    all_results = []

    for cluster in clusters:
        genes = extract_top_genes(marker_df, cluster, top_n=top_n)

        enr = run_enrichment_for_cluster(
            gene_list=genes,
            cluster_name=cluster,
            outdir=output_dir
        )

        if enr is not None:
            # Collect results
            for gs in enr.results:
                df = enr.results[gs]
                df["cluster"] = cluster
                df["gene_set"] = gs
                all_results.append(df)

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(os.path.join(output_dir, "combined_pathway_results.csv"), index=False)
        print("✅ Pathway analysis completed and saved.")

        return combined
    else:
        print("⚠️ No pathway results generated.")
        return None