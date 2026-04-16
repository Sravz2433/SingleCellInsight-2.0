# src/report.py

import os
from datetime import datetime


def generate_report(output_dir,
                    figures_dir="results/figures",
                    markers_file=None,
                    pathway_file=None,
                    title="SingleCellInsight Report"):
    """
    Generate an HTML report summarizing results.
    """

    report_path = os.path.join(output_dir, "report.html")
    os.makedirs(output_dir, exist_ok=True)

    def get_image_tag(img_path):
        if os.path.exists(img_path):
            return f'<img src="{img_path}" width="600">'
        return "<p>Image not found</p>"

    # Define expected figures
    umap_img = os.path.join(figures_dir, "umap.png")
    qc_img = os.path.join(figures_dir, "qc.png")
    marker_img = os.path.join(figures_dir, "markers.png")

    # Load marker preview
    marker_section = "<p>No marker data available</p>"
    if markers_file and os.path.exists(markers_file):
        import pandas as pd
        df = pd.read_csv(markers_file)
        marker_section = df.head(10).to_html(index=False)

    # Load pathway preview
    pathway_section = "<p>No pathway data available</p>"
    if pathway_file and os.path.exists(pathway_file):
        import pandas as pd
        df = pd.read_csv(pathway_file)
        pathway_section = df.head(10).to_html(index=False)

    # HTML Template
    html_content = f"""
    <html>
    <head>
        <title>{title}</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 40px;
                background-color: #f8f9fa;
            }}
            h1, h2 {{
                color: #2c3e50;
            }}
            .section {{
                background: white;
                padding: 20px;
                margin-bottom: 20px;
                border-radius: 10px;
                box-shadow: 0px 2px 5px rgba(0,0,0,0.1);
            }}
            table {{
                border-collapse: collapse;
                width: 100%;
            }}
            table, th, td {{
                border: 1px solid #ddd;
            }}
            th, td {{
                padding: 8px;
                text-align: left;
            }}
        </style>
    </head>
    <body>

    <h1>{title}</h1>
    <p><b>Generated on:</b> {datetime.now().strftime("%Y-%m-%d %H:%M")}</p>

    <div class="section">
        <h2>🧬 Pipeline Overview</h2>
        <p>
        Data → QC → Normalization → PCA → UMAP → Clustering → Annotation → Marker Detection → Pathway Analysis
        </p>
    </div>

    <div class="section">
        <h2>📊 Quality Control</h2>
        {get_image_tag(qc_img)}
    </div>

    <div class="section">
        <h2>📉 UMAP Visualization</h2>
        {get_image_tag(umap_img)}
    </div>

    <div class="section">
        <h2>🧠 Marker Genes (Top)</h2>
        {marker_section}
    </div>

    <div class="section">
        <h2>🧪 Pathway Enrichment (Top)</h2>
        {pathway_section}
    </div>

    <div class="section">
        <h2>📌 Summary</h2>
        <p>
        This analysis identifies distinct cell populations and their associated marker genes.
        Pathway enrichment suggests functional roles of these clusters, enabling biological interpretation.
        </p>
    </div>

    </body>
    </html>
    """

    with open(report_path, "w") as f:
        f.write(html_content)

    print(f"📄 Report generated at: {report_path}")

    return report_path