# 🧬 SingleCellInsight

### Automated scRNA-seq Analysis Toolkit

**SingleCellInsight** is a modular, reproducible, and Dockerized pipeline for single-cell RNA-seq (scRNA-seq) analysis built using **Scanpy**.
It automates the complete workflow—from quality control to clustering, cell type annotation, and marker gene detection—enabling rapid and consistent exploratory analysis of single-cell datasets.

---

## 🚀 Key Features

* 📊 **Automated Quality Control (QC)**

  * Filtering low-quality cells and genes
  * Mitochondrial gene percentage analysis

* 📉 **Dimensionality Reduction**

  * Principal Component Analysis (PCA)
  * UMAP visualization

* 🔍 **Clustering**

  * Leiden clustering for identifying cell populations

* 🧠 **Cell Type Annotation**

  * Integrated with CellTypist for automated biological interpretation

* 🧬 **Marker Gene Detection**

  * Differential gene expression analysis for clusters

* 🔁 **Reproducibility**

  * Config-driven pipeline (`config.yaml`)
  * Fully Dockerized execution

* 🧩 **Modular Design**

  * Easily extendable for custom workflows

---

## 🧠 Pipeline Workflow

Data → QC → Normalization → PCA → UMAP → Clustering → Annotation → Marker Detection

---

## 📊 Results

### 🔬 UMAP Clustering

![UMAP](results/figures/umap_umap.png)

### 🧠 Cell Type Annotation

![Cell Types](results/figures/umap_cell_types.png)

### 🧬 Marker Genes

![Markers](results/figures/top_markers.png)

---

## 🏗️ Project Structure

```
SingleCellInsight/
│── src/
│   ├── preprocessing.py
│   ├── dimensionality.py
│   ├── clustering.py
│   ├── markers.py
│   ├── annotation.py
│── results/
│   └── figures/
│── Dockerfile
│── config.yaml
│── main.py
│── requirements.txt
│── README.md
```

---

## ⚙️ Installation

```bash
git clone https://github.com/Sravz2433/SingleCellInsight.git
cd SingleCellInsight
pip install -r requirements.txt
```

---

## ▶️ Usage

### Run full pipeline

```bash
python main.py --config config.yaml
```

---

## 🐳 Docker (Recommended)

```bash
docker build -t singlecellinsight .
docker run -v E:/SingleCellInsight/data:/app/data -v E:/SingleCellInsight/results:/app/results singlecellinsight
```

---

## 📥 Input Data

Supports:

* 10x Genomics format (`mtx`, `barcodes`, `features`)
* `.h5ad` files

---

## 📤 Output

* Processed dataset (`.h5ad`)
* Marker genes (`.csv`)
* QC & visualization plots

---

## 🔬 Technologies

* Python
* Scanpy
* AnnData
* CellTypist
* NumPy / Pandas
* Matplotlib / Seaborn
* Docker

---

## 🎯 Motivation

Single-cell RNA-seq analysis involves multiple complex steps and tools.
This project simplifies and standardizes the workflow into a single reproducible pipeline, enabling faster and more reliable biological insights.

---

## 🚀 Future Work

* Batch correction (BBKNN / Harmony)
* Marker-based validation
* Streamlit-based visualization UI
* Multi-dataset support

---

## 📄 License

MIT License

---

## 👩‍💻 Author

**Sravya Sri Mallampalli**
B.Tech CSE (Bioinformatics)
Aspiring Computational Biologist
