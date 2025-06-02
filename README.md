# RNA-Seq-using-filteredData
📘 Overview
This repository contains an R script for single-cell RNA-seq (scRNA-seq) analysis using the 3k Peripheral Blood Mononuclear Cell (PBMC) dataset from a healthy donor. The dataset is publicly available via Figshare.

The analysis pipeline is built using the Seurat framework and performs:

Quality control

Normalization

Dimensionality reduction

Clustering

Cell type annotation

Differential expression analysis

The workflow is optimized to run on personal computers with modest specifications (~8–16 GB RAM, 4–8 CPU cores).

📂 Dataset

The dataset consists of ~2,700 PBMCs, processed using 10x Genomics' Cell Ranger pipeline. It includes filtered feature-barcode matrices and can be downloaded from:

👉 https://figshare.com/articles/dataset/3k_PBMCs_from_a_healthy_donor/28414916

Expected directory structure after download:

RNA-Seq-using-filteredData/

├── 28414916/

│   └── filtered_feature_bc_matrix/

├── pbmc_analysis.R

└── README.md

🔧 Prerequisites

📦 R Environment
R version 4.0 or higher

📚 Required R Packages

Install the required packages using the commands below:

install.packages(c("Seurat", "ggplot2", "patchwork"))
install.packages("BiocManager")
BiocManager::install(c("DoubletFinder", "SingleR", "celldex"))
🚀 Installation & Setup
Clone the repository:


git clone https://github.com/your-username/RNA-Seq-using-filteredData.git

cd RNA-Seq-using-filteredData

Download the dataset:

Go to Figshare

Download and extract the filtered_feature_bc_matrix folder

Place it inside the 28414916/ directory in the repo

▶️ Usage
Open pbmc_analysis.R in R or RStudio

Update the path to the dataset in the Read10X function if needed:


pbmc.data <- Read10X(data.dir = "28414916/filtered_feature_bc_matrix")
Run the script:


source("pbmc_analysis.R")

📈 Output
The script will generate:

UMAP plots for visualizing cell clusters (umap_clusters.png, umap_celltypes.png)

Cluster annotation results

Differential expression results

🧬 Tools Used

Seurat

DoubletFinder

SingleR

celldex

📫 Contact

For questions or suggestions, feel free to open an issue or reach out via GitHub.
