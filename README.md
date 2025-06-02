# RNA-Seq-using-filteredData
Overview
This repository contains an R script for single-cell RNA-seq (scRNA-seq) analysis of the 3k PBMC dataset from a healthy donor, sourced from Figshare (https://figshare.com/articles/dataset/3k_PBMCs_from_a_healthy_donor/28414916). The pipeline uses Seurat to process ~2,700 cells, performing quality control, clustering, cell type annotation, and differential expression analysis. The code is optimized to run on a personal computer (PC) with limited resources (~8-16 GB RAM, 4-8 cores).

Dataset
The dataset is a preprocessed 10x Genomics dataset of Peripheral Blood Mononuclear Cells (PBMCs) with ~2,700 cells, available at the link above. It includes filtered feature-barcode matrices (filtered_feature_bc_matrix).

Prerequisites
R: Version 4.0 or higher

R Packages:

Seurat (>= 4.0)

ggplot2

patchwork

DoubletFinder

SingleR

celldex

Install the packages using:

R


install.packages(c("Seurat", "ggplot2", "patchwork"))
install.packages("BiocManager")
BiocManager::install(c("DoubletFinder", "SingleR", "celldex"))
Installation
Clone this repository to your local machine:

bash


git clone https://github.com/<your-username>/<your-repo-name>.git
cd <your-repo-name>
Download the 3k PBMC dataset from Figshare and place the filtered_feature_bc_matrix folder in the 28414916/ directory within the repository.

Usage
Open the R script pbmc_analysis.R in R or RStudio.

Update the path in the Read10X function to point to the filtered_feature_bc_matrix directory:

R


pbmc.data <- Read10X(data.dir = "28414916/filtered_feature_bc_matrix")
Run the script:

R


source("pbmc_analysis.R")
The script will generate:

UMAP plots (umap_clusters.png, `umap_cell
