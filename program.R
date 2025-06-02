# Load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(SingleR)
library(celldex)

# 1. Load and Preprocess PBMC Dataset
pbmc.data <- Read10X(data.dir = "28414916/filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "PBMC3k", min.cells = 3, min.features = 200)
print(dim(pbmc))  # Initial dimensions: ~2,700 cells, ~20,000-30,000 genes

# 2. Quality Control (QC)
# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("qc_violin_plots.png", width = 10, height = 4)

# Filter cells based on QC thresholds
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)
print(dim(pbmc))  # Dimensions after QC filtering

# 3. Preprocessing for Dimensionality Reduction and Doublet Detection
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs = 30)

# Visualize PCA elbow plot to assess dimensionality
ElbowPlot(pbmc, ndims = 30)
ggsave("elbow_plot.png", width = 6, height = 4)

# Run UMAP with 15 PCs (adjusted from original 10 based on later code)
pbmc <- RunUMAP(pbmc, dims = 1:15)

# 4. Doublet Detection with DoubletFinder
ls("package:DoubletFinder")  # List DoubletFinder functions
sweep.res <- paramSweep(pbmc, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print(bcmvn)
pK_opt <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
print(pK_opt)
nExp <- round(0.075 * ncol(pbmc))  # Expected doublets: ~7.5% of cells
pbmc <- doubletFinder(pbmc, PCs = 1:10, pN = 0.25, pK = pK_opt, nExp = nExp, sct = FALSE)
pbmc <- subset(pbmc, subset = DF.classifications_0.25_0.01_201 == "Singlet")
print(dim(pbmc))  # Dimensions after doublet removal

# 5. Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.8)  # Higher resolution for refined clusters
DimPlot(pbmc, reduction = "umap", label = TRUE)
ggsave("umap_clusters.png", width = 6, height = 6)

# View cluster sizes
table(pbmc$seurat_clusters)

# 6. Cell Type Annotation with SingleR
ref <- HumanPrimaryCellAtlasData()
pbmc.counts <- GetAssayData(pbmc, slot = "counts")
pred <- SingleR(test = pbmc.counts, ref = ref, labels = ref$label.main)
pbmc$celltype <- pred$labels
DimPlot(pbmc, reduction = "umap", group.by = "celltype", label = TRUE)
ggsave("umap_celltypes.png", width = 6, height = 6)

# 7. Marker Gene Identification
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers, 5)
write.csv(markers, "cluster_markers.csv")

# Differential expression between clusters 0 and 4
de_genes <- FindMarkers(pbmc, ident.1 = 0, ident.2 = 4, min.pct = 0.25)
head(de_genes, 5)
write.csv(de_genes, "de_cluster0_vs_cluster4.csv")

# Visualize key marker genes
FeaturePlot(pbmc, features = c("CD3D", "CD14", "CD79A", "NKG7", "PPBP", "CD4", "CD8A", "FCGR3A"), reduction = "umap")
ggsave("feature_plots.png", width = 10, height = 8)

# 8. Subclustering T-cells (Cluster 0)
t_cells <- subset(pbmc, idents = 0)
t_cells <- FindNeighbors(t_cells, dims = 1:10)
t_cells <- FindClusters(t_cells, resolution = 0.5)
DimPlot(t_cells, reduction = "umap", label = TRUE)
ggsave("umap_tcell_subclusters.png", width = 6, height = 6)

# Identify T-cell subcluster markers
tcell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(tcell_markers, 5)
write.csv(tcell_markers, "tcell_subcluster_markers.csv")

# 9. Subclustering Monocytes (Cluster 2)
monocytes <- subset(pbmc, idents = 2)
monocytes <- FindNeighbors(monocytes, dims = 1:10)
monocytes <- FindClusters(monocytes, resolution = 0.5)
DimPlot(monocytes, reduction = "umap", label = TRUE)
ggsave("umap_monocyte_subclusters.png", width = 6, height = 6)

# 10. Save Processed Data
saveRDS(pbmc, file = "pbmc_processed.rds")

