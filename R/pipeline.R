# Single-cell RNA-seq of 10x Genomics PBMCs with Seurat v5.

# Canonical PBMC markers (Seurat PBMC3k tutorial and common practice).
PBMC_MARKERS <- list(
  `Naive CD4 T` = c("IL7R", "CCR7", "LEF1"),
  `Memory CD4 T` = c("IL7R", "S100A4", "CD40LG"),
  `CD8 T` = c("CD8A", "CD8B", "CD3D"),
  NK = c("GNLY", "NKG7", "KLRF1"),
  B = c("MS4A1", "CD79A", "CD79B"),
  `CD14 Mono` = c("CD14", "LYZ", "S100A8"),
  `FCGR3A Mono` = c("FCGR3A", "MS4A7", "LST1"),
  DC = c("FCER1A", "CST3", "CLEC10A"),
  Platelet = c("PPBP", "PF4", "GP9")
)

qc_filter <- function(seu, min_features = 200, max_features = 2500, max_mt = 5) {
  seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
  keep <- seu$nFeature_RNA > min_features & seu$nFeature_RNA < max_features & seu$percent.mt < max_mt
  seu[, keep]
}

preprocess <- function(seu, n_pcs = 30, n_hvg = 2000) {
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = n_hvg, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  Seurat::RunPCA(seu, npcs = n_pcs, verbose = FALSE)
}

# 10x Genomics guidance: roughly 0.8% multiplets per 1,000 cells recovered
# (v3 chemistry). The original assumed 7.5% for ~2,700 cells, about three
# times too many, which removes real singlets.
expected_doublet_rate <- function(n_cells, per_thousand = 0.008) per_thousand * n_cells / 1000

remove_doublets <- function(seu, dims, clusters, rate = NULL, seed = 42) {
  set.seed(seed)
  rate <- if (is.null(rate)) expected_doublet_rate(ncol(seu)) else rate
  sweep <- suppressMessages(DoubletFinder::paramSweep(seu, PCs = dims, sct = FALSE))
  stats <- DoubletFinder::summarizeSweep(sweep, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(stats)
  pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  homotypic <- DoubletFinder::modelHomotypic(clusters)       # homotypic doublets are undetectable
  n_exp <- round(rate * ncol(seu) * (1 - homotypic))
  seu <- suppressMessages(DoubletFinder::doubletFinder(seu, PCs = dims, pN = 0.25, pK = pk, nExp = n_exp, sct = FALSE))
  # Find the classification column by pattern; the original hard-coded
  # "DF.classifications_0.25_0.01_201", which breaks whenever pK or nExp change.
  col <- tail(grep("^DF.classifications", colnames(seu[[]]), value = TRUE), 1)
  seu$doublet_call <- seu[[col]][, 1]
  attr(seu, "doublet_params") <- list(pK = pk, nExp = n_exp, rate = rate, homotypic = homotypic)
  seu
}

cluster_cells <- function(seu, dims, resolution = 0.8, seed = 42) {
  seu <- Seurat::FindNeighbors(seu, dims = dims, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, resolution = resolution, random.seed = seed, verbose = FALSE)
  Seurat::RunUMAP(seu, dims = dims, seed.use = seed, verbose = FALSE)
}

# Cluster-level annotation from marker module scores. Each cluster gets the
# cell type whose markers score highest on average; the margin to the
# runner-up is reported so ambiguous clusters are visible.
annotate_by_markers <- function(seu, markers = PBMC_MARKERS) {
  markers <- lapply(markers, intersect, rownames(seu))
  markers <- markers[lengths(markers) > 0]
  ctrl <- min(100, max(5, floor(nrow(seu) / 30)))   # control genes per expression bin (small panels need fewer)
  seu <- Seurat::AddModuleScore(seu, features = markers, name = "ms_", ctrl = ctrl, seed = 1)
  score_cols <- paste0("ms_", seq_along(markers))
  scores <- as.matrix(seu[[score_cols]]); colnames(scores) <- names(markers)
  by_cluster <- apply(scores, 2, function(s) tapply(s, seu$seurat_clusters, mean))
  if (is.null(dim(by_cluster))) by_cluster <- matrix(by_cluster, nrow = 1, dimnames = list(levels(seu$seurat_clusters), names(markers)))
  top <- apply(by_cluster, 1, function(r) names(r)[which.max(r)])
  margin <- apply(by_cluster, 1, function(r) { s <- sort(r, decreasing = TRUE); s[1] - s[2] })
  seu$celltype <- unname(top[as.character(seu$seurat_clusters)])
  seu[[score_cols]] <- NULL
  attr(seu, "annotation") <- data.frame(cluster = rownames(by_cluster), celltype = top, margin = round(margin, 3),
                                        n_cells = as.integer(table(seu$seurat_clusters)[rownames(by_cluster)]))
  seu
}

# Optional reference-based annotation. Pass a SummarizedExperiment reference
# saved to disk, e.g. celldex::MonacoImmuneData() (immune-specific; better
# suited to PBMCs than the HumanPrimaryCellAtlas used before).
annotate_singler <- function(seu, ref, label = "label.main") {
  pred <- SingleR::SingleR(test = Seurat::GetAssayData(seu, layer = "data"), ref = ref, labels = ref[[label]],
                           clusters = seu$seurat_clusters)
  seu$celltype_singler <- pred$pruned.labels[match(as.character(seu$seurat_clusters), rownames(pred))]
  seu
}

# Subclustering re-selects variable genes and recomputes PCA/UMAP WITHIN the
# subset. The original reused the parent PCA/UMAP and picked subsets by
# hard-coded cluster numbers ("cluster 0 = T cells"), which change between runs.
subcluster <- function(seu, celltypes, dims = 1:10, resolution = 0.5, seed = 42) {
  sub <- seu[, seu$celltype %in% celltypes]
  if (ncol(sub) < 50) stop("too few cells to subcluster: ", ncol(sub))
  sub <- Seurat::FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
  sub <- Seurat::ScaleData(sub, verbose = FALSE)
  sub <- Seurat::RunPCA(sub, npcs = max(dims), verbose = FALSE)
  cluster_cells(sub, dims, resolution, seed)
}
