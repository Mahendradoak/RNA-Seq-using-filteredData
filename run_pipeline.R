#!/usr/bin/env Rscript
# Usage:
#   Rscript run_pipeline.R --tenx 28414916/filtered_feature_bc_matrix --out results
#   Rscript run_pipeline.R --rds pbmc_processed.rds --out results_from_rds   (cells already QC'd/doublet-filtered)
# Options: --dims 15  --resolution 0.8  --singler-ref monaco.rds  --seed 42

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
if (!exists("remove_doublets", mode = "function")) {
  here <- local({
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(f)) dirname(normalizePath(f)) else getwd()
  })
  source(file.path(here, "R", "pipeline.R"))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  a <- as.list(setNames(args[seq(2, length(args), 2)], sub("^--", "", args[seq(1, length(args), 2)])))
  if (is.null(a$out) || (is.null(a$tenx) && is.null(a$rds))) stop("--out and one of --tenx / --rds are required")
  dims <- seq_len(as.integer(if (is.null(a$dims)) 15 else a$dims))
  res <- as.numeric(if (is.null(a$resolution)) 0.8 else a$resolution)
  seed <- as.integer(if (is.null(a$seed)) 42 else a$seed)
  dir.create(a$out, showWarnings = FALSE, recursive = TRUE)
  log <- list()

  if (!is.null(a$tenx)) {
    seu <- CreateSeuratObject(Read10X(a$tenx), min.cells = 3, min.features = 200)
    log$cells_raw <- ncol(seu)
    seu <- qc_filter(seu)
    log$cells_after_qc <- ncol(seu)
    ggsave(file.path(a$out, "qc_violin.png"), VlnPlot(seu, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0), width = 10, height = 4)
    seu <- cluster_cells(preprocess(seu), dims, res, seed)
    seu <- remove_doublets(seu, dims, seu$seurat_clusters, seed = seed)       # same PCs as clustering
    log$doublets <- attr(seu, "doublet_params"); log$doublets$n_called <- sum(seu$doublet_call == "Doublet")
    seu <- seu[, seu$doublet_call == "Singlet"]
  } else {
    old <- readRDS(a$rds)
    seu <- CreateSeuratObject(GetAssayData(old, layer = "counts"))
    log$note <- "started from saved counts; QC and doublet removal were applied upstream"
  }
  log$cells_analysed <- ncol(seu)

  seu <- cluster_cells(preprocess(seu), dims, res, seed)
  ggsave(file.path(a$out, "elbow.png"), ElbowPlot(seu, ndims = 30), width = 6, height = 4)
  seu <- annotate_by_markers(seu)
  ann <- attr(seu, "annotation")
  utils::write.csv(ann, file.path(a$out, "cluster_annotation.csv"), row.names = FALSE)
  if (!is.null(a$`singler-ref`)) seu <- annotate_singler(seu, readRDS(a$`singler-ref`))
  ggsave(file.path(a$out, "umap_clusters.png"), DimPlot(seu, label = TRUE) + NoLegend(), width = 6, height = 6)
  ggsave(file.path(a$out, "umap_celltypes.png"), DimPlot(seu, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend(), width = 7, height = 6)
  ggsave(file.path(a$out, "marker_dotplot.png"),
         DotPlot(seu, features = unique(unlist(PBMC_MARKERS)), group.by = "celltype") + RotatedAxis(), width = 12, height = 5)

  markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  utils::write.csv(markers, file.path(a$out, "cluster_markers.csv"), row.names = FALSE)

  # A biologically defined contrast instead of "cluster 0 vs cluster 4".
  Idents(seu) <- "celltype"
  if (all(c("CD14 Mono", "FCGR3A Mono") %in% seu$celltype)) {
    de <- FindMarkers(seu, ident.1 = "FCGR3A Mono", ident.2 = "CD14 Mono", min.pct = 0.25, verbose = FALSE)
    utils::write.csv(cbind(gene = rownames(de), de), file.path(a$out, "de_FCGR3A_vs_CD14_monocytes.csv"), row.names = FALSE)
  }
  t_types <- intersect(c("Naive CD4 T", "Memory CD4 T", "CD8 T"), unique(seu$celltype))
  if (length(t_types) && sum(seu$celltype %in% t_types) >= 50) {
    tc <- subcluster(seu, t_types, seed = seed)
    ggsave(file.path(a$out, "umap_tcell_subclusters.png"), DimPlot(tc, label = TRUE), width = 6, height = 6)
    utils::write.csv(FindAllMarkers(tc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE),
                     file.path(a$out, "tcell_subcluster_markers.csv"), row.names = FALSE)
    log$tcell_subclusters <- length(unique(tc$seurat_clusters))
  }
  log$clusters <- length(unique(seu$seurat_clusters))
  log$celltype_counts <- as.list(table(seu$celltype))
  log$ambiguous_clusters <- ann$cluster[ann$margin < 0.1]
  saveRDS(seu, file.path(a$out, "pbmc_annotated.rds"))
  jsonlite::write_json(log, file.path(a$out, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 4)
  invisible(seu)
}

if (sys.nframe() == 0) main()
