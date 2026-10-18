suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(root, "R", "pipeline.R"))
source(file.path(root, "run_pipeline.R"))

write_10x <- function(dir, n_per = 150, seed = 1) {
  set.seed(seed)
  types <- list(T = c("CD3D", "IL7R", "CCR7", "LEF1", "CD8A", "CD8B"), B = c("MS4A1", "CD79A", "CD79B"),
                Mono = c("CD14", "LYZ", "S100A8", "CST3"))
  genes <- c(unique(unlist(types)), paste0("MT-CO", 1:3), sprintf("GENE%04d", 1:1200))
  cells <- list()
  for (t in names(types)) {
    m <- matrix(rpois(length(genes) * n_per, 0.6), length(genes), n_per, dimnames = list(genes, NULL))
    m[types[[t]], ] <- rpois(length(types[[t]]) * n_per, 12)
    m[sprintf("GENE%04d", sample(1200, 60)), ] <- rpois(60 * n_per, 4)     # type-specific programme
    cells[[t]] <- m
  }
  x <- do.call(cbind, cells)
  colnames(x) <- sprintf("CELL%04d-1", seq_len(ncol(x)))
  dir.create(dir, recursive = TRUE)
  Matrix::writeMM(methods::as(x, "CsparseMatrix"), file.path(dir, "matrix.mtx"))
  writeLines(colnames(x), file.path(dir, "barcodes.tsv"))
  utils::write.table(data.frame(genes, genes, "Gene Expression"), file.path(dir, "features.tsv"),
                     sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  for (f in c("matrix.mtx", "barcodes.tsv", "features.tsv")) system2("gzip", file.path(dir, f))
  rep(names(types), each = n_per)
}

test_that("expected doublet rate follows 10x guidance", {
  expect_equal(expected_doublet_rate(2700), 0.0216, tolerance = 1e-6)
  expect_lt(expected_doublet_rate(2700), 0.075)
})

test_that("full 10x run: QC, doublets with discovered column, clustering, annotation", {
  d <- withr::local_tempdir()
  truth <- write_10x(file.path(d, "tenx"))
  out <- file.path(d, "out")
  seu <- suppressWarnings(suppressMessages(main(c("--tenx", file.path(d, "tenx"), "--out", out, "--dims", "10", "--resolution", "0.3"))))
  s <- jsonlite::read_json(file.path(out, "summary.json"))
  expect_equal(s$cells_raw, 450)
  expect_true(s$doublets$nExp <= ceiling(0.008 * 450 / 1000 * 450))   # ~1.6 doublets expected, not 34
  expect_equal(s$cells_analysed, s$cells_after_qc - s$doublets$n_called)
  ann <- read.csv(file.path(out, "cluster_annotation.csv"))
  expect_true(all(c("B", "CD14 Mono") %in% ann$celltype))
  expect_true(file.exists(file.path(out, "pbmc_annotated.rds")))
})

test_that("subclustering recomputes PCA inside the subset", {
  d <- withr::local_tempdir()
  write_10x(file.path(d, "tenx"), n_per = 120)
  seu <- CreateSeuratObject(Read10X(file.path(d, "tenx")))
  seu <- cluster_cells(preprocess(seu, n_pcs = 15), 1:10, 0.3)
  seu$celltype <- ifelse(seq_len(ncol(seu)) <= 240, "T", "other")
  sub <- suppressWarnings(subcluster(seu, "T", dims = 1:5))
  expect_equal(ncol(sub), 240)
  expect_false(isTRUE(all.equal(Embeddings(sub, "pca")[, 1], Embeddings(seu, "pca")[colnames(sub), 1])))
})
