# Tested with R 4.3.3, Seurat 5.0.1, DoubletFinder 2.0.6, SingleR 2.4.1.
install.packages(c("Seurat", "ggplot2", "jsonlite", "testthat", "withr", "remotes", "BiocManager"),
                 repos = "https://cloud.r-project.org")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = "never")
BiocManager::install(c("SingleR", "celldex"), update = FALSE, ask = FALSE)   # optional reference annotation
# To build a reference once:  saveRDS(celldex::MonacoImmuneData(), "monaco.rds")
