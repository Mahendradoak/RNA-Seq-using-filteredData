# Single-cell RNA-seq of 3k PBMCs (Seurat v5)

> Despite the repository name, this is a **single-cell** RNA-seq analysis, not bulk RNA-seq. Consider renaming the repository to `pbmc3k-scrnaseq`.

A full analysis of the 10x Genomics 3k peripheral blood mononuclear cell (PBMC) dataset from a healthy donor ([figshare 28414916](https://figshare.com/articles/dataset/3k_PBMCs_from_a_healthy_donor/28414916)).

## What it does

1. **Quality control.** Keeps cells with 200–2,500 detected genes and less than 5% mitochondrial reads.
2. **Preprocessing.** Log-normalisation, 2,000 variable genes, scaling and PCA.
3. **Doublet removal (DoubletFinder):**
   - pK is optimised by parameter sweep;
   - the expected doublet rate follows 10x guidance (~0.8% per 1,000 cells, about 2% here);
   - the homotypic-doublet adjustment is applied;
   - the same PCs are used as for clustering;
   - the classification column is found by name rather than hard-coded.
4. **Clustering.** Louvain clustering and UMAP on 15 PCs.
5. **Cell-type annotation from canonical marker module scores** at the cluster level: naive/memory CD4 T, CD8 T, NK, B, CD14 and FCGR3A monocytes, DC and platelets. Each cluster's score margin over the runner-up is reported, and clusters with a margin below 0.1 are flagged as ambiguous. Reference-based annotation with SingleR is optional (use an immune reference such as Monaco).
6. **Markers and differential expression:**
   - all-cluster markers;
   - a biologically defined contrast, FCGR3A⁺ versus CD14⁺ monocytes.
7. **T-cell subclustering.** Variable genes, PCA and UMAP are recomputed **inside** the T-cell subset.

## Run

```bash
Rscript install_packages.R
# From the raw 10x matrix (full pipeline):
Rscript run_pipeline.R --tenx 28414916/filtered_feature_bc_matrix --out results
# From the saved counts in pbmc_processed.rds (cells already QC'd and doublet-filtered):
Rscript run_pipeline.R --rds pbmc_processed.rds --out results_from_saved_counts
Rscript tests/testthat.R      # 9 expectations on a synthetic 10x dataset
```

## Results (`results_from_saved_counts/`, real data)

The pipeline was run (Seurat 5.0.1) from the counts stored in `pbmc_processed.rds`: 2,483 cells, which had already been through QC and the original doublet removal. Cluster numbers and margins can shift slightly with other Seurat versions.

| Cluster | Cell type | Cells | Marker margin |
|---|---|---|---|
| 0 | Naive CD4 T | 546 | 0.100 |
| 1 | Memory CD4 T | 504 | 0.472 |
| 2 | CD14 Mono | 456 | 1.578 |
| 3 | B | 335 | 1.878 |
| 4 | CD8 T | 318 | **0.008 (ambiguous: CD8 effector vs NK)** |
| 5 | FCGR3A Mono | 166 | 1.232 |
| 6 | NK | 141 | 2.378 |
| 7 | Platelet | 17 | 3.975 |

Other outputs:
- **Plots:** UMAP by cluster and by cell type, a marker dot plot, and T-cell subclusters (4 found).
- **Tables:** marker tables and the FCGR3A vs CD14 monocyte differential expression.

Dendritic cells (about 1% of PBMCs) did not form their own cluster at this resolution.

The **doublet step should be re-run from the raw matrix** with the corrected settings. The saved object had already lost about 200 cells to the 7.5% assumption, where about 2% (~58 cells) was expected.

## What was wrong before (`archive/`)

| Problem | Effect | Fix |
|---|---|---|
| Expected doublet rate 7.5% (`nExp = 0.075 × cells`) | About 3× too many cells removed as doublets (201 vs ~58) | 10x rate (~0.8%/1,000 cells) with homotypic adjustment |
| DoubletFinder column hard-coded as `DF.classifications_0.25_0.01_201` while pK was computed | Script breaks, or uses the wrong column, whenever the optimal pK differs from 0.01 | Column found by pattern |
| Doublets detected on PCs 1:10, clustering on PCs 1:15 | Inconsistent embeddings | Same PCs throughout |
| "T cells = cluster 0", "monocytes = cluster 2" hard-coded | Cluster numbers change between runs and versions | Subsets chosen by annotation |
| Subclusters reused the parent PCA and UMAP | Subcluster structure driven by genes that separate the major lineages | Variable genes and PCA recomputed within the subset |
| SingleR with HumanPrimaryCellAtlas `label.main` per cell | Coarse, partly non-PBMC labels ("CMP", "Myelocyte", "Pro-B_cell_CD34+") | Cluster-level marker scores; optional immune reference |
| Differential expression "cluster 0 vs cluster 4" | No biological question | FCGR3A vs CD14 monocytes |
| README referred to `pbmc_analysis.R`; the script was `program.R`; Seurat v4 `slot=` API | Instructions did not work | `run_pipeline.R`, Seurat v5 `layer=` API |

**Statistical note:** with a single donor, cell-level Wilcoxon tests describe this sample only. Population-level claims need several donors and pseudobulk testing.
