# ============================================================================
# Figure 7A | Monocle3 pseudotime computation for 4 cell types (CM, Endocardial, Fibroblast, Macrophage)
# Zebrafish heart-injury scRNA-seq meta-analysis
#
# INPUT(S):  ~/zebrafish/results/rds/Seurat_Integrated_reclustered.rds
# OUTPUT:    pseudotime/<CellType>_Pseudotime.png + <CellType>_pseudotime_values.csv
# NOTE: Compute step; its *_pseudotime_values.csv feed Figure7C and Figure8. Run with system R 4.5.1 (monocle3).
# ============================================================================
#!/usr/bin/env Rscript
# ============================================================================
# Pseudotime trajectory inference (Monocle3) for zebrafish heart injury
# scRNA-seq meta-analysis.
#
# Per Methods Section 7: for each annotated cell type containing >=50 cells
# (Cardiomyocytes, Endocardial cells, Fibroblasts, Macrophages), the
# corresponding subset is extracted from the integrated Seurat object and
# converted into a Monocle3 CellDataSet via as.cell_data_set(). Trajectories
# are inferred with learn_graph() after cluster_cells(). Pseudotime ordering
# roots the trajectory at 0 dpi (uninjured/healthy) cells. Output is
# visualized with plot_cells(), mapping the transition across 0, 1, 3, 7,
# and 14 dpi.
#
# Run with: /usr/bin/Rscript run_monocle3_pseudotime.R
# (uses system R 4.5.1 personal library, where monocle3/SeuratWrappers/
#  Seurat are installed -- the scrna conda env, R 4.3.3, lacks monocle3
#  and SeuratWrappers and must not be used for this step.)
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(ggplot2)
})

set.seed(1234)

rds_path     <- "/home/sbau/zebrafish/results/rds/Seurat_Integrated_reclustered.rds"
cluster_col  <- "integrated_snn_res.0.5"
out_dir      <- "/home/sbau/zebrafish/results/pseudotime"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading integrated Seurat object...\n")
obj <- readRDS(rds_path)
stopifnot(cluster_col %in% colnames(obj@meta.data))
cat("Loaded object:", ncol(obj), "cells\n")

# Join RNA layers ONCE on the full object before any subsetting. Joining
# per-subset instead (as in an earlier version of this script) hit a
# Seurat v5 bug for small clusters (e.g. 876-cell Fibroblasts), where a
# per-sample layer with zero cells in the subset caused JoinLayers() to
# fail with "incorrect number of dimensions" inside StitchMatrix/LayerData.
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- NormalizeData(obj, verbose = FALSE)
cat("RNA layers joined and normalized on full object.\n")

# Cell-type -> cluster mapping, from new_cluster_annotation_FINAL.csv
celltype_clusters <- list(
  "Cardiomyocytes"    = c(2, 5),
  "Endocardial_cell"  = c(0, 8, 12),   # includes cluster 12, "Endocardial cell (frzb+)"
  "Fibroblasts"       = c(19),
  "Macrophages"       = c(6)
)

run_pseudotime_for_celltype <- function(name, clusters) {
  cat("\n=== Processing:", name, "(clusters", paste(clusters, collapse=","), ") ===\n")

  cells_keep <- colnames(obj)[obj@meta.data[[cluster_col]] %in% clusters]
  n_cells <- length(cells_keep)
  cat("  n cells:", n_cells, "\n")
  if (n_cells < 50) {
    cat("  SKIPPED: fewer than 50 cells\n")
    return(invisible(NULL))
  }

  sub <- subset(obj, cells = cells_keep)
  DefaultAssay(sub) <- "RNA"

  if (!"dpi" %in% colnames(sub@meta.data)) {
    stop("Column 'dpi' not found in metadata -- cannot root trajectory at 0 dpi.")
  }
  n_root <- sum(sub$dpi == "0dpi", na.rm = TRUE)
  cat("  n 0-dpi (root candidate) cells:", n_root, "\n")
  if (n_root == 0) {
    cat("  SKIPPED: no 0 dpi cells present in this subset -- cannot root trajectory\n")
    return(invisible(NULL))
  }

  cds <- as.cell_data_set(sub)
  # reuse Seurat's existing UMAP embedding (already computed during integration)
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = TRUE)

  root_cells <- colnames(sub)[sub$dpi == "0dpi"]
  cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)

  p <- plot_cells(cds,
                   color_cells_by = "pseudotime",
                   label_groups_by_cluster = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   graph_label_size = 3) +
       ggtitle(paste0(gsub("_", " ", name), " (n=", n_cells, ")"))

  out_png <- file.path(out_dir, paste0(name, "_Pseudotime.png"))
  ggsave(out_png, p, width = 7, height = 6, dpi = 300)
  cat("  saved:", out_png, "\n")

  # also save the pseudotime values + dpi/condition for downstream stats
  pt_df <- data.frame(
    cell = colnames(cds),
    pseudotime = monocle3::pseudotime(cds),
    dpi = sub$dpi,
    Condition = sub$Condition,
    Dataset = sub$Dataset
  )
  write.csv(pt_df, file.path(out_dir, paste0(name, "_pseudotime_values.csv")), row.names = FALSE)

  invisible(cds)
}

for (nm in names(celltype_clusters)) {
  tryCatch(
    run_pseudotime_for_celltype(nm, celltype_clusters[[nm]]),
    error = function(e) cat("  ERROR for", nm, ":", conditionMessage(e), "\n")
  )
}

cat("\n=== DONE: Monocle3 pseudotime for all 4 cell types ===\n")
