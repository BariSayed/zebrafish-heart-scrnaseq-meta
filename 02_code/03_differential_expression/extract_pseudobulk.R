suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

obj <- readRDS("~/zebrafish/results/rds/Seurat_Integrated_reclustered.rds")
cat("Loaded object:", ncol(obj), "cells x", nrow(obj), "genes\n")
cat("Available metadata columns:\n"); print(colnames(obj@meta.data))

required_cols <- c("GSM", "Dataset", "Condition", "dpi")
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("STOP: expected metadata columns missing: ", paste(missing_cols, collapse=", "))
}

# HARDCODED (verified 2026-10-17): 'seurat_clusters' and 'integrated_snn_res.0.5' are confirmed
# IDENTICAL (the correct 26-cluster integrated-assay solution). 'SCT_snn_res.0.5' (34 levels) is
# the OLD/WRONG SCT-based clustering and must NOT be used.
cluster_col <- "seurat_clusters"
n_levels <- length(unique(obj@meta.data[[cluster_col]]))
cat("Using cluster column:", cluster_col, "- levels:", n_levels, "\n")
if (n_levels != 26) {
  stop("STOP: expected 26 cluster levels in '", cluster_col, "', found ", n_levels,
       ". Do not proceed -- re-verify which column holds the integrated-assay solution.")
}

annot <- read.csv("new_cluster_annotation_FINAL.csv")
cluster_to_celltype <- setNames(annot$final_cell_type, as.character(annot$new_cluster))

obj$cluster_id <- as.character(obj@meta.data[[cluster_col]])
obj$cell_type  <- unname(cluster_to_celltype[obj$cluster_id])

obj <- subset(obj, subset = cluster_id != "7")
cat("After excluding cluster 7:", ncol(obj), "cells\n")

DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
cat("RNA counts matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "cells\n")

meta <- obj@meta.data
meta$pseudo_group <- paste(meta$GSM, meta$cluster_id, sep = "__")

MIN_CELLS_PER_PSEUDOBULK <- 10
group_sizes <- table(meta$pseudo_group)
keep_groups <- names(group_sizes[group_sizes >= MIN_CELLS_PER_PSEUDOBULK])
cat("Pseudobulk groups retained (>=", MIN_CELLS_PER_PSEUDOBULK, "cells):", length(keep_groups), "of", length(group_sizes), "\n")

pb_list <- list(); pb_meta_list <- list()
for (g in keep_groups) {
  cells_in_g <- rownames(meta)[meta$pseudo_group == g]
  pb_list[[g]] <- Matrix::rowSums(counts_mat[, cells_in_g, drop = FALSE])
  parts <- strsplit(g, "__")[[1]]; gsm_g <- parts[1]; cl_g <- parts[2]
  m_row <- meta[meta$GSM == gsm_g, ][1, ]
  pb_meta_list[[g]] <- data.frame(
    pseudo_group = g, GSM = gsm_g, cluster_id = cl_g,
    cell_type = cluster_to_celltype[cl_g],
    Dataset = m_row$Dataset, Condition = m_row$Condition, dpi = m_row$dpi,
    n_cells = length(cells_in_g)
  )
}
pb_counts <- do.call(cbind, pb_list)
pb_meta <- do.call(rbind, pb_meta_list)
rownames(pb_meta) <- pb_meta$pseudo_group

write.csv(as.data.frame(as.matrix(pb_counts)), "pseudobulk_counts_matrix.csv", row.names = TRUE)
write.csv(pb_meta, "pseudobulk_sample_metadata.csv", row.names = FALSE)
cat("STAGE1_DONE: saved pseudobulk_counts_matrix.csv (", nrow(pb_counts), "genes x", ncol(pb_counts),
    "pseudo-samples) and pseudobulk_sample_metadata.csv (", nrow(pb_meta), "rows)\n")
