# =====================================================================================
# Pseudobulk differential expression: Injured vs Healthy, per cell type, ~Dataset+Condition
# Built on VERIFIED metadata columns from the original script's metadata tribble:
#   GSM, GSE, dpi, Condition (Healthy/Injured), Replicate, Dataset (== GSE), Platform
# Uses the RE-CLUSTERED, integrated-assay object (26 clusters) -- NOT the old SCT-based
# 34-cluster object. RAW COUNTS are taken from the "RNA" assay (DESeq2 requires raw/integer
# counts; the "SCT" and "integrated" assays are NOT appropriate inputs for DESeq2).
# =====================================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(dplyr)
  library(Matrix)
})

# ---- 1. Load the reclustered object ----
obj <- readRDS("~/zebrafish/results/rds/Seurat_Integrated_reclustered.rds")
cat("Loaded object:", ncol(obj), "cells x", nrow(obj), "genes\n")
cat("Available metadata columns:\n")
print(colnames(obj@meta.data))

# ---- 2. Sanity-check required columns exist BEFORE proceeding (no assuming) ----
required_cols <- c("GSM", "Dataset", "Condition", "dpi")
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("STOP: expected metadata columns missing: ", paste(missing_cols, collapse=", "),
       ". Re-check obj@meta.data column names above and edit this script before re-running.")
}

# Identify the cluster column (may be 'seurat_clusters' or 'integrated_snn_res.0.5' etc.)
cluster_col_candidates <- grep("seurat_clusters|snn_res", colnames(obj@meta.data), value = TRUE)
cat("Candidate cluster columns found:", paste(cluster_col_candidates, collapse=", "), "\n")
if (length(cluster_col_candidates) == 0) {
  stop("STOP: no cluster assignment column found in meta.data. Inspect colnames(obj@meta.data) and set cluster_col manually.")
}
cluster_col <- cluster_col_candidates[1]  # EDIT if the wrong one is picked -- check the printed list above first
cat("Using cluster column:", cluster_col, "\n")

# Attach the cell-type labels from new_cluster_annotation_FINAL.csv (must be uploaded alongside this script)
annot <- read.csv("new_cluster_annotation_FINAL.csv")
cluster_to_celltype <- setNames(annot$final_cell_type, as.character(annot$new_cluster))

obj$cluster_id <- as.character(obj@meta.data[[cluster_col]])
obj$cell_type  <- cluster_to_celltype[obj$cluster_id]

# Drop cluster 7 (confirmed non-cell-type / mitochondrial-high) from DE entirely
obj <- subset(obj, subset = cluster_id != "7")
cat("After excluding cluster 7: ", ncol(obj), "cells remain\n")

# ---- 3. Raw counts from the RNA assay (NOT SCT, NOT integrated) ----
DefaultAssay(obj) <- "RNA"
counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
cat("RNA counts matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "cells\n")

meta <- obj@meta.data
meta$pseudo_group <- paste(meta$GSM, meta$cluster_id, sep = "__")

# ---- 4. Aggregate raw counts into pseudobulk (sum per GSM x cluster) ----
groups <- unique(meta$pseudo_group)
cat("Total GSM x cluster pseudobulk groups (before size filter):", length(groups), "\n")

MIN_CELLS_PER_PSEUDOBULK <- 10  # standard minimum; document in Methods if changed
group_sizes <- table(meta$pseudo_group)
keep_groups <- names(group_sizes[group_sizes >= MIN_CELLS_PER_PSEUDOBULK])
cat("Pseudobulk groups retained (>=", MIN_CELLS_PER_PSEUDOBULK, "cells):", length(keep_groups), "of", length(groups), "\n")

pb_list <- list()
pb_meta_list <- list()
for (g in keep_groups) {
  cells_in_g <- rownames(meta)[meta$pseudo_group == g]
  pb_list[[g]] <- Matrix::rowSums(counts_mat[, cells_in_g, drop = FALSE])
  parts <- strsplit(g, "__")[[1]]
  gsm_g <- parts[1]; cl_g <- parts[2]
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

write.csv(pb_meta, "pseudobulk_sample_metadata.csv", row.names = FALSE)
cat("Saved pseudobulk_sample_metadata.csv (", nrow(pb_meta), "pseudo-samples)\n")

# ---- 5. Run DESeq2 PER CELL TYPE (cluster), design ~ Dataset + Condition ----
all_results <- list()
clusters_present <- sort(unique(pb_meta$cluster_id))

for (cl in clusters_present) {
  cl_meta <- pb_meta[pb_meta$cluster_id == cl, ]
  cl_counts <- pb_counts[, cl_meta$pseudo_group, drop = FALSE]

  n_datasets <- length(unique(cl_meta$Dataset))
  cond_table <- table(cl_meta$Dataset, cl_meta$Condition)
  if (length(unique(cl_meta$Condition)) < 2) {
    cat("SKIP cluster", cl, "(", unique(cl_meta$cell_type), "): only one Condition level present, cannot test Injured vs Healthy\n")
    next
  }
  if (ncol(cl_counts) < 4) {
    cat("SKIP cluster", cl, "(", unique(cl_meta$cell_type), "): fewer than 4 pseudobulk samples (n=", ncol(cl_counts), "), underpowered\n")
    next
  }

  dataset_has_both <- sapply(rownames(cond_table), function(d) sum(cond_table[d, ] > 0) >= 2)
  if (n_datasets > 1 && any(dataset_has_both)) {
    design_formula <- ~ Dataset + Condition
    design_used <- "~Dataset+Condition"
  } else {
    design_formula <- ~ Condition
    design_used <- "~Condition (Dataset dropped: confounded or single-dataset cluster)"
  }

  cl_meta$Condition <- factor(cl_meta$Condition, levels = c("Healthy", "Injured"))
  cl_meta$Dataset <- factor(cl_meta$Dataset)

  dds <- tryCatch({
    DESeqDataSetFromMatrix(countData = round(cl_counts), colData = cl_meta, design = design_formula)
  }, error = function(e) { cat("SKIP cluster", cl, "- DESeq2 design matrix error:", conditionMessage(e), "\n"); return(NULL) })
  if (is.null(dds)) next

  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  if (nrow(dds) < 10) {
    cat("SKIP cluster", cl, "- too few genes pass expression filter (", nrow(dds), ")\n")
    next
  }

  dds <- tryCatch(DESeq(dds), error = function(e) { cat("SKIP cluster", cl, "- DESeq() failed:", conditionMessage(e), "\n"); return(NULL) })
  if (is.null(dds)) next

  res <- results(dds, contrast = c("Condition", "Injured", "Healthy"), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$cluster_id <- cl
  res_df$cell_type <- unique(cl_meta$cell_type)
  res_df$design_used <- design_used
  res_df$n_pseudosamples <- ncol(cl_counts)
  all_results[[cl]] <- res_df
  cat("DONE cluster", cl, "(", unique(cl_meta$cell_type), "): n=", ncol(cl_counts),
      "design=", design_used, "sig genes (padj<0.05):", sum(res_df$padj < 0.05, na.rm=TRUE), "\n")
}

final_results <- do.call(rbind, all_results)
write.csv(final_results, "pseudobulk_DESeq2_results_all_clusters.csv", row.names = FALSE)
cat("\n=== FINAL: saved pseudobulk_DESeq2_results_all_clusters.csv,", nrow(final_results), "gene x cluster rows across", length(all_results), "clusters tested ===\n")
