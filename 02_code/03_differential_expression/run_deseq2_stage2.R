# =====================================================================================
# STAGE 2: Pseudobulk differential expression (DESeq2), Injured vs Healthy, per cell type
# Reads the Stage-1 outputs produced by extract_pseudobulk.R:
#   pseudobulk_counts_matrix.csv   (genes x pseudo-samples, raw summed RNA-assay counts)
#   pseudobulk_sample_metadata.csv (one row per pseudo-sample: pseudo_group, GSM,
#                                    cluster_id, cell_type, Dataset, Condition, dpi, n_cells)
# Design: ~ Dataset + Condition when >=2 datasets both contain Healthy & Injured samples
# for that cluster; otherwise falls back to ~ Condition (documented per-cluster in output).
# =====================================================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

setwd("~/zebrafish/results")

# ---- 1. Load Stage-1 outputs ----
counts_mat <- read.csv("pseudobulk_counts_matrix.csv", row.names = 1, check.names = FALSE)
counts_mat <- as.matrix(counts_mat)
pb_meta <- read.csv("pseudobulk_sample_metadata.csv", stringsAsFactors = FALSE)

cat("Loaded pseudobulk counts:", nrow(counts_mat), "genes x", ncol(counts_mat), "pseudo-samples\n")
cat("Loaded pseudobulk metadata:", nrow(pb_meta), "rows\n")

# ---- 2. Sanity check: counts columns and metadata pseudo_group must match exactly ----
if (!identical(sort(colnames(counts_mat)), sort(pb_meta$pseudo_group))) {
  stop("STOP: pseudobulk_counts_matrix.csv columns do not match pseudobulk_sample_metadata.csv ",
       "pseudo_group values. Re-check Stage 1 outputs before proceeding.")
}
# reorder metadata rows to match counts column order exactly
pb_meta <- pb_meta[match(colnames(counts_mat), pb_meta$pseudo_group), ]
rownames(pb_meta) <- pb_meta$pseudo_group
cat("Verified: counts matrix columns and metadata rows are aligned.\n")

# ---- 3. Run DESeq2 per cluster (cell type), design chosen per cluster ----
all_results <- list()
clusters_present <- sort(unique(pb_meta$cluster_id))
cat("Clusters to test:", paste(clusters_present, collapse = ", "), "\n\n")

for (cl in clusters_present) {
  cl_meta <- pb_meta[pb_meta$cluster_id == cl, ]
  cl_counts <- counts_mat[, cl_meta$pseudo_group, drop = FALSE]
  cell_type_label <- unique(cl_meta$cell_type)

  if (length(unique(cl_meta$Condition)) < 2) {
    cat("SKIP cluster", cl, "(", cell_type_label, "): only one Condition level present\n")
    next
  }
  if (ncol(cl_counts) < 4) {
    cat("SKIP cluster", cl, "(", cell_type_label, "): fewer than 4 pseudobulk samples (n=",
        ncol(cl_counts), "), underpowered\n")
    next
  }

  cond_table <- table(cl_meta$Dataset, cl_meta$Condition)
  n_datasets <- length(unique(cl_meta$Dataset))
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
  }, error = function(e) {
    cat("SKIP cluster", cl, "- DESeq2 design matrix error:", conditionMessage(e), "\n")
    return(NULL)
  })
  if (is.null(dds)) next

  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  if (nrow(dds) < 10) {
    cat("SKIP cluster", cl, "- too few genes pass expression filter (", nrow(dds), ")\n")
    next
  }

  dds <- tryCatch(DESeq(dds), error = function(e) {
    cat("SKIP cluster", cl, "- DESeq() failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  if (is.null(dds)) next

  res <- results(dds, contrast = c("Condition", "Injured", "Healthy"), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$cluster_id <- cl
  res_df$cell_type <- cell_type_label
  res_df$design_used <- design_used
  res_df$n_pseudosamples <- ncol(cl_counts)
  all_results[[as.character(cl)]] <- res_df

  cat("DONE cluster", cl, "(", cell_type_label, "): n=", ncol(cl_counts),
      "design=", design_used, "sig genes (padj<0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
}

final_results <- do.call(rbind, all_results)
write.csv(final_results, "pseudobulk_DESeq2_results_all_clusters.csv", row.names = FALSE)
cat("\n=== STAGE2_DONE: saved pseudobulk_DESeq2_results_all_clusters.csv,",
    nrow(final_results), "gene x cluster rows across", length(all_results), "clusters tested ===\n")
