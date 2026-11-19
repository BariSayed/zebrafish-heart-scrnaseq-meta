# ============================================================================
# Figure 7B | Monocle3 pseudotime computation for Epicardial cells (clusters 3+22)
# Zebrafish heart-injury scRNA-seq meta-analysis
#
# INPUT(S):  ~/zebrafish/results/rds/Seurat_Integrated_reclustered.rds
# OUTPUT:    pseudotime/Epicardial_cell_Pseudotime.png + Epicardial_cell_pseudotime_values.csv
# NOTE: Same pipeline as 6A, separate because epicardial clusters are handled on their own. Run with system R 4.5.1.
# ============================================================================
#!/usr/bin/env Rscript
# Monocle3 pseudotime for Epicardial cells (clusters 3 + 22), identical pipeline
# to run_monocle3_pseudotime.R (global JoinLayers fix; root at 0 dpi).
suppressPackageStartupMessages({
  library(Seurat); library(SeuratWrappers); library(monocle3); library(ggplot2)
})
set.seed(1234)
rds_path    <- "/home/sbau/zebrafish/results/rds/Seurat_Integrated_reclustered.rds"
cluster_col <- "integrated_snn_res.0.5"
out_dir     <- "/home/sbau/zebrafish/results/pseudotime"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(rds_path)
stopifnot(cluster_col %in% colnames(obj@meta.data))
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- NormalizeData(obj, verbose = FALSE)

name <- "Epicardial_cell"; clusters <- c(3, 22)
cells_keep <- colnames(obj)[obj@meta.data[[cluster_col]] %in% clusters]
cat("n cells:", length(cells_keep), "\n")
sub <- subset(obj, cells = cells_keep)
DefaultAssay(sub) <- "RNA"
n_root <- sum(sub$dpi == "0dpi", na.rm = TRUE)
cat("n 0-dpi root cells:", n_root, "\n")
stopifnot(n_root > 0)

cds <- as.cell_data_set(sub)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)
root_cells <- colnames(sub)[sub$dpi == "0dpi"]
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)

p <- plot_cells(cds, color_cells_by = "pseudotime",
                label_groups_by_cluster = FALSE, label_leaves = FALSE,
                label_branch_points = FALSE, graph_label_size = 3) +
     ggtitle(paste0("Epicardial cell (n=", length(cells_keep), ")"))
ggsave(file.path(out_dir, "Epicardial_cell_Pseudotime.png"), p, width = 7, height = 6, dpi = 300)

pt_df <- data.frame(cell = colnames(cds), pseudotime = monocle3::pseudotime(cds),
                    dpi = sub$dpi, Condition = sub$Condition, Dataset = sub$Dataset)
write.csv(pt_df, file.path(out_dir, "Epicardial_cell_pseudotime_values.csv"), row.names = FALSE)

# quick stats to stdout/log
pt <- pt_df$pseudotime; pt <- pt[is.finite(pt)]
cat("pseudotime summary:\n"); print(summary(pt))
agg <- aggregate(pseudotime ~ dpi, data = pt_df[is.finite(pt_df$pseudotime),], FUN = median)
cat("median pseudotime by dpi:\n"); print(agg)
cat("median by condition:\n"); print(aggregate(pseudotime ~ Condition, data=pt_df[is.finite(pt_df$pseudotime),], FUN=median))
cat("DONE\n")
