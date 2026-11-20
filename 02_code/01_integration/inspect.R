# Run from the repository root; all paths below are repository-relative.
# NOTE: Seurat .rds objects under 04_seurat_objects/ are not distributed due to size.
# Regenerate them with the pipeline in 02_code/ (from the GEO raw data in the README),
# or place your own copy at the path shown before running.
library(Seurat)
o <- readRDS("04_seurat_objects/Seurat_Integrated_SCT_with_UMAP_Clusters.rds")
sink("~/inspect_out.txt")
cat("assays:", names(o@assays), "\n")
cat("DefaultAssay:", DefaultAssay(o), "\n")
cat("reductions:", names(o@reductions), "\n")
cat("PCA ran on assay:", o@reductions$pca@assay.used, "\n")
cat("UMAP ran on assay:", o@reductions$umap@assay.used, "\n")
cat("cells:", ncol(o), " genes:", nrow(o), "\n")
cat("n clusters:", length(levels(o$seurat_clusters)), "\n")
cat("meta cols:", paste(colnames(o@meta.data), collapse=", "), "\n")
if ("integrated" %in% names(o@assays)) {
  cat("integrated assay dims:", nrow(o@assays$integrated), ncol(o@assays$integrated), "\n")
}
sink()
cat("DONE\n")
