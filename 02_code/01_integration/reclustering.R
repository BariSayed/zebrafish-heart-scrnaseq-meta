# Run from the repository root; all paths below are repository-relative.
# NOTE: Seurat .rds objects under 04_seurat_objects/ are not distributed due to size.
# Regenerate them with the pipeline in 02_code/ (from the GEO raw data in the README),
# or place your own copy at the path shown before running.
library(Seurat)
o <- readRDS("04_seurat_objects/Seurat_Integrated_SCT_with_UMAP_Clusters.rds")

o$old_cluster_SCT <- o$seurat_clusters

DefaultAssay(o) <- "integrated"
o <- RunPCA(o, npcs = 50, verbose = FALSE)
o <- FindNeighbors(o, reduction = "pca", dims = 1:50)
o <- FindClusters(o, resolution = 0.5)
o$new_cluster_integrated <- o$seurat_clusters

cat("old clusters (SCT):", length(levels(o$old_cluster_SCT)), "\n")
cat("new clusters (integrated):", length(levels(o$new_cluster_integrated)), "\n\n")

tab <- table(Old = o$old_cluster_SCT, New = o$new_cluster_integrated)
write.csv(as.data.frame.matrix(tab), "03_processed_data/clustering/old_vs_new_cluster_crosstab.csv")
cat("saved: old_vs_new_cluster_crosstab.csv\n")

o <- RunUMAP(o, reduction = "pca", dims = 1:50)
saveRDS(o, "04_seurat_objects/Seurat_Integrated_reclustered.rds")
cat("saved: Seurat_Integrated_reclustered.rds\n")
cat("ALL DONE\n")
