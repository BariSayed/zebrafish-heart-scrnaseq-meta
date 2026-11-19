library(Seurat)
o <- readRDS("~/zebrafish/results/rds/Seurat_Integrated_SCT_with_UMAP_Clusters.rds")

o$old_cluster_SCT <- o$seurat_clusters

DefaultAssay(o) <- "integrated"
o <- RunPCA(o, npcs = 50, verbose = FALSE)
o <- FindNeighbors(o, reduction = "pca", dims = 1:50)
o <- FindClusters(o, resolution = 0.5)
o$new_cluster_integrated <- o$seurat_clusters

cat("old clusters (SCT):", length(levels(o$old_cluster_SCT)), "\n")
cat("new clusters (integrated):", length(levels(o$new_cluster_integrated)), "\n\n")

tab <- table(Old = o$old_cluster_SCT, New = o$new_cluster_integrated)
write.csv(as.data.frame.matrix(tab), "~/zebrafish/results/old_vs_new_cluster_crosstab.csv")
cat("saved: old_vs_new_cluster_crosstab.csv\n")

o <- RunUMAP(o, reduction = "pca", dims = 1:50)
saveRDS(o, "~/zebrafish/results/rds/Seurat_Integrated_reclustered.rds")
cat("saved: Seurat_Integrated_reclustered.rds\n")
cat("ALL DONE\n")
