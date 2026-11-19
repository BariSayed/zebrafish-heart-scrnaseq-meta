library(Seurat)
library(dplyr)

o <- readRDS("~/zebrafish/results/rds/Seurat_Integrated_reclustered.rds")
DefaultAssay(o) <- "SCT"
o <- PrepSCTFindMarkers(o)

markers <- FindAllMarkers(o, assay = "SCT", slot = "data",
                           only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)
write.csv(top10, "~/zebrafish/results/new_cluster_top10_markers.csv", row.names = FALSE)
write.csv(markers, "~/zebrafish/results/new_cluster_all_markers.csv", row.names = FALSE)
cat("ALL DONE\n")
