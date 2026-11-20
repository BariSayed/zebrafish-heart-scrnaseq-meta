# Run from the repository root; all paths below are repository-relative.
# NOTE: Seurat .rds objects under 04_seurat_objects/ are not distributed due to size.
# Regenerate them with the pipeline in 02_code/ (from the GEO raw data in the README),
# or place your own copy at the path shown before running.
library(Seurat)
library(dplyr)

o <- readRDS("04_seurat_objects/Seurat_Integrated_reclustered.rds")
DefaultAssay(o) <- "SCT"
o <- PrepSCTFindMarkers(o)

markers <- FindAllMarkers(o, assay = "SCT", slot = "data",
                           only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)
write.csv(top10, "03_processed_data/clustering/new_cluster_top10_markers.csv", row.names = FALSE)
write.csv(markers, "03_processed_data/clustering/new_cluster_all_markers.csv", row.names = FALSE)
cat("ALL DONE\n")
