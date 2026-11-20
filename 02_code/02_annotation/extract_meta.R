# Run from the repository root; all paths below are repository-relative.
# NOTE: Seurat .rds objects under 04_seurat_objects/ are not distributed due to size.
# Regenerate them with the pipeline in 02_code/ (from the GEO raw data in the README),
# or place your own copy at the path shown before running.
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("04_seurat_objects/Seurat_Integrated_reclustered.rds")
cat("ncells:", ncol(obj), "\n")
cat("meta cols:", paste(colnames(obj@meta.data), collapse=", "), "\n")
md <- obj@meta.data
write.csv(md, "03_processed_data/clustering/full_metadata_reclustered.csv", row.names=TRUE)
cat("saved metadata csv, nrow=", nrow(md), "\n")
