suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("/home/sbau/zebrafish/results/rds/Seurat_Integrated_reclustered.rds")
cat("ncells:", ncol(obj), "\n")
cat("meta cols:", paste(colnames(obj@meta.data), collapse=", "), "\n")
md <- obj@meta.data
write.csv(md, "full_metadata_reclustered.csv", row.names=TRUE)
cat("saved metadata csv, nrow=", nrow(md), "\n")
