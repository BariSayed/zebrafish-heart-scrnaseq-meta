#!/usr/bin/env Rscript

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(Seurat)
  library(glmGamPoi) # speeds up SCTransform
})

message("🚀 Starting SCT normalization and integration ...")

# ---- Load filtered Seurat list ----
filtered_list <- readRDS("~/zebrafish/results/rds/Seurat_List_QC_filter.rds")
message("✅ Loaded filtered Seurat list: ", length(filtered_list), " samples")

# ---- SCT normalization per sample ----
sct_list <- lapply(filtered_list, function(obj) {
  SCTransform(
    obj,
    assay = "RNA",
    new.assay.name = "SCT",
    variable.features.n = 3000,
    vars.to.regress = "percent.mt",
    method = "glmGamPoi",
    verbose = FALSE
  )
})

message("✅ SCT normalization completed for all samples")

# ---- Integration ----
integration_features <- SelectIntegrationFeatures(
  object.list = sct_list,
  nfeatures = 2500
)

sct_list <- PrepSCTIntegration(
  object.list = sct_list,
  anchor.features = integration_features
)

anchors <- FindIntegrationAnchors(
  object.list = sct_list,
  normalization.method = "SCT",
  anchor.features = integration_features,
  reduction = "rpca" # RPCA = faster and more memory efficient
)

seurat_integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

message("✅ Integration completed")

# ---- Save output ----
saveRDS(seurat_integrated, "~/zebrafish/results/rds/Seurat_Integrated_SCT.rds")
message("💾 Saved integrated object to ~/zebrafish/results/rds/Seurat_Integrated_SCT.rds")

message("🎉 All done!")

