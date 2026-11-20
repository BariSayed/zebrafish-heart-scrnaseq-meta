# Run from the repository root; all paths below are repository-relative.
# NOTE: Seurat .rds objects under 04_seurat_objects/ are not distributed due to size.
# Regenerate them with the pipeline in 02_code/ (from the GEO raw data in the README),
# or place your own copy at the path shown before running.
#!/usr/bin/env Rscript

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(Seurat)
  library(glmGamPoi) # speeds up SCTransform
})

message("🚀 Starting SCT normalization and integration ...")

# ---- Load filtered Seurat list ----
filtered_list <- readRDS("04_seurat_objects/Seurat_List_QC_filter.rds")
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
saveRDS(seurat_integrated, "04_seurat_objects/Seurat_Integrated_SCT.rds")
message("💾 Saved integrated object to 04_seurat_objects/Seurat_Integrated_SCT.rds")

message("🎉 All done!")

