# 04_seurat_objects  (not included in this repository)

The Seurat objects (`.rds`) are not redistributed here because of their size.
They are produced by the pipeline in `02_code/` and consumed by the figure and
differential-expression scripts. Main objects in the chain:

- `all_samples_seurat_list_raw.rds`   - per-sample raw list (QC input)
- `Seurat_List_QC_filter.rds`         - after QC filtering (Figure2_QC_violins.R)
- `Seurat_Integrated_SCT.rds`         - after SCT/rPCA integration (01_integration/Seurat_SCT_Integration.R)
- `Seurat_Integrated_reclustered.rds` - final integrated + reclustered object (01_integration/reclustering.R)

Regenerate them by running the scripts in `02_code/01_integration/` (and
`02_code/02_annotation/`) starting from the raw data in `01_raw_data/`, then
place them here. Scripts that read these objects expect them at
`04_seurat_objects/<name>.rds`.
