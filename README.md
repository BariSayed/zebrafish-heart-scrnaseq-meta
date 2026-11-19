# Single-cell meta-analysis of zebrafish heart regeneration

Analysis code, processed data, figures, and tables for a meta-analysis of
single-cell transcriptomic (scRNA-seq) responses to cardiac injury in zebrafish
(*Danio rerio*). Three public GEO datasets — **GSE227190**, **GSE159032**,
**GSE262689** — are integrated into a single atlas (22 samples, 0–14 days
post-injury; **77,375 cells, 26 clusters** after SCT/rPCA integration).

## What is in this repository
```
00_metadata/              sample sheet (22 samples: GSM, GSE, dpi, condition)
02_code/
   01_integration/        SCT integration (rPCA) + reclustering
   02_annotation/         cluster markers, metadata export, composition
   03_differential_expression/  pseudobulk DESeq2 + mt/ribo flagging + checks
   figures/               one script per manuscript figure (Figure2..8, README inside)
   tables/                supplementary-table builder
03_processed_data/
   qc/                    per-sample QC summaries, before & after filtering
   clustering/            full cell metadata, 26-cluster annotation, marker tables
   differential_expression/  pseudobulk counts + DESeq2 FINAL (+ intermediate/)
   enrichment/            GO-BP & KEGG per-cluster (UP/DOWN)
   pseudotime/            Monocle3 pseudotime value tables (per cell type)
05_figures/               main/ (Figure_1..8) + pseudotime/ panels
06_tables/                manuscript tables (Table 2 QC metrics, Table 3 annotation)
07_supplementary/         submitted supplementary workbook
```

**Not in this repository (by design):**
- **Raw data** — public in NCBI GEO: GSE227190, GSE159032, GSE262689.
- **Integrated Seurat objects** (`.rds`, ~4.6 GB) — exceed GitHub's file limit;
  archived on Zenodo (DOI in the manuscript's Data Availability) or available on request.
- **Manuscript file** — not deposited here.

## Reproduction order
1. `02_code/01_integration/Seurat_SCT_Integration.R` — QC-filtered per-sample
   list → SCTransform (3000 variable features, regress `percent.mt`, glmGamPoi) →
   rPCA integration (2500 anchor features) → `Seurat_Integrated_SCT.rds`.
2. `02_code/01_integration/reclustering.R` — PCA/UMAP + clustering at
   resolution 0.5 → `Seurat_Integrated_reclustered.rds` (77,375 cells, 26 clusters).
3. `02_code/02_annotation/markers.R` — `PrepSCTFindMarkers` + `FindAllMarkers`
   (SCT assay, min.pct 0.25, logfc 0.25) → cluster marker tables.
4. `02_code/03_differential_expression/` — pseudobulk aggregation → DESeq2
   (`~ Dataset + Condition`) → `flag_mt_ribo.R` → `pseudobulk_DESeq2_results.csv`.
5. `02_code/figures/Figure6A_enrichment_analysis.R` — GO-BP & KEGG per cluster.
6. `02_code/figures/Figure7A_pseudotime_compute_4celltypes.R` and
   `Figure7B_pseudotime_compute_epicardial.R` — Monocle3 trajectories.
7. `02_code/figures/Figure2..8` — all manuscript figures (see
   `02_code/figures/README.md` for the figure↔script map and inputs).

## Software environments
- **scrna** — R 4.3.3; Seurat 5.4.0, SeuratObject 5.3.0, sctransform 0.4.3,
  glmGamPoi (integration, clustering, markers, Figures 2 & 4).
- **metaseq** — R 4.4.3; DESeq2 1.46.0, clusterProfiler 4.14.0,
  org.Dr.eg.db 3.20.0, AnnotationDbi 1.68.0 (pseudobulk DE, enrichment).
- **system R** 4.5.1 + Monocle3 (pseudotime, Figures 7A/7B).
- **ComplexHeatmap** 2.26.1 + circlize (Figure 5 composition heatmap).
- **Python 3** + matplotlib/pandas for the composite figure scripts
  (Figures 3, 6B, 7C, 8).

Full package versions are recorded in `environment.txt`.

## Figures
Figure scripts live in `02_code/figures/`, named by manuscript figure number
(a script per panel where a figure combines several). Each was traced and
verified against the final figure image; details in `02_code/figures/README.md`.
Figure 5 (cell-type composition heatmap) uses a 0–10 % colour cap so the many
mid/low-abundance clusters remain visible (a few dominant clusters reach ~42 %
and are drawn at the darkest shade).

## Data & licensing
- Code: **MIT** (`LICENSE`).
- Processed data, figures, and tables: **CC-BY-4.0** (`LICENSE_data_CC-BY-4.0`).
- Raw data: NCBI GEO (original dataset licences apply).

## Citation
See `CITATION.cff`. Please cite the associated manuscript and this repository if
you use the code or processed data.
