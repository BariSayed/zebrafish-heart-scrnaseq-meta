# Figure-generation code

One script per manuscript figure (A/B/C suffix when a figure needs more than one
script). File names follow the **current manuscript** figure numbers.

| File | Manuscript figure | Language |
|---|---|---|
| `Figure2_QC_violins.R` | Fig 2 - batch QC violins (before/after × GSE) | R (Seurat) |
| `Figure3A_UMAP_overview.py` | Fig 3 (upper) - UMAP overview, 6 panels | Python |
| `Figure3B_UMAP_healthy_vs_injured.py` | Fig 3 (lower) - Healthy vs Injured split | Python |
| `Figure4_marker_featureplots.R` | Fig 4 - canonical-marker FeaturePlot grid (16 panels) | R (Seurat) |
| `Figure5_composition_heatmap.R` | Fig 5 - cell-type composition heatmap | R (ComplexHeatmap) |
| `Figure6A_enrichment_analysis.R` | Fig 6 - GO/KEGG enrichment per cluster (data) | R (clusterProfiler) |
| `Figure6B_enrichment_dotplot.py` | Fig 6 - curated GO/KEGG dot-plot, 4 panels | Python |
| `Figure7A_pseudotime_compute_4celltypes.R` | Fig 7 - Monocle3 pseudotime compute (CM/Endo/Fibro/Macro) | R (Monocle3) |
| `Figure7B_pseudotime_compute_epicardial.R` | Fig 7 - Monocle3 pseudotime compute (epicardial) | R (Monocle3) |
| `Figure7C_pseudotime_UMAP_plot.py` | Fig 7 - pseudotime-on-UMAP, 5 panels | Python |
| `Figure8_pseudotime_boxplot_by_dpi.py` | Fig 8 - pseudotime boxplots by dpi | Python |

Figure 1 is a BioRender schematic (no code).

## Output file names
Each plotting script writes an output named by its manuscript figure number:
`Figure_2.png`, `Figure_3A.png`, `Figure_3B.png`, `Figure_4.png`, `Figure_5.png`,
`Figure_6.png`, `Figure_7.png`, `Figure_8.png`. Figure7A/7B write per-cell-type
intermediates (`<CellType>_Pseudotime.png` + `<CellType>_pseudotime_values.csv`)
consumed by Figure7C and Figure8.

## Dependencies between scripts
```
Figure7A ┐
Figure7B ┼─► *_pseudotime_values.csv ─► Figure7C (Fig 7)
         └                            └─► Figure8  (Fig 8)
Figure6A ──► GO/KEGG_*_per_cluster.csv ─► Figure6B (Fig 6)
```

## Inputs each script expects
- Figs 4, 7A, 7B: `Seurat_Integrated_reclustered.rds` (77,375-cell object).
- Fig 2: `all_samples_seurat_list_raw.rds` (pre-QC per-sample list).
- Figs 3A/3B: `umap_full_export.csv(.gz)` + `new_cluster_annotation.csv`.
- Figs 7C/8: the five `*_pseudotime_values.csv` from 7A/7B.
- Fig 6B: `GO_BP/KEGG_enrichment_{UP,DOWN}_per_cluster.csv` from 6A.
- Fig 6A: `pseudobulk_DESeq2_results_primary.csv`.
- Fig 5: `comp_matrix_for_heatmap.csv` + `heatmap_col_annotation.csv`
  (built from `full_metadata_reclustered.csv`).

## Software
R scripts: scrna env (R 4.3.3, Seurat 5.4.0) for Figs 2/4; metaseq env
(R 4.4.3, clusterProfiler 4.14.0) for Fig 6A; system R 4.5.1 (Monocle3 1.4.26)
for Figs 7A/7B. Python scripts: matplotlib (+ figure-style helper for 3A & 6B).

## Note on numbering
Figures were renumbered after the manuscript figure order was rearranged. The
mapping from the previous numbering was: old 3->5 (composition), old 4->3 (UMAP),
old 5->4 (markers), old 6->7 (pseudotime UMAP), old 7->8 (pseudotime boxplot),
old 8->6 (enrichment); Figs 1 and 2 unchanged. File names, in-script output
names, and headers all reflect the current numbers.
