# ============================================================================
# Figure 4 | Canonical marker FeaturePlot grid, 16 panels (8 cell types x 2 markers)
# Zebrafish heart-injury scRNA-seq meta-analysis
#
# INPUT(S):  ~/zebrafish/results/rds/Seurat_Integrated_reclustered.rds
# OUTPUT:    Figure_4.png
# NOTE: Build-time file name says 'Figure6' but content = manuscript Figure 5.
# ============================================================================
library(Seurat)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

obj <- readRDS("~/zebrafish/results/rds/Seurat_Integrated_reclustered.rds")
cat("Loaded object:", ncol(obj), "cells\n")

cluster_col <- "seurat_clusters"
n_levels <- length(unique(obj@meta.data[[cluster_col]]))
cat("Using cluster column:", cluster_col, "- levels:", n_levels, "\n")
if (n_levels != 26) stop("STOP: expected 26 levels, got ", n_levels, ". Check cluster column before proceeding.")

DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)

gene_panels <- list(
  "Cardiomyocytes"      = c("myl7", "cmlc1"),
  "Endocardial cell"    = c("f8", "notch1b"),
  "Epicardial cells"    = c("tcf21", "tbx18"),
  "Macrophages"         = c("c1qa", "marco"),
  "Fibroblasts"         = c("rspo1", "ptx3a"),
  "T cells"             = c("ptprc", "ccr9a"),
  "Neutrophils"         = c("mpx", "lyz"),
  "Perivascular cell"   = c("pdgfrb", "rgs4")
)

all_genes <- unlist(gene_panels)
missing_genes <- setdiff(all_genes, rownames(obj))
if (length(missing_genes) > 0) {
  cat("WARNING: genes not found in object, will be skipped:", paste(missing_genes, collapse=", "), "\n")
}
all_genes <- intersect(all_genes, rownames(obj))
cat("Genes to plot (", length(all_genes), "of 16 ):", paste(all_genes, collapse=", "), "\n")

panel_letters <- letters[1:length(all_genes)]
plots <- list()
idx <- 1
for (ct in names(gene_panels)) {
  for (g in gene_panels[[ct]]) {
    if (!(g %in% rownames(obj))) next
    p <- FeaturePlot(obj, features = g, reduction = "umap", pt.size = 0.15, order = TRUE) +
      ggtitle(paste0(panel_letters[idx], "  ", ct, ": ", g)) +
      theme(plot.title = element_text(size = 13, face = "plain", hjust = 0),
            axis.title = element_text(size = 11),
            axis.text = element_text(size = 10),
            legend.text = element_text(size = 9),
            legend.title = element_text(size = 9),
            legend.key.width = unit(0.15, "cm"))
    plots[[idx]] <- p
    idx <- idx + 1
  }
}

combined <- wrap_plots(plots, ncol = 4)
ggsave("Figure_4.png", combined, width = 16, height = 16, dpi = 300, limitsize = FALSE)
cat("DONE: saved Figure_4.png with", length(plots), "panels\n")
