# Run from the repository root; all paths below are repository-relative.
# ============================================================================
# Figure 5 | Cell-type composition heatmap (ComplexHeatmap; Condition/GSE/DPI bars, 26 clusters)
# Zebrafish heart-injury scRNA-seq meta-analysis
#
# INPUT:   03_processed_data/clustering/full_metadata_reclustered.csv
#          (self-contained: rebuilds the cluster x GSM % matrix from metadata,
#           so no intermediate prep CSVs are required)
# OUTPUT:  Figure_5.png
#
# Blue gradient; top Condition/GSE/DPI annotation bars; cell-type labels on the
# right; legends stacked on the right. Colour scale capped at 10% so the many
# mid/low-abundance clusters stay visible (a few dominant clusters reach ~42%
# and are drawn at the darkest blue).
# Run with an R that has ComplexHeatmap + circlize, e.g.:
#   Rscript 02_code/figures/Figure5_composition_heatmap.R
# ============================================================================
#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ComplexHeatmap); library(circlize); library(grid)})

N   <- "."
md  <- read.csv(file.path(N,"03_processed_data/clustering/full_metadata_reclustered.csv"),
                stringsAsFactors=FALSE, check.names=FALSE)
cl  <- as.character(md[["integrated_snn_res.0.5"]])
gsm <- as.character(md$GSM)

gsm_order <- c("GSM4817932","GSM4817939","GSM4817941",
               "GSM7092936","GSM7092937","GSM7092944","GSM7092945","GSM7092952","GSM7092953",
               paste0("GSM", 8173831:8173843))

# cluster -> label exactly as in the published figure (data come from the cluster
# column; labels are the figure's short names appended with the cluster number)
labs <- c("Endocardial cell","Immature erythroid progenitor","Cardiomyocytes","Epicardial cells",
 "Erythroid-thrombocyte progenitor","Cardiomyocytes ventricular","Macrophages",
 "Low-quality / mitochondrial-high","Endocardial cell","Smooth muscle cell / pericyte",
 "T cells / lymphocytes","Thrombocytes","Endocardial cell (frzb+)","Erythroid cells",
 "Vascular endothelial fenestrated","Erythroid cells","Proliferating / cycling cells",
 "Erythroid-thrombocyte progenitor","Neutrophils","Fibroblasts","Erythroid-thrombocyte progenitor",
 "Perivascular cell / pericyte","Epicardial cell","Erythroid cells",
 "apln+ injury-responsive endothelial","Lymphatic endothelial cells")
rowlab <- paste0(labs, " (", 0:25, ")")

tab <- table(factor(cl, levels=as.character(0:25)), factor(gsm, levels=gsm_order))
pct <- sweep(tab, 2, colSums(tab), "/") * 100
mat <- as.matrix(pct); rownames(mat) <- rowlab; colnames(mat) <- gsm_order

u <- function(col) sapply(gsm_order, function(g) md[[col]][md$GSM==g][1])
anno <- data.frame(GSM=gsm_order, Condition=u("Condition"), GSE=u("GSE"),
                   dpi=gsub("dpi","",as.character(u("dpi"))), stringsAsFactors=FALSE)

cond_col <- c(Healthy="#2CA02C", Injured="#D62728")
gse_col  <- c(GSE227190="#7FD4C6", GSE159032="#FBF6B4", GSE262689="#C7BBDD")
dpi_col  <- c(`0`="#1F77B4", `1`="#FF7F0E", `3`="#2CA02C", `7`="#D62728", `14`="#9467BD")

top_anno <- HeatmapAnnotation(
  Condition=anno$Condition, GSE=anno$GSE, DPI=anno$dpi,
  col=list(Condition=cond_col, GSE=gse_col, DPI=dpi_col),
  annotation_name_side="right", annotation_name_gp=gpar(fontsize=11, fontface="bold"),
  simple_anno_size=unit(4,"mm"), gap=unit(0.6,"mm"),
  annotation_legend_param=list(
    Condition=list(title="Condition", title_gp=gpar(fontface="bold",fontsize=11), labels_gp=gpar(fontsize=10)),
    GSE=list(title="GSE", title_gp=gpar(fontface="bold",fontsize=11), labels_gp=gpar(fontsize=10)),
    DPI=list(title="DPI", at=c("0","1","3","7","14"), labels=c("0dpi","1dpi","3dpi","7dpi","14dpi"),
             title_gp=gpar(fontface="bold",fontsize=11), labels_gp=gpar(fontsize=10))))

cap <- 10
col_fun <- colorRamp2(c(0, cap/2, cap), c("#FFFFFF","#9ECAE1","#2171B5"))

ht <- Heatmap(mat, name="% of sample", col=col_fun, top_annotation=top_anno,
  cluster_rows=FALSE, cluster_columns=FALSE,
  row_names_side="right", row_names_gp=gpar(fontsize=10.5),
  column_names_side="bottom", column_names_gp=gpar(fontsize=9), column_names_rot=90,
  rect_gp=gpar(col="white", lwd=0.5),
  heatmap_legend_param=list(title="% of sample\n(capped at 10)", title_gp=gpar(fontface="bold",fontsize=11),
                            labels_gp=gpar(fontsize=10), legend_height=unit(3,"cm")),
  width=unit(ncol(mat)*0.55,"cm"), height=unit(nrow(mat)*0.55,"cm"))

png("05_figures/main/Figure_5.png", width=3200, height=2000, res=230)
draw(ht, heatmap_legend_side="right", annotation_legend_side="right",
     merge_legend=TRUE, padding=unit(c(2,2,2,8),"mm"))
dev.off()
cat("saved Figure_5.png; cap =", cap, "; matrix max% =", round(max(mat),1), "\n")
