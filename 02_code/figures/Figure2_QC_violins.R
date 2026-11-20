# Run from the repository root; all paths below are repository-relative.
# NOTE: Seurat .rds objects under 04_seurat_objects/ are not distributed due to size.
# Regenerate them with the pipeline in 02_code/ (from the GEO raw data in the README),
# or place your own copy at the path shown before running.
# ============================================================================
# Figure 2 | Batch-level QC violin plots (before vs after filtering, by GSE dataset)
# Zebrafish heart-injury scRNA-seq meta-analysis
#
# INPUT(S):  04_seurat_objects/all_samples_seurat_list_raw.rds
# OUTPUT:    Figure_2.png
# NOTE: Extracted verbatim from the master script scRNAseq_Meta_Analysis_Github_2026.r (QC section). Six stray note lines were commented out so the file runs standalone.
# ============================================================================
# ---- QC and Filtering ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)  # for comma formatting

# ---- Load cleaned Seurat list ----
# ---- Load cleaned Seurat list with metadata ----
seurat_list <- readRDS("04_seurat_objects/all_samples_seurat_list_raw.rds")

# ---- QC thresholds ----
min_features <- 200
max_features <- 3500
max_percent_mt <- 30   # stricter for clean data
max_percent_ribo <- 50 # optional, for plotting

# ---- Font size settings ----
axis_text_size <- 14
axis_title_size <- 16
subfigure_title_size <- 16
main_title_size <- 18

# ---- Create folder for QC figures if not exists ----
dir.create("05_figures/qc", recursive = TRUE, showWarnings = FALSE)

# ---- Compute QC metrics and filter samples ----
qc_before <- list()
qc_after  <- list()
filtered_list <- list()

for (gsm in names(seurat_list)) {
  obj <- seurat_list[[gsm]]
  
  # % mitochondrial genes
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-|^MT-")
  
  # % ribosomal genes (rps/rpl)
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^rps|^rpl")
  
  # Save QC before filtering
  qc_before[[gsm]] <- data.frame(
    GSM = gsm,
    nFeature_RNA = obj$nFeature_RNA,
    nCount_RNA = obj$nCount_RNA,
    percent.mt = obj$percent.mt,
    percent.ribo = obj$percent.ribo,
    Batch = obj$GSE
  )
  
  # Filter cells
  obj_filt <- subset(obj, 
                     subset = nFeature_RNA > min_features &
                              nFeature_RNA < max_features &
                              percent.mt < max_percent_mt)
  
  # Save QC after filtering
  qc_after[[gsm]] <- data.frame(
    GSM = gsm,
    nFeature_RNA = obj_filt$nFeature_RNA,
    nCount_RNA = obj_filt$nCount_RNA,
    percent.mt = obj_filt$percent.mt,
    percent.ribo = obj_filt$percent.ribo,
    Batch = obj_filt$GSE
  )
  
  filtered_list[[gsm]] <- obj_filt
}

# Combine QC data for plotting
qc_before_df <- do.call(rbind, qc_before)
qc_after_df  <- do.call(rbind, qc_after)

# ---- Save QC metrics to CSV ----
library(dplyr)

# ---- Summarize QC metrics ----
qc_before_summary <- qc_before_df %>%
  group_by(GSM) %>%
  summarise(
    nCells = n(),
    Mean_nFeature = mean(nFeature_RNA),
    Mean_nCount = mean(nCount_RNA),
    Mean_percentMT = mean(percent.mt),
    Mean_percentRibo = mean(percent.ribo)
  )

qc_after_summary <- qc_after_df %>%
  group_by(GSM) %>%
  summarise(
    nCells = n(),
    Mean_nFeature = mean(nFeature_RNA),
    Mean_nCount = mean(nCount_RNA),
    Mean_percentMT = mean(percent.mt),
    Mean_percentRibo = mean(percent.ribo)
  )

# ---- Save to CSV ----
write.csv(qc_before_summary, "03_processed_data/qc/qc_metrics_before_summary.csv", row.names = FALSE)
write.csv(qc_after_summary, "03_processed_data/qc/qc_metrics_after_summary.csv", row.names = FALSE)

# ---- Preview head ----
head(qc_before_summary)
head(qc_after_summary)

# This will give you exactly the table you showed:
# nCells -> number of cells in that sample
# Mean_nFeature -> mean detected genes per cell
# Mean_nCount -> mean UMI counts per cell
# Mean_percentMT -> mean mitochondrial gene percentage per cell
# Mean_percentRibo -> mean ribosomal gene percentage per cell

# ---- Metrics to plot ----
metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")

# ---- Violin plot function ----
make_violin <- function(df, metric){
  p <- ggplot(df, aes(x = Batch, y = .data[[metric]], fill = Batch)) +
       geom_violin(trim = TRUE, scale = "width", color = "black", draw_quantiles = c(0.25, 0.5, 0.75)) +
       geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
       labs(y = metric, x = "") +
       theme_minimal() +
       theme(
         axis.text = element_text(size = axis_text_size, face = "bold"),
         axis.title = element_text(size = axis_title_size, face = "bold"),
         plot.title = element_text(size = subfigure_title_size, face = "bold"),
         legend.position = "none"
       )
  if (metric == "nCount_RNA") p <- p + scale_y_continuous(labels = comma)
  return(p)
}

# ---- Generate plots ----
vln_before <- lapply(metrics, function(m) make_violin(qc_before_df, m))
vln_after  <- lapply(metrics, function(m) make_violin(qc_after_df, m))

# ---- Combine plots: top = before, bottom = after ----
full_plot <- wrap_plots(vln_before, ncol = length(metrics)) /
             wrap_plots(vln_after, ncol = length(metrics)) +
             plot_annotation(title = "Batch-level QC: Top = Before QC | Bottom = After QC") &
             theme(plot.title = element_text(size = main_title_size, face = "bold"))

# ---- Save figure ----
ggsave(filename = "05_figures/main/Figure_2.png",
       plot = full_plot, width = 24, height = 14, dpi = 600)

# ---- Save filtered Seurat objects ----
saveRDS(filtered_list, file = "04_seurat_objects/Seurat_List_QC_filter.rds")
message("✅ QC completed, violin plots saved, and filtered Seurat objects saved!")
