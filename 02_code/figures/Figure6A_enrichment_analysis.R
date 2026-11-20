# Run from the repository root; all paths below are repository-relative.
# ============================================================================
# Figure 6A | GO-BP & KEGG enrichment per cluster (clusterProfiler) + per-cluster dotplots
# Zebrafish heart-injury scRNA-seq meta-analysis
#
# INPUT(S):  pseudobulk_DESeq2_results_primary.csv
# OUTPUT:    GO_BP/KEGG_enrichment_{UP,DOWN}_per_cluster.csv + per-cluster dotplot PNGs
# NOTE: Produces the enrichment tables that Figure6B curates into the 4-panel figure.
# ============================================================================
# =====================================================================================
# STEP 9: GO (Biological Process) and KEGG pathway enrichment, PER CELL TYPE / CLUSTER
# Input: pseudobulk_DESeq2_results_primary.csv (mt-/rps-rpl genes already excluded, Step
#        8 pseudobulk DESeq2 output). Up/down gene sets built separately per cluster
#        (padj < 0.05), using a SHARED detected-gene universe (all genes tested by DESeq2
#        in at least one cluster) converted to Entrez IDs via org.Dr.eg.db. No hand-
#        filtering beyond the padj<0.05 / mt-/ribosomal exclusion already documented.
# =====================================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Dr.eg.db)
  library(enrichplot)
})

# Run from the repository root; paths below are repository-relative.
de <- read.csv("03_processed_data/differential_expression/intermediate/pseudobulk_DESeq2_results_primary.csv", stringsAsFactors = FALSE)
cat("Loaded PRIMARY DE table:", nrow(de), "rows,", length(unique(de$gene)), "unique genes,",
    length(unique(de$cluster_id)), "clusters\n")

# ---- 1. Build the shared detected-gene universe (all genes DESeq2 actually tested) ----
universe_symbols <- unique(de$gene)
cat("Universe (all tested genes, pre-Entrez-mapping):", length(universe_symbols), "\n")

map_universe <- bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dr.eg.db, drop = TRUE)
cat("Universe genes mapped to Entrez ID:", nrow(map_universe), "of", length(universe_symbols),
    sprintf("(%.1f%%)\n", 100 * nrow(map_universe) / length(universe_symbols)))
universe_entrez <- unique(map_universe$ENTREZID)

# ---- 2. Build per-cluster up/down Entrez gene lists (padj < 0.05) ----
clusters <- sort(unique(de$cluster_id))
up_list <- list()
down_list <- list()
cluster_labels <- c()

for (cl in clusters) {
  sub <- de[de$cluster_id == cl & !is.na(de$padj), ]
  label <- unique(sub$cell_type)[1]
  cluster_labels[as.character(cl)] <- label

  up_sym <- sub$gene[sub$padj < 0.05 & sub$log2FoldChange > 0]
  down_sym <- sub$gene[sub$padj < 0.05 & sub$log2FoldChange < 0]

  up_entrez <- unique(map_universe$ENTREZID[map_universe$SYMBOL %in% up_sym])
  down_entrez <- unique(map_universe$ENTREZID[map_universe$SYMBOL %in% down_sym])

  cat("Cluster", cl, "(", label, "): up =", length(up_sym), "sym /", length(up_entrez), "entrez;",
      "down =", length(down_sym), "sym /", length(down_entrez), "entrez\n")

  if (length(up_entrez) >= 5) up_list[[paste0(cl, "_", label)]] <- up_entrez
  if (length(down_entrez) >= 5) down_list[[paste0(cl, "_", label)]] <- down_entrez
}

cat("\nClusters with >=5 mapped UP genes (eligible for enrichment):", length(up_list), "\n")
cat("Clusters with >=5 mapped DOWN genes (eligible for enrichment):", length(down_list), "\n")

# ---- 3. GO (Biological Process) enrichment, compareCluster, shared universe ----
run_compare <- function(gene_list, label) {
  if (length(gene_list) == 0) {
    cat("SKIP", label, "- no eligible clusters\n")
    return(NULL)
  }
  tryCatch({
    compareCluster(geneCluster = gene_list, fun = "enrichGO",
                    OrgDb = org.Dr.eg.db, keyType = "ENTREZID", ont = "BP",
                    universe = universe_entrez, pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) {
    cat("ERROR in GO compareCluster for", label, ":", conditionMessage(e), "\n")
    NULL
  })
}

go_up <- run_compare(up_list, "GO_UP")
go_down <- run_compare(down_list, "GO_DOWN")

if (!is.null(go_up)) {
  write.csv(as.data.frame(go_up), "03_processed_data/enrichment/GO_BP_enrichment_UP_per_cluster.csv", row.names = FALSE)
  cat("Saved GO_BP_enrichment_UP_per_cluster.csv,", nrow(as.data.frame(go_up)), "rows\n")
}
if (!is.null(go_down)) {
  write.csv(as.data.frame(go_down), "03_processed_data/enrichment/GO_BP_enrichment_DOWN_per_cluster.csv", row.names = FALSE)
  cat("Saved GO_BP_enrichment_DOWN_per_cluster.csv,", nrow(as.data.frame(go_down)), "rows\n")
}

# ---- 4. KEGG pathway enrichment, compareCluster, organism = dre (zebrafish) ----
run_compare_kegg <- function(gene_list, label) {
  if (length(gene_list) == 0) {
    cat("SKIP", label, "- no eligible clusters\n")
    return(NULL)
  }
  tryCatch({
    compareCluster(geneCluster = gene_list, fun = "enrichKEGG",
                    organism = "dre", keyType = "ncbi-geneid",
                    universe = universe_entrez, pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) {
    cat("ERROR in KEGG compareCluster for", label, ":", conditionMessage(e), "\n")
    NULL
  })
}

kegg_up <- run_compare_kegg(up_list, "KEGG_UP")
kegg_down <- run_compare_kegg(down_list, "KEGG_DOWN")

if (!is.null(kegg_up)) {
  write.csv(as.data.frame(kegg_up), "03_processed_data/enrichment/KEGG_enrichment_UP_per_cluster.csv", row.names = FALSE)
  cat("Saved KEGG_enrichment_UP_per_cluster.csv,", nrow(as.data.frame(kegg_up)), "rows\n")
}
if (!is.null(kegg_down)) {
  write.csv(as.data.frame(kegg_down), "03_processed_data/enrichment/KEGG_enrichment_DOWN_per_cluster.csv", row.names = FALSE)
  cat("Saved KEGG_enrichment_DOWN_per_cluster.csv,", nrow(as.data.frame(kegg_down)), "rows\n")
}

# ---- 5. Dotplot figures (only for non-null results) ----
save_dotplot <- function(obj, fname, title) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) {
    cat("SKIP dotplot", fname, "- no results\n")
    return(invisible(NULL))
  }
  tryCatch({
    p <- dotplot(obj, showCategory = 5) + ggplot2::ggtitle(title)
    ggplot2::ggsave(fname, p, width = 12, height = 10, dpi = 300, limitsize = FALSE)
    cat("Saved", fname, "\n")
  }, error = function(e) cat("ERROR plotting", fname, ":", conditionMessage(e), "\n"))
}

save_dotplot(go_up, "05_figures/Figure_GO_BP_UP_per_cluster.png", "GO (BP) enrichment - Upregulated genes per cell type")
save_dotplot(go_down, "05_figures/Figure_GO_BP_DOWN_per_cluster.png", "GO (BP) enrichment - Downregulated genes per cell type")
save_dotplot(kegg_up, "05_figures/Figure_KEGG_UP_per_cluster.png", "KEGG enrichment - Upregulated genes per cell type")
save_dotplot(kegg_down, "05_figures/Figure_KEGG_DOWN_per_cluster.png", "KEGG enrichment - Downregulated genes per cell type")

cat("\nSTEP9_DONE\n")
