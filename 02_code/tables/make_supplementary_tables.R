# Run from the repository root; all paths below are repository-relative.
suppressPackageStartupMessages({
  library(org.Dr.eg.db)
  library(AnnotationDbi)
})

annot_path  <- "00_metadata/cluster_annotation_simple.txt"
marker_path <- "03_processed_data/clustering/new_cluster_top10_markers.csv"
out_path    <- "06_tables/Supplementary_Table_ClusterMarkers_detailed.csv"

annot <- read.delim(annot_path, sep = "|", stringsAsFactors = FALSE, check.names = FALSE)
markers <- read.csv(marker_path, stringsAsFactors = FALSE, check.names = FALSE)

cat("Annotation rows:", nrow(annot), "| Marker rows:", nrow(markers), "\n")

genes <- unique(markers$gene)
desc_map <- rep(NA_character_, length(genes))
names(desc_map) <- genes

lookup <- tryCatch(
  AnnotationDbi::mapIds(org.Dr.eg.db, keys = genes, keytype = "SYMBOL",
                        column = "GENENAME", multiVals = "first"),
  error = function(e) { cat("SYMBOL lookup note:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(lookup)) desc_map[names(lookup)] <- lookup

missing <- names(desc_map)[is.na(desc_map)]
if (length(missing) > 0) {
  alias_try <- tryCatch(
    AnnotationDbi::mapIds(org.Dr.eg.db, keys = missing, keytype = "ALIAS",
                          column = "GENENAME", multiVals = "first"),
    error = function(e) NULL
  )
  if (!is.null(alias_try)) {
    for (g in names(alias_try)) if (!is.na(alias_try[g])) desc_map[g] <- alias_try[g]
  }
}
desc_map[is.na(desc_map)] <- "not annotated in org.Dr.eg.db"

out_list <- list()
for (cl in sort(annot$cluster)) {
  a <- annot[annot$cluster == cl, ][1, ]
  g <- markers[markers$cluster == cl, ]
  g <- g[order(-g$avg_log2FC), ]
  g <- head(g, 10)
  if (nrow(g) == 0) next
  out_list[[as.character(cl)]] <- data.frame(
    cluster         = cl,
    final_cell_type = a$final_cell_type,
    n_cells         = a$n_cells,
    pct_total       = a$pct_total,
    gene            = g$gene,
    avg_log2FC      = round(g$avg_log2FC, 3),
    pct.1           = g$pct.1,
    pct.2           = g$pct.2,
    p_val_adj       = g$p_val_adj,
    gene_function   = desc_map[g$gene],
    reference       = a$reference,
    stringsAsFactors = FALSE, check.names = FALSE
  )
}
out <- do.call(rbind, out_list)
rownames(out) <- NULL

write.csv(out, out_path, row.names = FALSE)

n_annotated   <- sum(out$gene_function != "not annotated in org.Dr.eg.db")
n_unannotated <- sum(out$gene_function == "not annotated in org.Dr.eg.db")
cat("\n=== DONE ===\n")
cat("Total gene rows:", nrow(out), "across", length(out_list), "clusters\n")
cat("Genes WITH database function:", n_annotated, "\n")
cat("Genes NOT annotated (flagged):", n_unannotated, "\n")
cat("Output:", out_path, "\n")
