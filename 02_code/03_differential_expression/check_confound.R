# Confound check: percent.mt / percent.ribo by Condition & Dataset (from lightweight CSV, no RDS reload)
meta <- read.csv("~/zebrafish/results/full_metadata_reclustered.csv", stringsAsFactors = FALSE)
cat("Metadata rows:", nrow(meta), "\n")
cat("Columns:", paste(colnames(meta), collapse=", "), "\n\n")

if (all(c("percent.mt","percent.ribo","Condition","Dataset") %in% colnames(meta))) {
  cat("Mean percent.mt by Condition x Dataset:\n")
  print(aggregate(percent.mt ~ Dataset + Condition, data = meta, FUN = mean))
  cat("\nMean percent.ribo by Condition x Dataset:\n")
  print(aggregate(percent.ribo ~ Dataset + Condition, data = meta, FUN = mean))
} else {
  cat("STOP: expected columns not found, got:", paste(colnames(meta), collapse=", "), "\n")
}

# Check DESeq2 results for mt-/rps/rpl gene over-representation among significant hits
de <- read.csv("~/zebrafish/results/pseudobulk_DESeq2_results_all_clusters.csv", stringsAsFactors = FALSE)
is_mt  <- grepl("^mt-", de$gene, ignore.case = TRUE)
is_rib <- grepl("^rps|^rpl", de$gene, ignore.case = TRUE)
cat("\nTotal mt- genes in results (rows, all clusters combined):", sum(is_mt), "\n")
cat("Total rps/rpl genes in results (rows, all clusters combined):", sum(is_rib), "\n")

sig <- !is.na(de$padj) & de$padj < 0.05
cat("\nOverall significance rate (all genes):", round(mean(sig), 4), "\n")
cat("Significance rate among mt- genes:", round(mean(sig[is_mt]), 4), "(n=", sum(is_mt), ")\n")
cat("Significance rate among rps/rpl genes:", round(mean(sig[is_rib]), 4), "(n=", sum(is_rib), ")\n")

cat("\nAll significant mt- gene hits (gene, cluster_id, cell_type, log2FC, padj):\n")
print(de[is_mt & sig, c("gene","cluster_id","cell_type","log2FoldChange","padj")])

cat("\nAll significant rps/rpl gene hits (gene, cluster_id, cell_type, log2FC, padj):\n")
print(de[is_rib & sig, c("gene","cluster_id","cell_type","log2FoldChange","padj")])

cat("\nCONFOUND_CHECK_DONE\n")
