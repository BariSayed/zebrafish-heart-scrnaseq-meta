# Build final flagged DESeq2 results table: mt-/rps/rpl genes flagged (not deleted) for transparency.
de <- read.csv("~/zebrafish/results/pseudobulk_DESeq2_results_all_clusters.csv", stringsAsFactors = FALSE)

de$gene_class <- ifelse(grepl("^mt-", de$gene, ignore.case = TRUE), "mitochondrial",
                  ifelse(grepl("^rps|^rpl", de$gene, ignore.case = TRUE), "ribosomal_protein",
                         "other"))
de$exclude_from_primary <- de$gene_class %in% c("mitochondrial", "ribosomal_protein")

cat("Gene class breakdown (rows, all clusters combined):\n")
print(table(de$gene_class))

write.csv(de, "~/zebrafish/results/pseudobulk_DESeq2_results_FINAL.csv", row.names = FALSE)
cat("\nSaved pseudobulk_DESeq2_results_FINAL.csv with gene_class and exclude_from_primary columns,", nrow(de), "rows\n")

# Also save a "primary" (clean) subset for GO/KEGG universe + headline tables
de_primary <- de[!de$exclude_from_primary, ]
write.csv(de_primary, "~/zebrafish/results/pseudobulk_DESeq2_results_PRIMARY.csv", row.names = FALSE)
cat("Saved pseudobulk_DESeq2_results_PRIMARY.csv (mt-/rps/rpl excluded),", nrow(de_primary), "rows\n")

cat("\nFLAG_DONE\n")
