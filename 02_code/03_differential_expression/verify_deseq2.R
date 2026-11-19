# Verification checks on pseudobulk_DESeq2_results_all_clusters.csv
df <- read.csv("~/zebrafish/results/pseudobulk_DESeq2_results_all_clusters.csv", stringsAsFactors = FALSE)
cat("Total rows:", nrow(df), "\n")
cat("Columns:", paste(colnames(df), collapse=", "), "\n\n")

cat("NA counts per column:\n")
print(sapply(df, function(x) sum(is.na(x))))

cat("\nlog2FoldChange summary:\n")
print(summary(df$log2FoldChange))

cat("\npadj summary:\n")
print(summary(df$padj))

cat("\nAny Inf/-Inf in log2FoldChange:", sum(is.infinite(df$log2FoldChange)), "\n")

cat("\nSignificant genes (padj<0.05) per cluster:\n")
print(table(df$cell_type[df$padj < 0.05 & !is.na(df$padj)]))

cat("\nTop 5 upregulated genes overall (padj<0.05, highest log2FC):\n")
sig <- df[!is.na(df$padj) & df$padj < 0.05, ]
top_up <- sig[order(-sig$log2FoldChange), ][1:5, c("gene","cluster_id","cell_type","log2FoldChange","padj","baseMean")]
print(top_up)

cat("\nTop 5 downregulated genes overall (padj<0.05, lowest log2FC):\n")
top_down <- sig[order(sig$log2FoldChange), ][1:5, c("gene","cluster_id","cell_type","log2FoldChange","padj","baseMean")]
print(top_down)

cat("\nCardiomyocyte cluster (2) top 10 by padj:\n")
cm <- df[df$cluster_id == 2 & !is.na(df$padj), ]
cm <- cm[order(cm$padj), ][1:10, c("gene","log2FoldChange","padj","baseMean")]
print(cm)

cat("\nVERIFY_DONE\n")
