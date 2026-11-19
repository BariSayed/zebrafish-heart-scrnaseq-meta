
d <- read.csv("/home/sbau/zebrafish/results/full_metadata_reclustered.csv")
cat("total cells:", nrow(d), "\n")
cat("Condition table:\n"); print(table(d$Condition))

# group clusters
endo <- c(0,8,12); cm <- c(2,5); epi <- c(3,22)
mac <- 6; tcell <- 10; neut <- 18

pct_by_cond <- function(clusters, label){
  for(cond in c("Healthy","Injured")){
    sub <- d[d$Condition==cond,]
    n_total <- nrow(sub)
    n_clust <- sum(sub$new_cluster_integrated %in% clusters)
    cat(sprintf("%s | %s: %d/%d = %.2f%%\n", label, cond, n_clust, n_total, 100*n_clust/n_total))
  }
}
pct_by_cond(endo, "Endocardial(0,8,12)")
pct_by_cond(cm, "CM(2,5)")
pct_by_cond(epi, "Epicardial(3,22)")
pct_by_cond(mac, "Macrophage(6)")
pct_by_cond(tcell, "Tcell(10)")
pct_by_cond(neut, "Neutrophil(18)")
