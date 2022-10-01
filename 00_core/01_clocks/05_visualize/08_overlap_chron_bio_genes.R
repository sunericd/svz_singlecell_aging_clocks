library(tidyverse)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/05_visualize")

chron_genes <- readRDS("data/clock_genes_all.rds")
bio_genes <- readRDS("data/clock_genes_all_bio.rds")

for (ct in unique(bio_genes$Celltype)){
  print(unname(as.data.frame(ct)),quote = FALSE, row.names = FALSE)
  ct_chron_genes <- chron_genes %>% filter(Celltype == ct)
  print(unname(as.data.frame(paste(dim(ct_chron_genes)[1], "chron. genes", sep=" "))),quote = FALSE, row.names = FALSE)
  ct_bio_genes <- bio_genes %>% filter(Celltype == ct)
  print(unname(as.data.frame(paste(dim(ct_bio_genes)[1], "bio. genes", sep=" "))),quote = FALSE, row.names = FALSE)
  print(unname(as.data.frame(paste(length(intersect(ct_chron_genes$gene, ct_bio_genes$gene)), "overlapping genes", sep=" "))),quote = FALSE, row.names = FALSE)
}