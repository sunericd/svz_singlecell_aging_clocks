
library(tidyverse)
library(Seurat)
library(glmnet)
library(glmnetUtils)
library(cowplot)
library(ggpubr)
set.seed(7)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/04_bulk")
d <- readRDS("data/celltype_specific_bulk_leave1batch_results.rds")



sumStats <- c()
for (cell in unique(d$Celltype)) {
    print(cell)
    cell_df <- filter(d, Celltype == cell)
    mae <- round(median(abs(cell_df$Age - cell_df$Prediction)), 3)
    meanae <- round(mean(abs(cell_df$Age - cell_df$Prediction)), 3)
    rho <- round(cor(cell_df$Age, cell_df$Prediction, method = "spearman"), 3)
    r2 <- round(cor(cell_df$Age, cell_df$Prediction, method = "pearson") ** 2, 3)
    print(paste0("MedianAbsErr : ", mae))
    print(paste0("MeanAbsErr : ", meanae))
    print(paste0("r2 : ", r2))
    print(paste0("spearman's rho : ", rho))
    print("------------------------- ")
    sumStats <- rbind(sumStats, c(cell, mae, meanae, rho, r2))
}

sumStats <- as.data.frame(sumStats)
colnames(sumStats) <- c("Celltype", "MedianAbsErr", "MeanAbsErr", "rho", "r2")
# sumStats
#         Celltype MedianAbsErr MeanAbsErr   rho    r2
# 1      Microglia        3.172      3.576 0.883 0.757
# 2    Endothelial         5.57      7.385 0.201 0.015
# 3    Oligodendro        2.936      3.014 0.872 0.803
# 4 Astrocyte_qNSC        3.692      4.935 0.676 0.495
# 5       aNSC_NPC         6.94      7.173 0.142 0.039
# 6     Neuroblast        7.997      8.334 0.019 0.004



