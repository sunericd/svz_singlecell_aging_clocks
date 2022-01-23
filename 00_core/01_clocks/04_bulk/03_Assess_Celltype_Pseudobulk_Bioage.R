
library(tidyverse)
library(Seurat)
library(glmnet)
library(glmnetUtils)
library(cowplot)
library(ggpubr)
set.seed(7)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/04_bulk")

d <- readRDS("data/celltype_specific_bulk_leave1batch_results_bioage.rds")


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
sumStats

        Celltype MedianAbsErr MeanAbsErr    rho    r2
1      Microglia        0.031      0.046  0.569 0.376
2    Endothelial        0.056      0.065  0.119 0.015
3    Oligodendro        0.039      0.042  0.594 0.467
4 Astrocyte_qNSC        0.035       0.04  0.816  0.62
5       aNSC_NPC        0.059      0.064 -0.155     0
6     Neuroblast        0.049      0.062 -0.244 0.064




