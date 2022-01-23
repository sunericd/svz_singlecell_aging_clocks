# Predict Bootstrap Pseudocells using batch cross validation

library(tidyverse)
library(glmnet)
library(viridis)
library(ggpubr)
library(scales)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/01_bootstrap_bioage")
df <- readRDS("data/bootstrap_cell_15_bio.rds")
df <- df %>% ungroup


lognorm <- function(input) {
    norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
    log1p(norm * 10000)
}

df2 <- df %>% mutate(lognorm = map(pseudocell_all, lognorm))
df2$pseudocell_all <- NULL
df2 <- unnest(df2, lognorm)
df <- NULL # Clear memory



#==================================================================================================
# Leave One Batch Out Version
#==================================================================================================

# Initialize main df
full_df <- c()
batches <- unique(df2$orig.ident)
Celltypes <- c("Endothelial", "Astrocyte_qNSC", "aNSC_NPC",
               "Neuroblast", "Oligodendro", "Microglia")

for (CT in Celltypes) {
    print(CT)
    df_CT <- filter(df2, Celltype.LowRes == CT)
    for (batch in batches) {
        print(batch)
        # Split object
        df_fmo <- filter(df_CT, orig.ident != batch)
        df_ss <- filter(df_CT, orig.ident == batch)
        # Train
        model <- cv.glmnet(x = as.matrix(df_fmo[,-c(1:4)]),
                           y = as.matrix(df_fmo[,2]),
                           type.measure="mae",
                           standardize=F,
                           relax=F,
                           nfolds = 5)
        # Predict
        testPredictions <- predict(model, newx=as.matrix(df_ss[,-c(1:4)]), s="lambda.min")
        df_ss$Predictions <- testPredictions[,1]
        test_df <- df_ss %>% select("hash.ID", "Prolif_Lineage_Fraction_of_SVZ", "Celltype.LowRes",
                                    "orig.ident", "Predictions")
        test_df <- as.data.frame(test_df)
        colnames(test_df) <- c("Mouse", "Prolif_Lineage_Fraction_of_SVZ", "Celltype", "Batch", "Prediction")
        # Update
        full_df <- rbind(full_df, test_df)
    }
}

dim(full_df)
# Don't overwrite
#saveRDS(full_df, "data/leaveOneBatchPredictions_bootstrapPC_bio.rds")

full_df <- readRDS("data/leaveOneBatchPredictions_bootstrapPC_bio.rds")
#==================================================================================================
# Viz and stats
#==================================================================================================

d <- as_tibble(full_df)

sumStats <- c()
for (cell in unique(d$Celltype)) {
    print(cell)
    cell_df <- filter(d, Celltype == cell)
    mae <- round(median(abs(cell_df$Prolif_Lineage_Fraction_of_SVZ - cell_df$Prediction)), 3)
    meanae <- round(mean(abs(cell_df$Prolif_Lineage_Fraction_of_SVZ - cell_df$Prediction)), 3)
    rho <- round(cor(cell_df$Prolif_Lineage_Fraction_of_SVZ, cell_df$Prediction, method = "spearman"), 3)
    r2 <- round(cor(cell_df$Prolif_Lineage_Fraction_of_SVZ, cell_df$Prediction, method = "pearson") ** 2, 3)
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

        Celltype MedianAbsErr MeanAbsErr   rho    r2
1    Endothelial        0.049      0.053 0.399 0.149
2 Astrocyte_qNSC        0.034      0.038 0.748 0.591
3       aNSC_NPC        0.029      0.039 0.657 0.452
4     Neuroblast        0.036      0.045 0.562 0.334
5    Oligodendro        0.026       0.03 0.816 0.723
6      Microglia         0.03      0.037  0.76 0.559

# For comparison, leave one batch out stats for age:
        Celltype MedianAbsErr MeanAbsErr   rho    r2
1    Endothelial        4.646      5.473 0.587 0.375
2 Astrocyte_qNSC        3.644      4.208 0.837 0.708
3       aNSC_NPC        3.937      4.494 0.826 0.699
4     Neuroblast        5.071      5.091 0.778  0.53
5    Oligodendro        2.133      2.903 0.876 0.775
6      Microglia        2.633      3.183 0.884 0.768




#==================================================================================================
# Median and viz for median values

d2 <- d %>% group_by(Celltype, Prolif_Lineage_Fraction_of_SVZ) %>% mutate(med = median(Prediction)) %>% select(-Prediction) %>% distinct()

sumStats <- c()
for (cell in unique(d2$Celltype)) {
    print(cell)
    cell_df <- filter(d2, Celltype == cell)
    mae <- round(median(abs(cell_df$Prolif_Lineage_Fraction_of_SVZ - cell_df$med)), 3)
    meanae <- round(mean(abs(cell_df$Prolif_Lineage_Fraction_of_SVZ - cell_df$med)), 3)
    rho <- round(cor(cell_df$Prolif_Lineage_Fraction_of_SVZ, cell_df$med, method = "spearman"), 3)
    r2 <- round(cor(cell_df$Prolif_Lineage_Fraction_of_SVZ, cell_df$med, method = "pearson") ** 2, 3)
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

        Celltype MedianAbsErr MeanAbsErr   rho    r2
1    Endothelial        0.046      0.052 0.431 0.168
2 Astrocyte_qNSC        0.031      0.036 0.795 0.655
3       aNSC_NPC        0.028      0.038 0.667 0.508
4     Neuroblast        0.029      0.043 0.592 0.375
5    Oligodendro        0.029      0.028 0.867 0.788
6      Microglia        0.023      0.034 0.822 0.651

# For comparison, leave one batch out stats for age:
        Celltype MedianAbsErr MeanAbsErr   rho    r2
1    Endothelial        4.805       5.24 0.688 0.499
2 Astrocyte_qNSC        3.302      4.072 0.904 0.829
3       aNSC_NPC        3.743      4.295 0.878 0.768
4     Neuroblast        5.424      4.908 0.863 0.632
5    Oligodendro        1.619      2.474 0.909 0.826
6      Microglia        2.072      2.757 0.916 0.841


#==================================================================================================

score_func <- function(input) {
    return(35 - 100 * input)
}

d$Prolif_Lineage_Fraction_of_SVZ_Score <- score_func(d$Prolif_Lineage_Fraction_of_SVZ)
d$Prediction_Score <- score_func(d$Prediction)
d2$Prolif_Lineage_Fraction_of_SVZ_Score <- score_func(d2$Prolif_Lineage_Fraction_of_SVZ)
d2$med_Score <- score_func(d2$med)

sumStats <- c()
for (cell in unique(d2$Celltype)) {
    print(cell)
    cell_df <- filter(d2, Celltype == cell)
    mae <- round(median(abs(cell_df$Prolif_Lineage_Fraction_of_SVZ_Score - cell_df$med_Score)), 3)
    meanae <- round(mean(abs(cell_df$Prolif_Lineage_Fraction_of_SVZ_Score - cell_df$med_Score)), 3)
    r <- round(cor(cell_df$Prolif_Lineage_Fraction_of_SVZ_Score, cell_df$med_Score, method = "pearson"), 3)
    r2 <- round(cor(cell_df$Prolif_Lineage_Fraction_of_SVZ_Score, cell_df$med_Score, method = "pearson") ** 2, 3)
    print(paste0("MedianAbsErr : ", mae))
    print(paste0("MeanAbsErr : ", meanae))
    print(paste0("r2 : ", r2))
    print(paste0("spearman's r : ", r))
    print("------------------------- ")
    sumStats <- rbind(sumStats, c(cell, mae, meanae, r, r2))
}
sumStats <- as.data.frame(sumStats)
colnames(sumStats) <- c("Celltype", "MedianAbsErr", "MeanAbsErr", "r", "r2")
sumStats

#         Celltype MedianAbsErr MeanAbsErr     r    r2
# 1    Endothelial         4.61        5.2  0.41 0.168
# 2 Astrocyte_qNSC        3.065      3.617  0.81 0.655
# 3       aNSC_NPC        2.751      3.755 0.713 0.508
# 4     Neuroblast        2.865      4.312 0.612 0.375
# 5    Oligodendro        2.876      2.752 0.888 0.788
# 6      Microglia        2.299      3.448 0.807 0.651


celltype_order <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")
d$Celltype <- factor(d$Celltype, levels = celltype_order)
d2$Celltype <- factor(d2$Celltype, levels = celltype_order)

ggplot(d, aes(x = Prolif_Lineage_Fraction_of_SVZ_Score, y = Prediction_Score)) +
    geom_bin2d(binwidth = 1.5) +
    scale_fill_distiller(palette="PuBuGn", direction = 1, limits = c(0, 70), oob = scales::squish)  +
    facet_wrap(Celltype~.) +
    theme_classic() +
    geom_point(data = d2, aes(x=Prolif_Lineage_Fraction_of_SVZ_Score, y=med_Score), size = 1.5) +
    geom_smooth(data = d2, aes(x=Prolif_Lineage_Fraction_of_SVZ_Score, y=med_Score), method = "lm", color = "black", size = .7) +
    stat_cor(data = d2, aes(x=Prolif_Lineage_Fraction_of_SVZ_Score, y=med_Score, label = ..r.label..), label.y = 31, size = 3) +
    ylab("Predicted Bio Age Score") + xlab("Bio Age Score: 35 - (Proliferative %)")

ggsave("plots/batch_cross_val_PlusMedianStats_bio_Fig1i.pdf", width = 6.36, height = 4.29, useDingbats = F)

#==================================================================================================
library(Seurat)
svz <- readRDS("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Integrate/data/multi_intergrated_seurat_Dec2020.rds")
svz$BioAge <- score_func(svz@meta.data$Prolif_Lineage_Fraction_of_SVZ)
plot(x = svz$Age, y = svz$BioAge)
plot(x = svz$Age, y = svz$Prolif_Lineage_Fraction_of_SVZ)
