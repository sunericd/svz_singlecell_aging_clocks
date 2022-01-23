# Predict Bootstrap Pseudocells using leave-one-batch-out cross validation
library(tidyverse)
library(glmnet)
library(viridis)
library(ggpubr)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")
df <- readRDS("data/bootstrap_pseudocell_15.rds")
df <- df %>% ungroup %>% select(-c(data))

#   hash.ID            Age Celltype.LowRes orig.ident pseudocell_all         
#   <chr>            <dbl> <fct>           <chr>      <list>                 
# 1 Sample4-ACCAATGC  29   Oligodendro     Batch-1    <tibble [100 x 20,948]>
# 2 Sample4-ACCAATGC  29   Endothelial     Batch-1    <tibble [100 x 20,948]>
# 3 Sample4-ACCAATGC  29   Microglia       Batch-1    <tibble [100 x 20,948]>
# 4 Sample1-TGTGATGG   4.7 Oligodendro     Batch-1    <tibble [100 x 20,948]>
# 5 Sample3-CTCTAGAC  20.8 Microglia       Batch-1    <tibble [100 x 20,948]>
# 6 Sample2-TCAATGGC   6.7 Oligodendro     Batch-1    <tibble [100 x 20,948]>


# Lognormalize counts, including a pseudocount
lognorm <- function(input) {
    norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
    log1p(norm * 10000)
}

df2 <- df %>% mutate(lognorm = map(pseudocell_all, lognorm)) 
df2$pseudocell_all <- NULL
df2 <- unnest(df2, lognorm)
dim(df2) # 16800 20952
df <- NULL # Clear memory

#==================================================================================================
# Leave One Batch Out 
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
        test_df <- df_ss %>% select("hash.ID", "Age", "Celltype.LowRes",
                                    "orig.ident", "Predictions")
        test_df <- as.data.frame(test_df)
        colnames(test_df) <- c("Mouse", "Age", "Celltype", "Batch", "Prediction")
        # Update
        full_df <- rbind(full_df, test_df)
    }
}


dim(full_df)
saveRDS(full_df, "data/leaveOneBatchPredictions_bootstrapPseudocell.rds")
full_df <- readRDS("data/leaveOneBatchPredictions_bootstrapPseudocell.rds")

#==================================================================================================
# Viz and stats
#==================================================================================================

d <- as_tibble(full_df)

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

#         Celltype MedianAbsErr MeanAbsErr   rho    r2
# 1    Endothelial        4.646      5.473 0.587 0.375
# 2 Astrocyte_qNSC        3.644      4.208 0.837 0.708
# 3       aNSC_NPC        3.937      4.494 0.826 0.699
# 4     Neuroblast        5.071      5.091 0.778  0.53
# 5    Oligodendro        2.133      2.903 0.876 0.775
# 6      Microglia        2.633      3.183 0.884 0.768


#==================================================================================================
# Stats on the median values

d2 <- d %>% group_by(Celltype, Age) %>% mutate(med = median(Prediction)) %>% select(-Prediction) %>% distinct()
sumStats <- c()
for (cell in unique(d2$Celltype)) {
    print(cell)
    cell_df <- filter(d2, Celltype == cell)
    mae <- round(median(abs(cell_df$Age - cell_df$med)), 3)
    meanae <- round(mean(abs(cell_df$Age - cell_df$med)), 3)
    rho <- round(cor(cell_df$Age, cell_df$med, method = "spearman"), 3)
    r2 <- round(cor(cell_df$Age, cell_df$med, method = "pearson") ** 2, 3)
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

#         Celltype MedianAbsErr MeanAbsErr   rho    r2
# 1    Endothelial        4.805       5.24 0.688 0.499
# 2 Astrocyte_qNSC        3.302      4.072 0.904 0.829
# 3       aNSC_NPC        3.743      4.295 0.878 0.768
# 4     Neuroblast        5.424      4.908 0.863 0.632
# 5    Oligodendro        1.619      2.474 0.909 0.826
# 6      Microglia        2.072      2.757 0.916 0.841


celltype_order <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")
d$Celltype <- factor(d$Celltype, levels = celltype_order)
d2$Celltype <- factor(d2$Celltype, levels = celltype_order)

ggplot(d, aes(x = Age, y = Prediction)) +
    geom_bin2d(binwidth = 1.5) +
    scale_fill_distiller(palette="PuBuGn", direction = 1, limits = c(0, 70), oob = scales::squish) +
    facet_wrap(Celltype~.) +
    theme_classic() +
    stat_cor(label.y = 28, size = 2.5, color = "darkblue") +
    geom_point(data = d2, aes(x=Age, y=med), size = 1.5, color = "black") +
    geom_smooth(data = d2, aes(x=Age, y=med), method = "lm", color = "black", size = 1) +
    stat_cor(data = d2, aes(x=Age, y=med), label.y = 31, size = 2.7) # pearson r not rho
ggsave("plots/Figure1h_PredictedAges_BootstrapChronological.pdf", width = 6.36, height = 4.29, useDingbats = F)










