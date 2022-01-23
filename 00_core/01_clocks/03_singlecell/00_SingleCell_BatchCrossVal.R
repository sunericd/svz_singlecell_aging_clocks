# Pure Single Cell
# Batch Cross Validation

library(modelr)
library(tidyverse)
library(broom)
library(glmnet)
library(glmnetUtils)
library(Seurat)
library(viridis)
library(ggpubr)
library(ggthemes)

#============================================================================
# Load Seurat Object
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/03_singlecell")
svz <- readRDS("../../00_multiseq_aging/Integrate/data/multi_intergrated_seurat_Dec2020.rds")

#============================================================================
# Batch Cross Val
#===========================================================================

# Initialize main df
full_df_2 <- c()

batches <- unique(svz@meta.data$orig.ident)
Celltypes <- c("Endothelial", "Astrocyte_qNSC", "aNSC_NPC",
               "Neuroblast", "Oligodendro", "Microglia")
for (CT in Celltypes) {
    print(CT)
    svz_CT <- subset(svz, subset = Celltype.LowRes == CT)
    for (batch in batches) {
        print(batch)
        # Split object
        svz_fmo <- subset(svz_CT, subset = orig.ident == batch, invert = T)
        svz_ss <- subset(svz_CT, subset = orig.ident == batch, invert = F)

        # Convert to dataframe
        rna.data_fmo <- as_tibble(t(as.matrix((svz_fmo[["RNA"]]@data))))
        df_fmo <- tibble("Age"=as.numeric(svz_fmo@meta.data$Age), rna.data_fmo)
        rna.data_ss <- as_tibble(t(as.matrix((svz_ss[["RNA"]]@data))))
        df_ss <- tibble("Age"=as.numeric(svz_ss@meta.data$Age), rna.data_ss)

        # Train
        model <- cv.glmnet(x = as.matrix(df_fmo[,-1]), y = as.matrix(df_fmo[,1]),
                           type.measure="mae",
                           standardize=F,
                           relax=F,
                           nfolds = 5)

        # Predict
        testAges <- df_ss[,1]
        testPredictions <- predict(model, newx=as.matrix(df_ss[,-1]), s="lambda.min")
        test_df <- data.frame("Age" = testAges$Age, "Prediction" = testPredictions[,1])
        test_df$Celltype <- CT
        test_df$Batch <- batch

        # Update the outer world
        full_df_2 <- rbind(full_df_2, test_df)
    }
}

saveRDS(full_df_2, "data/singlecell_BatchCrossVal.rds")

#============================================================================
# Viz and get stats
#============================================================================

d <- full_df_2

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
# 1    Endothelial         5.93      6.711 0.027 0.001
# 2 Astrocyte_qNSC        4.678      5.269 0.541 0.277
# 3       aNSC_NPC        4.211      4.884 0.418 0.217
# 4     Neuroblast         4.82      5.375 0.272 0.081
# 5    Oligodendro        3.992       5.14 0.662 0.396
# 6      Microglia        5.306      6.155 0.475 0.201


# For the curious, using SCT instead of log-normalized RNA:
#         Celltype MedianAbsErr MeanAbsErr    rho    r2
# 1    Endothelial        6.216      6.832 -0.069 0.011
# 2 Astrocyte_qNSC         4.84      5.413   0.52 0.255
# 3       aNSC_NPC        4.434      5.108  0.329 0.148
# 4     Neuroblast        4.975      5.363  0.267 0.082
# 5    Oligodendro        3.955      5.204  0.645 0.362
# 6      Microglia        5.684      6.444  0.391 0.137 
