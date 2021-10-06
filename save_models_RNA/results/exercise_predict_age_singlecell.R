# Script for predicting age based on the ensemble 
# models saved in "./models/".
# Creates distribution plots and reports mean and 
# rank-sum test on mean predicted ages in Exercise
# groups. This should be run on SCG cluster.

# load packages
library(Seurat)
library(tidyverse)
library(Matrix)
library(Matrix.utils)
library(glmnet)
library(caret) 
library(dplyr)
library(ggplot2)
library(Metrics)
library(stringr)

#setwd("C:/Users/edsun/Desktop/RESEARCH/Brunet")
rsq <- function (x, y) cor(x, y, use="na.or.complete") ^ 2

normalization = TRUE # TRUE or FALSE (normalize by library size after summing)
N = 15 # number per cluster

# Load training data
sct.data <- as_tibble(t(as.matrix((readRDS("../data/multi_intergrated_seurat_Dec2020.rds")[["RNA"]]@counts))))# 17879 18220
genes <- colnames(sct.data) # use these genes to subset parabiosis SCT data
rm(sct.data)
gc()

# Load Exercise Data
pb_all <- readRDS("../data/joint_harmony_SVZ_2020-11-30.rds") 
pb <- subset(pb_all, AgeCond == "O_Exercise")
print(unique(pb_all[[]][, "AgeCond"]))

# get genes
pb.sct.data <- as_tibble(t(as.matrix((pb[["RNA"]]@counts)))) # 25531 19103
para.genes <- colnames(pb.sct.data)

# Fix parabiosis data
para.missing <- setdiff(genes, para.genes)
missing_df <- matrix(0, nrow(pb.sct.data), length(para.missing)) # fill in missing SCT gene with zeros.
colnames(missing_df) <- para.missing
rownames(missing_df) <- rownames(pb.sct.data)
pb.sct.data2 <- cbind(pb.sct.data, missing_df)[, genes] # Add, reorder, and cut to size.

rm(para.missing)
rm(missing_df)
rm(pb.sct.data)
rm(para.genes)
gc()

# Run pseudo-bulking on parabiosis data
meta <- pb[[]] # Object Metadata. Age, Celltype, and Celltype.LowRes most important
uniq_celltypes <- c("Oligodendro","Microglia","aNSC_NPC","Astrocyte_qNSC","Neuroblast","Endothelial")
uniq_AgeCond <- c("O_Control","O_Exercise","Y_Control" ,"Y_Exercise")

# Run models on parabiosis data

cidx = 1
for (ct in uniq_celltypes) {
    for (ac in uniq_AgeCond){
        # chunk and fix data
        pb <- subset(pb_all, AgeCond == ac)
        meta <- pb[[]]
        pb.sct.data <- as_tibble(t(as.matrix((pb[["RNA"]]@counts)))) # 25531 19103
        rm(pb)
        gc()
        para.genes <- colnames(pb.sct.data)
        para.missing <- setdiff(genes, para.genes)
        missing_df <- matrix(0, nrow(pb.sct.data), length(para.missing)) # fill in missing SCT gene with zeros.
        colnames(missing_df) <- para.missing
        rownames(missing_df) <- rownames(pb.sct.data)
        pb.sct.data2 <- cbind(pb.sct.data, missing_df)[, genes] # Add, reorder, and cut to size.
        rm(para.missing)
        rm(missing_df)
        rm(pb.sct.data)
        rm(para.genes)
        gc()
    
        clustered_data <- pb.sct.data2[meta$Celltype.LowRes == ct, ]
        submeta <- meta[meta$Celltype.LowRes == ct, ]
        #subumap <- pb[["umap"]]@cell.embeddings[meta$Celltype == ct, ]
        preds <- rep(0,dim(clustered_data)[1])
        # normalize data
        if (normalization == TRUE){
          #clustered_data = t(scale(t(clustered_data), center=FALSE, scale=colSums(t(clustered_data))*1e-4))
          clustered_data <- sweep(clustered_data, MARGIN = 1, FUN = "/", STATS = rowSums(clustered_data))
          clustered_data <- log1p(clustered_data * 10000)
        }
        fnlist <- list.files("./models/")
        # replace colnames to match model's data
        colnames(clustered_data) <- rep("",dim(clustered_data)[2])
        psdcl_ids <- c(1:dim(clustered_data)[1])
        model_ct = 0
        for (fn in fnlist){
          # retrieve appropriate models
          model_key <- str_split(fn, "_")[[1]][2]
          # run age prediction
          if (ct == "Oligo"){
              ct2 = "Oligodendro"
          }
          else {
              ct2 = ct
          }
          if (grepl(model_key,ct2) & !grepl("fraction",fn)){
            fit = readRDS(paste("./models/",fn,sep=""))
            pred = predict(fit, clustered_data)
            preds <- preds+pred
            model_ct <- model_ct + 1
          }
        }
        predicted_age <- preds/model_ct
        celltypes <- submeta$Celltype.LowRes
        actual_age <- submeta$Age
        condition <- submeta$AgeCond
        hash_ids_ct <- submeta$orig.ident
        #batches <- submeta$Experiment
        #umap1 <- subumap[,1]
        #umap2 <- subumap[,2]
        # get results
        if (cidx == 1){
          results <- cbind.data.frame(celltypes, hash_ids_ct, condition, predicted_age)
        }
        else {
          result <- cbind.data.frame(celltypes, hash_ids_ct, condition, predicted_age)
          results <- rbind.data.frame(results, result)
        }
        cidx <- cidx+1
        rm(clustered_data)
        rm(submeta)
        rm(subumap)
        gc()
    }
}

# save results file
results.df <- as.data.frame(results)
colnames(results.df) <- c('celltype','hash_id','condition','predicted_age')
write.csv(results.df,"./exercise_predicted_age_singlecell.csv",row.names=FALSE)

