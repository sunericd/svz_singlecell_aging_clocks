# Script for predicting age based on the ensemble 
# models saved in "./models/".
# Creates distribution plots and reports mean and 
# rank-sum test on mean predicted ages in parabiosis
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

rsq <- function (x, y) cor(x, y, use="na.or.complete") ^ 2

normalization = TRUE # TRUE or FALSE (normalize by library size after summing)
N = 15 # number per cluster

# Load training data
sct <- readRDS("../data/multi_intergrated_seurat_Dec2020.rds")
sct.data <- as_tibble(t(as.matrix((sct[["RNA"]]@counts))))# 17879 18220
print(dim(sct.data))

# Run pseudo-bulking on parabiosis data
meta <- sct[[]] # Object Metadata. Age, Celltype, and Celltype.LowRes most important
uniq_celltypes <- c("Oligodendro","Microglia","aNSC_NPC","Astrocyte_qNSC","Neuroblast","Endothelial")
uniq_ages <- unique(meta[, "Age"])


cidx = 1
for (ct in uniq_celltypes) {
    clustered_data <- sct.data[meta$Celltype.LowRes == ct, ]
    submeta <- meta[meta$Celltype.LowRes == ct, ]
    subumap <- sct[["umap"]]@cell.embeddings[meta$Celltype.LowRes == ct, ]
    #print(dim(subumap))
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
    hash_ids_ct <- submeta$hash.ID
    umap1 <- subumap[,1]
    umap2 <- subumap[,2]
    #print(length(umap1))
    #print(length(umap2))
    # get results
    if (cidx == 1){
      results <- cbind.data.frame(celltypes, hash_ids_ct, actual_age, predicted_age, umap1, umap2)
    }
    else {
      result <- cbind.data.frame(celltypes, hash_ids_ct, actual_age, predicted_age, umap1, umap2)
      #print(dim(results))
      #print(dim(result))
      results <- rbind.data.frame(results, result)
    }
    cidx <- cidx+1
    rm(clustered_data)
    rm(submeta)
    gc()
}


# save results file
results.df <- as.data.frame(results)
colnames(results.df) <- c('celltype','hash_id','actual_age','predicted_age','umap1','umap2')
write.csv(results.df,"./results/baseline_predicted_age_singlecell.csv",row.names=FALSE)

# print out R2 for each celltype
for (ct in uniq_celltypes) {
  print(ct)
  df_ct <- results.df %>% filter(celltype==ct)
  df_ct <- df_ct %>% group_by(hash_id) %>% summarise_at(vars(predicted_age, actual_age),funs(mean(.,na.rm=TRUE)))
  print(rsq(df_ct$predicted_age,df_ct$actual_age))
}

