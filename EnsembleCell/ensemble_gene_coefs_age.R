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

# Load training data
sct.data <- as_tibble(t(as.matrix((readRDS("../data/multi_intergrated_seurat_Dec2020.rds")[["RNA"]]@counts))))# 17879 18220
genes <- colnames(sct.data) # use these genes to subset parabiosis SCT data
rm(sct.data)
gc()

# Load Parabiosis 10x Data (v3 and v2 Combined)
pb <- readRDS("../data/pb_combined_2.rds") # pb-svz_annotatedSeuratObject_2020-08-10.rds, pb_combined_2.rds

# Fix metadata
pb$AgeCond <- plyr::mapvalues(pb$AgeCond, 
          from=c("yy", "oo", "het_y", "het_o"),
          to=c("Young-Iso", "Old-Iso", "Young-Het", "Old-Het"))
Mouse <- pb@meta.data$hash.ID
Mouse[is.na(Mouse)] <- paste0(pb@meta.data$Experiment[is.na(Mouse)], "_",
                              pb@meta.data$AgeCond[is.na(Mouse)])
pb@meta.data$hash.ID <- Mouse

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

gene.names <- colnames(pb.sct.data2)

# Run pseudo-bulking on parabiosis data
meta <- pb[[]] # Object Metadata. Age, Celltype, and Celltype.LowRes most important
uniq_celltypes <- c("Oligodendro","Microglia","aNSC_NPC","Astrocyte_qNSC","Neuroblast","Endothelial")
uniq_AgeCond <- unique(na.omit(meta[, "AgeCond"]))
uniq_ages <- unique(meta[, "Months"])

# Run models on parabiosis data

cidx = 1
for (ct in uniq_celltypes) {
    fnlist <- list.files("./models/")
    coefs <- rep(0,length(gene.names))
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
        if (model_ct == 0){
          coefs <- coef(fit$finalModel, fit$bestTune$lambda)[2:length(coef(fit$finalModel, fit$bestTune$lambda))] # [1] is the intercept
        }
        else{
          coef = coef(fit$finalModel, fit$bestTune$lambda)[2:length(coef(fit$finalModel, fit$bestTune$lambda))] # [1] is the intercept
          coefs <- coefs+coef
        }
        model_ct <- model_ct + 1
      }
    }
    norm_coefs <- coefs/model_ct
    # get results
    if (cidx == 1){
      results <- norm_coefs
    }
    else {
      results <- cbind.data.frame(results, norm_coefs)
    }
    cidx <- cidx+1
}

# save results file
results <- cbind.data.frame(gene.names, results)
results.df <- as.data.frame(results)
colnames(results.df) <- c('gene.name',uniq_celltypes)
write.csv(results.df,"./results/ensemble_gene_coefs_age.csv",row.names=FALSE)

