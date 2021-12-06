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

# Run pseudo-bulking on parabiosis data
meta <- pb[[]] # Object Metadata. Age, Celltype, and Celltype.LowRes most important
uniq_celltypes <- c("Oligodendro","Microglia","aNSC_NPC","Astrocyte_qNSC","Neuroblast","Endothelial")
#uniq_celltypes <- c("Astrocyte_qNSC")
uniq_AgeCond <- unique(na.omit(meta[, "AgeCond"]))
print(uniq_AgeCond)
uniq_ages <- unique(meta[, "Months"])

# Run models on parabiosis data
cidx = 1
celltypes <- c()
predicted_age <- c()
actual_age <- c()
condition <- c()
pseudocell_id <- c()
hash_ids_ct <- c()
batches <- c()

for (ct in uniq_celltypes) { # look at each cell type
  print (ct)
  for (ac in uniq_AgeCond){
    data <- pb.sct.data2[meta$Celltype == ct & meta$AgeCond == ac, ]
    submeta <- meta[meta$Celltype == ct & meta$AgeCond == ac, ]
    uniq_ids <- unique(na.omit(submeta[, "hash.ID"]))

    preds <- c() # predicted ages
    actuals <- c() # actual ages
    psdcell_id <- c()
    hash_ids_labels <- c()
    batchs <- c()

    clustered_data <- c()
    ages <- c()
    hash_ids <- c()
    batch <- c()

    # cluster within each age
    for (h in uniq_ids) {
      a <- mean(submeta[submeta$hash.ID == h, ]$Months)
      b <- submeta[submeta$hash.ID == h, ]$Experiment
      b <- b[1]
      randdata <- data[submeta$hash.ID == h, ]
      #print(dim(submeta))
      #print(dim(submeta[submeta$hash.ID == h, ]))
      #print(h)
      if (dim(submeta[submeta$hash.ID == h, ])[1] > 1 & N > 1) {
        # shuffle cells
        rows <- sample(nrow(randdata))
        randdata <- randdata[rows, ]
        
        # aggregate values by sliding frame
        idxs <- c()
        for (id in 1:ceiling(dim(randdata)[1]/N)) {
          idxs <- c(idxs, rep(id,N))
        }
        idxs <- idxs[1:dim(randdata)[1]]
        randdata <- rowsum(t(apply(randdata,1,as.numeric)),idxs)
        ages <- c(ages, rep(a,dim(randdata)[1]))
        hash_ids <- c(hash_ids, rep(h,dim(randdata)[1]))
        batch <- c(batch, rep(b,dim(randdata)[1]))
      }
      # for N =1
      else if (dim(submeta[submeta$hash.ID == h, ])[1] > 1 & N == 1){
        ages <- c(ages, rep(a,dim(randdata)[1]))
        hash_ids <- c(hash_ids, rep(h,dim(randdata)[1]))
        batch <- c(batch, rep(b,dim(randdata)[1]))
      }
      # if only one sample
      else if (dim(submeta[submeta$hash.ID == h, ])[1] == 1){
        ages <- c(ages, a) # will be overwritten if more than one sample
        hash_ids <- c(hash_ids, h)
        batch <- c(batch, b)
      }
      clustered_data <- rbind(clustered_data, randdata)
      rm(randdata)
      gc()
    }
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
        #print(dim(clustered_data))
        #print(colnames(clustered_data))
        #print(dim(fit$trainingData))
        #print(colnames(fit$trainingData))
        pred = predict(fit, clustered_data)
        preds <- c(preds, pred)
        actuals <- c(actuals,ages)
        hash_ids_labels <- c(hash_ids_labels, hash_ids)
        psdcell_id <- c(psdcell_id,psdcl_ids)
        batchs <- c(batchs,batch)
        # clear
        rm(fit)
        gc()
      }
    }
    # clear
    rm(data)
    rm(submeta)
    rm(clustered_data)
    gc()
    
    # print stats
    #print (ac)
    #print (actuals)
    #print(preds)
    print (rsq(actuals,preds))
    # append to cell-type results
    predicted_age <- c(predicted_age,preds)
    actual_age <- c(actual_age,actuals)
    condition <- c(condition,rep(ac,length(preds)))
    pseudocell_id <- c(pseudocell_id,psdcell_id)
    celltypes <- c(celltypes,rep(ct,length(preds)))
    hash_ids_ct <- c(hash_ids_ct,hash_ids_labels)
    batches <- c(batches, batchs)
  }
  # get results
  if (cidx == 1){
      results <- cbind(celltypes, hash_ids_ct, batches, condition, actual_age, predicted_age)
  }
  else {
      result <- cbind(celltypes, hash_ids_ct, batches, condition, actual_age, predicted_age)
      results <- rbind(results, result)
  }
  cidx <- cidx+1
}

# save results file
results.df <- as.data.frame(results)
print(celltypes[1:10])
print(hash_ids_ct[1:10])
print(batches[1:10])
print(condition[1:10])
print(actual_age[1:10])
print(predicted_age[1:10])
colnames(results.df) <- c('celltype','hash_id','batch','condition','actual_age','predicted_age')
write.csv(results.df,"./results/parabiosis_predicted_age.csv",row.names=FALSE)

