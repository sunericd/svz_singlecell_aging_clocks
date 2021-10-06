#This script runs through training-validation of a ML method
#at a give value of N (number of cells per pseudocluster) and
#reports the R2/MAE of mean predicted ages with actual age of animal
#and saves a file with the mean predicted ages by cell type.


library(Seurat)
library(tidyverse)
library(Matrix)
library(Matrix.utils)
library(glmnet)
library(caret) 
library(dplyr)
library(ggplot2)
library(Metrics)

fileid = "Mural"
n_repeat = 20
clustering.method = "Random" # can be "Equal", "Random"
predictive.model = "glmnet" # can be any caret method: glmnet, knn, lssvmPoly, ranger, mlp, dnn
normalization = TRUE # TRUE or FALSE (normalize by library size after summing)
N = 15 # number per cluster
test.method = "LeaveOneBatch" # "Random" (80-20 split) or "LeaveOne" (take one out and avg)

#======================================================================================
#setwd("C:/Users/edsun/Desktop/RESEARCH/Brunet")

kmvar <- function(mat, clsize=10, method=c('random','maxd', 'mind', 'elki')){
  k = ceiling(nrow(mat)/clsize)
  km.o = kmeans(mat, k)#, algorithm=("Lloyd"))
  labs = rep(NA, nrow(mat))
  centd = lapply(1:k, function(kk){
    euc = t(mat)-km.o$centers[kk,]
    sqrt(apply(euc, 2, function(x) sum(x^2)))
  })
  centd = matrix(unlist(centd), ncol=k)
  clsizes = rep(0, k)
  if(method[1]=='random'){
    ptord = sample.int(nrow(mat))
  } else if(method[1]=='elki'){
    ptord = order(apply(centd, 1, min) - apply(centd, 1, max))
  } else if(method[1]=='maxd'){
    ptord = order(-apply(centd, 1, max))
  } else if(method[1]=='mind'){
    ptord = order(apply(centd, 1, min))
  } else {
    stop('unknown method')
  }
  for(ii in ptord){
    bestcl = which.max(centd[ii,])
    labs[ii] = bestcl
    clsizes[bestcl] = clsizes[bestcl] + 1
    if(clsizes[bestcl] >= clsize){
      centd[,bestcl] = NA
    }
  }
  return(labs)
}


# Load seurat object, contains all data
svz <- readRDS("../data/multi_intergrated_seurat_Dec2020.rds")

# retrieve info
meta <- svz[[]] # Object Metadata. Age, Celltype, and Celltype.LowRes most important
rawCounts <- svz[["RNA"]]@counts #svz[["SCT"]]@data # Raw UMI gene counts as sparse matrix
rm(svz)
gc()

# clustering by random pairing
uniq_celltypes <- rep(fileid,n_repeat)
uniq_ages <- unique(meta[, "Age"])
uniq_ids <- unique(meta[, "hash.ID"])
uniq_batchs <- unique(meta[, "orig.ident"])


# here we will iterate over repetitions of the same ct so each iteration is a new model
model_num <- 1
for (ct in uniq_celltypes) { # look at each cell type
  data <- t(rawCounts)[meta$Celltype.LowRes == ct, ]
  submeta <- meta[meta$Celltype.LowRes == ct, ]
  #uniq_ages <- unique(submeta1[, "Age"])
  
  mean_preds <- c() # only for LeaveOne
  preds <- c() # predicted ages
  actuals <- c() # actual ages
  hash_ids_labels <- c() # ids corresponding to predictions
  mae_vals <- c() # only for LeaveOne
  
  clustered_data <- c()
  ages <- c()
  hash_ids <- c()
  batch_ids <- c()
  
  # cluster within each age
  for (h in uniq_ids) { 
    
    randdata <- data[submeta$hash.ID == h, ]
    a <- mean(submeta[submeta$hash.ID == h, ]$Age)
    b <- submeta[submeta$hash.ID == h, ]$orig.ident[1]
    
    if (sum(submeta$hash.ID == h) > 1 & N > 1) {
      
      if (clustering.method == "Equal"){
        # PCA dim reduction
        pca <- prcomp(randdata, retx = TRUE) 
        
        # equal size clustering using k-means
        labels <- kmvar(pca$x[,1:2],N) # use first 2 PCs to cluster
        rm(pca)
        gc()
        
        # aggregate across cluster groups
        randdata <- aggregate.Matrix(randdata, 
                                     groupings = labels, fun = "sum")
        ages <- c(ages, rep(a,dim(randdata)[1]))
        hash_ids <- c(hash_ids, rep(h,dim(randdata)[1]))
        batch_ids <- c(batch_ids, rep(b,dim(randdata)[1]))
      }
      else if (clustering.method == "Random"){
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
        batch_ids <- c(batch_ids, rep(b,dim(randdata)[1]))
      }
      else {
        print ("clustering.method not defined!")
      }
    }
    # for N =1
    else if (sum(submeta$Age == a) > 1 & N == 1){
      ages <- c(ages, rep(a,dim(randdata)[1]))
      hash_ids <- c(hash_ids, rep(h,dim(randdata)[1]))
      batch_ids <- c(batch_ids, rep(b,dim(randdata)[1]))
    }
    # if only one sample
    else if (sum(submeta$Age == a) == 1){
      ages <- c(ages, a) # will be overwritten if more than one sample
      hash_ids <- c(hash_ids, h)
      batch_ids <- c(batch_ids, b)
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
  # separate train test
  if (test.method != "LeaveOneBatch"){
      for (hid in uniq_ids){
        if (test.method == "Random"){ # DEFUNCT FOR THIS SCRIPT
          clustered_data1 <- cbind(ages, clustered_data)
          data1 = sort(sample(nrow(clustered_data1), nrow(clustered_data1)*.8)) # 80-20 train-test split
          train_data<-clustered_data1[data1,-1]
          test_data<-clustered_data1[-data1,-1]
          train_labels<-clustered_data1[data1,1]
          test_labels<-clustered_data1[-data1,1]
        }
        else if (test.method == "LeaveOne"){
          clustered_data1 <- cbind(ages, clustered_data)
          train_data<-clustered_data1[hash_ids != hid,-1]
          test_data<-clustered_data1[hash_ids == hid,-1]
          train_labels<-clustered_data1[hash_ids != hid,1]
          test_labels<-clustered_data1[hash_ids == hid,1]
          train_ids<-hash_ids[hash_ids != hid]
          test_ids<-hash_ids[hash_ids == hid]
        }
        else {
          train_data<-clustered_data
          train_labels<-ages
        }
        # append number of eff samples
        nsamps <- dim(clustered_data)[1]
        
        # run age prediction
        fit = train(train_data, train_labels, method = predictive.model)
        if (length(test_labels)==1){ # to make dims match if only one sample
          pred = predict(fit, t(as.matrix(test_data)))
        }
        else {
          pred = predict(fit, test_data)
        }
        
        # compute stats
        rsq <- function (x, y) cor(x, y, use="na.or.complete") ^ 2
        preds <- c(preds, pred)
        actuals <- c(actuals,test_labels)
        hash_ids_labels <- c(hash_ids_labels,test_ids)
        mean_preds <- c(mean_preds, mean(pred))
        mae_vals <- c(mae_vals, mae(test_labels,pred))
        
        # clear
        #rm(train_data)
        #rm(test_data)
        #rm(train_labels)
        #rm(test_labels)
        #rm(fit)
        #gc()
      }
  }
  else {
    for (bid in uniq_batchs){
        clustered_data1 <- cbind(ages, clustered_data)
        train_data<-clustered_data1[batch_ids != bid,-1]
        test_data<-clustered_data1[batch_ids == bid,-1]
        train_labels<-clustered_data1[batch_ids != bid,1]
        test_labels<-clustered_data1[batch_ids == bid,1]
        train_ids<-batch_ids[batch_ids != bid]
        test_ids<-batch_ids[batch_ids == bid]
        
        # append number of eff samples
        nsamps <- dim(clustered_data)[1]
        
        # run age prediction
        fit = train(train_data, train_labels, method = predictive.model)
        if (length(test_labels)==1){ # to make dims match if only one sample
          pred = predict(fit, t(as.matrix(test_data)))
        }
        else {
          pred = predict(fit, test_data)
        }
        
        # compute stats
        rsq <- function (x, y) cor(x, y, use="na.or.complete") ^ 2
        preds <- c(preds, pred)
        actuals <- c(actuals,test_labels)
        hash_ids_labels <- c(hash_ids_labels,test_ids)
        mean_preds <- c(mean_preds, mean(pred))
        mae_vals <- c(mae_vals, mae(test_labels,pred))
    }
  }
  # clear
  rm(data)
  rm(submeta)
  rm(clustered_data)
  gc()
  
  print (ct)
  print (rsq(actuals,preds))
  print (mae(actuals,preds))
  
  # append results
  if (model_num==1){
    results <- cbind(rep(ct,length(preds)), hash_ids_labels, rep('normal',length(preds)), rep(model_num, length(preds)), actuals, preds)
  }
  else{
    result <- cbind(rep(ct,length(preds)), hash_ids_labels, rep('normal',length(preds)), rep(model_num, length(preds)), actuals, preds)
    results <- rbind(results, result)
  }
  model_num <- model_num + 1
}

# save final results for LeaveOne
results.df <- as.data.frame(results)
colnames(results.df) <- c('celltype','hash_id','condition','model_num','actual_age','predicted_age')
filename <- paste("./results/", fileid, test.method, clustering.method, predictive.model, N, "LeaveOneBatchOut", sep="_")
saveRDS(results.df,paste(filename,".rds",sep=""))