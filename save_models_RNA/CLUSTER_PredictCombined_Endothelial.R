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

fileid = "Endothelial"
n_repeat = 20
clustering.method = "Random" # can be "Equal", "Random"
predictive.model = "glmnet" # can be any caret method: glmnet, knn, lssvmPoly, ranger, mlp, dnn
normalization = TRUE # TRUE or FALSE (normalize by library size after summing)
N = 15 # number per cluster
test.method = "LeaveOne" # "Random" (80-20 split) or "LeaveOne" (take one out and avg)


#======================================================================================
#setwd("C:/Users/edsun/Desktop/RESEARCH/Brunet")
prolif_frac_fn <- "proliferative_fractions.csv"
rsq <- function (x, y) cor(x, y, use="na.or.complete") ^ 2

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
svz <- readRDS("data/multi_intergrated_seurat_Dec2020.rds")

# retrieve info
meta <- svz[[]] # Object Metadata. Age, Celltype, and Celltype.LowRes most important
rawCounts <- svz[["RNA"]]@counts # Raw UMI gene counts as sparse matrix
rm(svz)
gc()

# clustering by random pairing
uniq_celltypes <- rep(fileid,n_repeat)
uniq_ages <- unique(meta[, "Age"])

# get prolif fractions
#fracvec <- meta$Prolif_Lineage_Fraction_of_SVZ
#if (test.method == "LeaveOne"){
#  results <- c(fracvec)
#}

count_idx <- 0
for (ct in uniq_celltypes) { # look at each cell type
  data <- t(rawCounts)[meta$Celltype.LowRes == ct, ]
  submeta <- meta[meta$Celltype.LowRes == ct, ]
  #uniq_ages <- unique(submeta1[, "Age"])
  
  mean_preds <- c() # only for LeaveOne
  preds <- c() # predicted ages
  actuals <- c() # actual ages
  mae_vals <- c() # only for LeaveOne
  
  clustered_data <- c()
  ages <- c()
  fracs <- c()
  
  # cluster within each age
  idx = 1
  for (a in uniq_ages) { 
    
    randdata <- data[submeta$Age == a, ]
    f <- mean(submeta[submeta$Age == a, ]$Prolif_Lineage_Fraction_of_SVZ)
    
    if (sum(submeta$Age == a) > 1 & N > 1) {
      
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
        fracs <- c(fracs, rep(f,dim(randdata)[1]))
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
        fracs <- c(fracs, rep(f,dim(randdata)[1]))
        ages <- c(ages, rep(a,dim(randdata)[1]))
      }
      else {
        print ("clustering.method not defined!")
      }
    }
    # for N =1
    else if (sum(submeta$Age == a) > 1 & N == 1){
      fracs <- c(fracs, rep(f,dim(randdata)[1]))
      ages <- c(ages, rep(a,dim(randdata)[1]))
    }
    # if only one sample
    else if (sum(submeta$Age == a) == 1){
      fracs <- c(fracs, f) # will be overwritten if more than one sample
      ages <- c(ages, a) # will be overwritten if more than one sample
    }
    clustered_data <- rbind(clustered_data, randdata)
    rm(randdata)
    gc()
    idx = idx + 1
  }
  # normalize data
  if (normalization == TRUE){
    #clustered_data = t(scale(t(clustered_data), center=FALSE, scale=colSums(t(clustered_data))*1e-4))
    clustered_data <- sweep(clustered_data, MARGIN = 1, FUN = "/", STATS = rowSums(clustered_data))
    clustered_data <- log1p(clustered_data * 10000)
  }
  # separate train test
  clustered_data1 <- cbind(fracs, clustered_data)
  train_data<-clustered_data1[,-1]
  train_labels<-clustered_data1[,1]
  
  # run age prediction
  fit = train(train_data, train_labels, method = predictive.model)
  pred = predict(fit, train_data)
  
  # clear
  rm(data)
  rm(submeta)
  rm(clustered_data)
  gc()
  
  print (ct)
  if (test.method == "LeaveOne"){
    print ("Mean Predicted R2 and MAE:")
    print (rsq(pred,train_labels))
    print (mae(pred,train_labels))
  }
  
  # save model
  count_idx <- count_idx + 1
  savename <- paste("./results/models/", fileid, test.method, clustering.method, predictive.model, N, count_idx, 'fraction', sep="_")
  saveRDS(fit, paste(savename,".rds",sep=""))
}
