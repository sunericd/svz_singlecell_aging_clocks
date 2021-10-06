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
sct <- readRDS("data/multiseq_intergrated_seurat.rds")
sct.data <- as_tibble(t(as.matrix((sct[["SCT"]]@data))))# 17879 18220
genes <- colnames(sct.data) # use these genes to subset parabiosis SCT data
meta <- sct[[]]
table(meta$Celltype.LowRes)
rm(sct.data)
gc()

# Load Parabiosis 10x Data (v3 and v2 Combined)
pb <- readRDS("data/pb_combined_2.rds") # pb-svz_annotatedSeuratObject_2020-08-10.rds, pb_combined_2.rds

meta <- pb[[]]
table(meta$Celltype)

rm(para.missing)
rm(missing_df)
rm(pb.sct.data)
rm(para.genes)
gc()

