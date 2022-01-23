# Seurat Process Exercise Data
# Analysis Part 1
# Filter and Cluster

rm(list = ls())
library(Seurat)
library(tidyverse)
# R version 3.6.0
# Seurat_3.1.0.9007

#setwd("~/Dropbox/Exercise/Run2/") # local
setwd("/labs/abrunet1/Buckley/10.Exercise") # cluster

#tissues <- c("SVZ", "BM", "MU", "BL")
tissues <- "SVZ"

# Data Directory
setwd("1.CellRanger/all_filtered_feature_bc_matrices") # local

for (tissue in tissues) {
    print(tissue)

    files <- list.files(pattern=tissue, recursive = FALSE)

    object_list <- c()
    for (i in (1:length(files))) {
        print(files[i])
        data <- Read10X(paste0(files[i], "/filtered_feature_bc_matrix"))
        obj <- CreateSeuratObject(counts = data, project = files[i], min.cells = 10)
        object_list <- c(object_list, obj)
    }

    # Combine Objects
    obj <- Seurat:::merge.Seurat(x = object_list[[1]], y = object_list[c(2:length(object_list))])
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    #======================================================================================
    # Filter & Redo QC
    #======================================================================================

    # cut off of 1000 UMI instead of elbow from cell ranger
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
    obj <- subset(obj, subset = nCount_RNA > 1000 & percent.mt < 10 & nFeature_RNA > 500) 

    #======================================================================================
    # Dim reduction and clustering
    #======================================================================================

    # Run PCA, select PCs for tSNE visualization and graph-based clustering
    obj <- RunPCA(obj, verbose = FALSE)
    obj <- SCTransform(obj)
    obj <- FindNeighbors(obj, dims = 1:30)
    obj <- FindClusters(obj, resolution = 0.4)
    obj <- RunUMAP(obj, dims = 1:30)
    saveRDS(obj, paste0("data/seurat.", tissue, ".", Sys.Date(), ".rds"))
}





































