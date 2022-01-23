# LMO + 8 SVZ Seurat Processing 
# Batch 4

# VERSION 3.0 Seurat

rm(list = ls())
library(Seurat)
library(dplyr)
library(MAST)
library(tidyverse)
library(Matrix)
library(sctransform)
library(scales)
library(ggthemes)
library(viridis)

sessionInfo()


# All subsequent paths will be relative.
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Batch4_8SVZ_2020_11/3.Seurat")
dir.create("plots")
dir.create("data")

#======================================================================================
# Create Seurat object and basic processing
#======================================================================================

# Read in Cellranger count matrix and create basic seurat object
svz.data <- Read10X("../1.CellRanger/filtered_feature_bc_matrix") # 31053 12429
lmo.umi.data <- Read10X("../2.LMO/umi_count/", gene.column=1) # 10 12413

dim(lmo.umi.data)
lmo.umi.data[1:10,1:10]

# Remove negative control barcode and unmapped row:
lmo.data <- lmo.umi.data[1:8,]
columns_lmo <- paste0(colnames(lmo.data), "-1")
colnames(lmo.data) <- columns_lmo

barcode.intersect <- intersect(colnames(svz.data), columns_lmo) # 5263
svz.data <- svz.data[, barcode.intersect]
lmo.umi.data <- lmo.data[, barcode.intersect]

svz <- CreateSeuratObject(counts = svz.data)
svz <- NormalizeData(svz)
svz <- FindVariableFeatures(svz)
svz <- ScaleData(svz)


# QC Filter
svz[["percent.mt"]] <- PercentageFeatureSet(svz, pattern = "^mt-")
svz <- subset(svz, subset = nCount_RNA > 1000 & percent.mt < 10 & nFeature_RNA > 500) 

#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select 21 PCs for tSNE visualization and graph-based clustering
svz <- RunPCA(svz, verbose = FALSE)

svz <- FindNeighbors(svz, dims = 1:21)
svz <- FindClusters(svz, resolution = 0.15)
svz <- RunUMAP(svz, dims = 1:21)
DimPlot(svz)

#======================================================================================
# Add LMO sample label data to object
#======================================================================================
barcode.intersect <- intersect(colnames(svz), columns_lmo)
lmo.data <- lmo.data[ , barcode.intersect]
svz[["LMO"]] <- CreateAssayObject(counts = lmo.data)
svz <- NormalizeData(object = svz, assay = "LMO", normalization.method = "CLR")
svz <- HTODemux(svz, assay = "LMO", positive.quantile = 0.95, nsamples = 1000)
table(svz$LMO_classification.global)
 # Doublet Negative  Singlet 
 #     565      826     3626 

color_pal.8 <- tableau_color_pal(palette = "Tableau 10")(8)

# Group cells based on the max HTO signal
Idents(svz) <- "LMO_maxID"

# Heatmap 
HTOHeatmap(svz, assay = "LMO", ncells = 3000) + scale_fill_viridis(option="viridis", direction=1)

# Violin plots of expression by single/doublet/negative catagories
Idents(svz) <- "LMO_classification.global"
VlnPlot(svz, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

#======================================================================================
# Reductions based on LMO counts

Idents(svz) <- "LMO_classification.global"
svz.subset <- subset(svz, idents = "Negative", invert = TRUE)

# Set two colors
color_pal.2 <- rev(tableau_color_pal(palette = "Tableau 10")(2))

# Calculate a distance matrix using HTO
lmo.dist.mtx <- as.matrix(dist(t(GetAssayData(object = svz.subset, assay = "LMO"))))

# Calculate tSNE embeddings with a distance matrix
svz.subset <- RunTSNE(svz.subset, distance.matrix = lmo.dist.mtx, perplexity = 30)

# tSNE based on LMO reduction, colored by doublet/singlet status
DimPlot(svz.subset, cols = color_pal.2, reduction = "tsne")
ggsave("plots/tsne.lmo.doub.sing.pdf", width = 5.9, height = 4.6)

#======================================================================================
# Subset to singlets
#======================================================================================

# Set colors
color_pal.1 <- tableau_color_pal(palette = "Tableau 10")(1)
color_pal.4 <- tableau_color_pal(palette = "Tableau 10")(4)

# Remove doublets and Negatives
svz_sing <- subset(svz, subset = LMO_classification.global == "Singlet")

# Calculate a distance matrix using HTO
lmo.dist.mtx <- as.matrix(dist(t(GetAssayData(object = svz_sing, assay = "LMO"))))

# Calculate tSNE embeddings with a distance matrix
svz_sing <- RunTSNE(svz_sing, distance.matrix = lmo.dist.mtx, perplexity = 30)


# tSNE based on LMO reduction, colored by sample
DimPlot(svz_sing, cols = color_pal.8, reduction = "tsne", group.by="hash.ID")
ggsave("plots/tsne.lmo.sing.pdf", width = 5.9, height = 4.6)

#======================================================================================

# Group cells based on the max HTO signal
Idents(svz_sing) <- "LMO_maxID"
svz_sing <- NormalizeData(object = svz_sing, assay = "LMO", normalization.method = "RC")
RidgePlot(svz_sing, assay = "LMO",
    features = rownames(svz[["LMO"]])[1:8], ncol = 2, cols=color_pal.8)
# Heatmap 
svz_sing <- NormalizeData(object = svz_sing, assay = "LMO", normalization.method = "CLR")
Idents(svz_sing) <- "LMO_maxID"
HTOHeatmap(svz_sing, assay = "LMO", ncells = 3000) + scale_fill_viridis(option="viridis", direction=1)
ggsave("plots/heatmap.png", width = 6.9, height=2.88)

#======================================================================================
# Run normalization and clustering on fully filtered subset
#======================================================================================
# Normalize, Standardize
svz_sing <- SCTransform(svz_sing)

# Reduce
svz_sing <- RunPCA(svz_sing, verbose = FALSE)
svz_sing <- FindNeighbors(svz_sing, dims = 1:20)
svz_sing <- FindClusters(svz_sing, resolution = 0.25)
svz_sing <- RunUMAP(svz_sing, dims = 1:20)

DimPlot(svz_sing, group.by = "hash.ID", pt.size = .6, cols = color_pal.8)
ggsave("plots/umap.select.sct.pdf", width = 5.9, height = 4.6)

svz_sing.markers <- FindAllMarkers(object=svz_sing)

#saveRDS(svz_sing.markers, paste0("data/svz_sing.markers.", Sys.Date(), ".rds"))
#saveRDS(svz_sing, paste0("data/svz_sing_", Sys.Date(), ".rds"))

# More analysis in seurat_batch4_part2.R

# Fin.