# LMO + 8 SVZ Seurat Processing
# Batch 2

# Seurat 3.0 
# R 3.6.3

library(Seurat)
library(dplyr)
library(MAST)
library(tidyverse)
library(Matrix)
library(sctransform)
library(scales)
library(ggthemes)
library(viridis)
# Package version info at end.

# All subsequent paths will be relative.
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Batch2_8SVZ_2019_12/3.Seurat")
dir.create("plots")
dir.create("data")

#======================================================================================
# Create Seurat object and basic processing
#======================================================================================

# Read in Cellranger count matrix and create basic seurat object
svz.data <- Read10X("../1.CellRanger/filtered_feature_bc_matrix")
lmo.umi.data <- Read10X("../2.LMO/umi_count/", gene.column=1)

lmo.data <- lmo.umi.data[1:8,]
columns_lmo <- paste0(colnames(lmo.data), "-1")
colnames(lmo.data) <- columns_lmo


barcode.intersect <- intersect(colnames(svz.data), columns_lmo) # 11964
svz.data <- svz.data[, barcode.intersect]
lmo.data <- lmo.data[, barcode.intersect]

svz <- CreateSeuratObject(counts = svz.data)
svz <- NormalizeData(svz)
svz <- FindVariableFeatures(svz)
svz <- ScaleData(svz)

# QC Filter
svz[["percent.mt"]] <- PercentageFeatureSet(svz, pattern = "^mt-")
svz <- subset(svz, subset = nCount_RNA > 1000 & percent.mt < 10 & nFeature_RNA > 500) # 11245


#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select 21 PCs for tSNE visualization and graph-based clustering
svz <- RunPCA(svz, verbose = FALSE)
ElbowPlot(svz, ndims = 50)
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
svz <- HTODemux(svz, assay = "LMO", positive.quantile = 0.99)
table(svz$LMO_classification.global)

color_pal.8 <- tableau_color_pal(palette = "Tableau 10")(8)

# Group cells based on the max HTO signal
Idents(svz) <- "LMO_maxID"

# Heatmap version of above ridgeplots (simple but less quant)
HTOHeatmap(svz, assay = "LMO", ncells = 3000) + scale_fill_viridis(option="viridis", direction=1)

# Violin plots of expression by single/doublet/negative catagories
Idents(svz) <- "LMO_classification.global"
VlnPlot(svz, features = "nCount_RNA", pt.size = 0.1, log = F)

#======================================================================================
# Reductions based on LMO counts
Idents(svz) <- "LMO_classification.global"
svz.subset <- subset(svz, idents = "Negative", invert = TRUE)

# Set two colors
color_pal.2 <- rev(tableau_color_pal(palette = "Tableau 10")(2))

# Calculate a distance matrix using HTO
lmo.dist.mtx <- as.matrix(dist(t(GetAssayData(object = svz.subset, assay = "LMO"))))

# Calculate tSNE embeddings with a distance matrix
svz.subset <- RunTSNE(svz.subset, distance.matrix = lmo.dist.mtx, perplexity = 50)

# tSNE based on LMO reduction, colored by doublet/singlet status
DimPlot(svz.subset, cols = color_pal.2, reduction = "tsne")
ggsave("plots/tsne.lmo.doub.sing.pdf", width = 5.9, height = 4.6)

# UMAP based on transcriptome reduction, colored by doublet/singlet status
DimPlot(svz.subset, cols = color_pal.2, reduction = "umap")
ggsave("plots/umap.doub.sing.pdf", width = 5.9, height = 4.6)


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
# Cell Selector
# Remove ambigous cells based on dimensionality reduction
#======================================================================================

plot <- DimPlot(svz_sing, cols = color_pal.8, reduction = "tsne", group.by="hash.ID")
cells.located <- CellSelector(plot = plot)
ambiguous <- cells.located # 464 Cells removed.
#svz_select <- subset(svz_sing, cells = ambiguous, invert = TRUE)

svz_select <- readRDS("data/svz_select_2020-04-12.rds")
DimPlot(svz_select, cols = color_pal.8, reduction = "tsne", group.by="hash.ID")
ggsave("plots/tsne.lmo.filtered.pdf", width = 5.9, height = 4.6)

#======================================================================================
# Group cells based on the max HTO signal
Idents(svz_select) <- "LMO_maxID"

svz_select <- NormalizeData(object = svz_select, assay = "LMO", normalization.method = "RC")
RidgePlot(svz_select, assay = "LMO", features = rownames(svz[["LMO"]])[1:8], ncol = 2, cols=color_pal.8)

# Heatmap version of above ridgep1 (simple but less quant)
svz_select <- NormalizeData(object = svz_select, assay = "LMO", normalization.method = "CLR")
Idents(svz_select) <- "LMO_maxID"
HTOHeatmap(svz_select, assay = "LMO", ncells = 3000) + scale_fill_viridis(option="viridis", direction=1)
ggsave("plots/heatmap.lmo.png", width = 6.9, height=2.88)

table(svz_select$LMO_maxID)
#  BC13-3.33-CAGTTAGG  BC14-5.40-AACCGAAC  BC15-9.47-AAGCAGTC BC16-14.50-GAATCAGG 
#                 498                1253                 836                 919 
# BC17-16.53-ACTCGAAG BC18-18.58-TTCCACGT BC19-20.60-AACTGCAG BC20-22.57-TAAGCAGC 
#                 131                 288                 116                 177 

#======================================================================================
# Run normalization and clustering on fully filtered subset
#======================================================================================

svz_select <- SCTransform(svz_select)

# Reduce
svz_select <- RunPCA(svz_select, verbose = FALSE)
svz_select <- FindNeighbors(svz_select, dims = 1:20)
svz_select <- FindClusters(svz_select, resolution = 0.5)
svz_select <- RunUMAP(svz_select, dims = 1:20)

DimPlot(svz_select, group.by = "hash.ID", pt.size = .6, cols = color_pal.8)
ggsave("plots/umap.select.sct.pdf", width = 5.9, height = 4.6)

svz_select.markers <- FindAllMarkers(object=svz_select)

saveRDS(svz_select.markers, paste0("data/svz_select.markers_", Sys.Date(), ".rds"))
saveRDS(svz_select, paste0("data/svz_select_", Sys.Date(), ".rds"))

# More analysis in seurat__batch2_part2.R
# Fin.


sessionInfo()
#   R version 3.6.0 (2019-04-26)
#  [1] viridis_0.5.1               viridisLite_0.3.0          
#  [3] ggthemes_4.2.0              scales_1.0.0               
#  [5] sctransform_0.2.0           Matrix_1.2-17              
#  [7] forcats_0.4.0               stringr_1.4.0              
#  [9] purrr_0.3.2                 readr_1.3.1                
# [11] tidyr_0.8.3                 tibble_2.1.3               
# [13] ggplot2_3.2.0               tidyverse_1.2.1            
# [15] MAST_1.10.0                 SingleCellExperiment_1.6.0 
# [17] SummarizedExperiment_1.14.0 DelayedArray_0.10.0        
# [19] BiocParallel_1.18.0         matrixStats_0.54.0         
# [21] Biobase_2.44.0              GenomicRanges_1.36.0       
# [23] GenomeInfoDb_1.20.0         IRanges_2.18.1             
# [25] S4Vectors_0.22.0            BiocGenerics_0.30.0        
# [27] dplyr_0.8.1                 Seurat_3.1.0.9007  
