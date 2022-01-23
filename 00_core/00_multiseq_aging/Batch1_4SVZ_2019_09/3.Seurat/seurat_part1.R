# LMO + 4 SVZ Seurat Processing
# Batch 1

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


# All subsequent paths will be relative.
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Batch1_4SVZ_2019_09/3.Seurat")
dir.create("plots")
dir.create("data")
#======================================================================================
# Create Seurat object and basic processing
#======================================================================================

# Read in Cellranger count matrix and create basic seurat object
svz.data <- Read10X("../1.CellRanger/filtered_feature_bc_matrix") # 31053 21556
lmo.umi.data <- Read10X("../2.LMO/umi_count/", gene.column=1) # 9 21509

dim(lmo.umi.data)
lmo.umi.data[1:9,1:10]

# Remove negative controls and unmapped row:
lmo.data <- lmo.umi.data[1:4,]

# Remove cells without barcodes and barcodes without cells (~400)
columns_lmo <- paste0(colnames(lmo.data), "-1")
colnames(lmo.data) <- columns_lmo
barcode.intersect <- intersect(colnames(svz.data), columns_lmo) # 21509
svz.data <- svz.data[, barcode.intersect]
lmo.umi.data <- lmo.umi.data[, barcode.intersect]

# Make Seurat Object
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
ElbowPlot(svz, ndims = 50)
ggsave("plots/elbow.pdf", width=7.3, heigh = 4.3)

svz <- FindNeighbors(svz, dims = 1:21)
svz <- FindClusters(svz, resolution = 0.15)
svz <- RunUMAP(svz, dims = 1:21)

#======================================================================================
# Add LMO sample label data to object
#======================================================================================
barcode.intersect <- intersect(colnames(svz), columns_lmo)
lmo.data <- lmo.data[ , barcode.intersect]
svz[["LMO"]] <- CreateAssayObject(counts = lmo.data)
svz <- NormalizeData(object = svz, assay = "LMO", normalization.method = "CLR")
svz <- HTODemux(svz, assay = "LMO", positive.quantile = 0.99)

table(svz$LMO_classification.global) # with more cell filtering early on so 14069 cell input.
 # Doublet Negative  Singlet 
 #    1900     5619     6550 


color_pal.4 <- tableau_color_pal(palette = "Tableau 10")(4)

# Group cells based on the max HTO signal
Idents(svz) <- "LMO_maxID"
RidgePlot(svz, assay = "LMO", features = rownames(svz[["LMO"]])[1:4], ncol = 2, cols=color_pal.4)

# Heatmap version of above ridgeplots (simple but less quant)
HTOHeatmap(svz, assay = "LMO", ncells = 3000) + scale_fill_viridis(option="viridis", direction=1)


# Violin plots of expression by single/doublet/negative catagories
Idents(svz) <- "LMO_classification.global"
VlnPlot(svz, features = "nCount_RNA", pt.size = 0.1, log = F)
# Unlike other datasets, doublets have less expression. Probably some doublets are are background.

#======================================================================================
# Reductions based on LMO counts
# Remove negative cells from object copy
Idents(svz) <- "LMO_classification.global"
svz.subset <- subset(svz, idents = "Negative", invert = TRUE)
color_pal.2 <- rev(tableau_color_pal(palette = "Tableau 10")(2))

# Calculate a distance matrix using HTO
lmo.dist.mtx <- as.matrix(dist(t(GetAssayData(object = svz.subset, assay = "LMO"))))

# Calculate tSNE embeddings with a distance matrix
svz.subset <- RunTSNE(svz.subset, distance.matrix = lmo.dist.mtx, perplexity = 50)

# tSNE based on LMO reduction, colored by doublet/singlet status
DimPlot(svz.subset, cols = color_pal.2, reduction = "tsne")
ggsave("plots/tsne.lmo.doub.sing.pdf", width = 5.9, height = 4.6)

#======================================================================================
# Subset to singlets
color_pal.1 <- tableau_color_pal(palette = "Tableau 10")(1)
color_pal.4 <- tableau_color_pal(palette = "Tableau 10")(4)

# Remove doublets and Negatives
svz_sing <- subset(svz, subset = LMO_classification.global == "Singlet")

# Calculate a distance matrix using HTO
lmo.dist.mtx <- as.matrix(dist(t(GetAssayData(object = svz_sing, assay = "LMO"))))

# Calculate tSNE embeddings with a distance matrix
svz_sing <- RunTSNE(svz_sing, distance.matrix = lmo.dist.mtx, perplexity = 50)

# tSNE based on LMO reduction, colored by sample
DimPlot(svz_sing, cols = color_pal.4, reduction = "tsne", group.by="hash.ID")
ggsave("plots/tsne.lmo.sing.age.pdf", width = 6.9, height = 4.6)


#======================================================================================
# Cell Selector
# Remove about 46 ambigous cells based on dimensionality reduction plot
# Unfortunately not reproducible but very minor change
#======================================================================================

#plot <- DimPlot(svz_sing, cols = color_pal.8, reduction = "tsne", group.by="hash.ID")
#cells.located <- CellSelector(plot = plot)
#ambiguous <- cells.located # 46 ambiguous cells to remove.
#saveRDS(ambiguous, "data/ambiguous_cells.rds")
ambiguous <- readRDS("data/ambiguous_cells.rds")
#  [1] "AACAACCCAACCCTCT-1" "AAGTCGTCACCAATTG-1" "AATGAAGCAGGTTCCG-1"
#  [4] "AATGGAACAGTGACCC-1" "ACACCAACACGACGCT-1" "ACTTATCTCTTACACT-1"
#  [7] "AGCCAATCAACCACAT-1" "AGGATCTCATGGCTGC-1" "AGGTCTAAGTATCTGC-1"
# [10] "ATGCATGCAAAGCTCT-1" "ATGTCCCTCGGTAGAG-1" "CAACCTCCAGGTGGAT-1"
# [13] "CAGATTGCAGGCACTC-1" "CAGCAATCACACGGAA-1" "CAGGCCAGTGCAACGA-1"
# [16] "CATCCACTCTGTGCTC-1" "CATCCCATCCATGCAA-1" "CCGATGGGTACGGTTT-1"
# [19] "CCGTAGGTCGGAATTC-1" "CTCCACAAGGAGACCT-1" "CTGAATGTCCTTTAGT-1"
# [22] "CTGAGCGCAAGGTACG-1" "CTGCATCCATATACCG-1" "CTGCATCTCCATTTCA-1"
# [25] "CTGCCTACAGCATACT-1" "CTGGTCTTCCGAACGC-1" "GACACGCCACCAGCGT-1"
# [28] "GAGCTGCCAGAGGAAA-1" "GAGTTTGAGCAACAGC-1" "GATCAGTTCGAAGAAT-1"
# [31] "GCATTAGGTTCTCTAT-1" "GGCTTTCGTGGATCAG-1" "GGTAATCTCGCAGATT-1"
# [34] "GGTCTGGCACTGCGAC-1" "GTCTGTCAGACTCCGC-1" "GTTCATTGTGTCATCA-1"
# [37] "TAAGCGTTCGGTAAGG-1" "TCAATCTGTCTACAAC-1" "TCACTCGAGCTCGTGC-1"
# [40] "TCATGGATCGTTATCT-1" "TCATTCAGTTACCCTC-1" "TCTATACCAACCGGAA-1"
# [43] "TGCACGGCATAAGCAA-1" "TTACCATCAATTGGTC-1" "TTGGGCGAGTTATGGA-1"
# [46] "TTTGGAGAGGTTCTTG-1"

svz_select <- subset(svz_sing, cells = ambiguous, invert = TRUE) # 6504 good cells
DimPlot(svz_select, cols = color_pal.4, reduction = "tsne", group.by="hash.ID")
ggsave("plots/tsne.lmo.filtered.pdf", width = 5.9, height = 4.6)

#======================================================================================
# Group cells based on the max HTO signal
Idents(svz_select) <- "LMO_maxID"

svz_select <- NormalizeData(object = svz_select, assay = "LMO", normalization.method = "RC")
RidgePlot(svz_select, assay = "LMO", features = rownames(svz[["LMO"]])[1:4], ncol = 2, cols=color_pal.4)

# Heatmap version of above ridgeplot
svz_select <- NormalizeData(object = svz_select, assay = "LMO", normalization.method = "CLR")
Idents(svz_select) <- "LMO_maxID"
HTOHeatmap(svz_select, assay = "LMO", ncells = 5000) + scale_fill_viridis(option="viridis", direction=1)
ggsave("plots/heatmap.lmo.png", width = 6.9, height=2.88)

table(svz_sing$LMO_maxID)
# Sample1-TGTGATGG Sample2-TCAATGGC Sample3-CTCTAGAC Sample4-ACCAATGC 
#             1201             1829             2047             1473 

#======================================================================================
# Run normalization and clustering on fully filtered subset
#======================================================================================

svz_select <- SCTransform(svz_select)

# Reduce
svz_select <- RunPCA(svz_select, verbose = FALSE)
svz_select <- FindNeighbors(svz_select, dims = 1:20)
svz_select <- FindClusters(svz_select, resolution = 0.25)
svz_select <- RunUMAP(svz_select, dims = 1:20)

DimPlot(svz_select, group.by = "hash.ID", pt.size = .6, cols = color_pal.4)
ggsave("plots/umap.select.sct.sample.pdf", width = 7.5, height = 6)

svz_select.markers <- FindAllMarkers(object=svz_select)

saveRDS(svz_select.markers, paste0("data/svz_select.markers_", Sys.Date(), ".rds"))
saveRDS(svz_select, paste0("data/svz_select_", Sys.Date(), ".rds"))

# More analysis in seurat_part2.R
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
