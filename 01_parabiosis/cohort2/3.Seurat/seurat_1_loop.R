# For loop LMO Demultiplexing GEX1 Parabiosis Cohort 2
# 4 lanes, but lane 3 labelling failed so only 1, 2, and 4 are used in paper.
# SVZ Seurat Processing
# Matthew Buckley
# Started: August 7, 2020

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
setwd("~/Dropbox/svz_singlecell_aging_clocks/01_parabiosis/cohort2/3.Seurat")
dir.create("plots")
dir.create("data")

#======================================================================================
# Create Seurat object and basic processing
#======================================================================================
runs <- c("GEX2","GEX4")
# Test
#run <- "GEX1"
# run <- "GEX4"
# run <- "GEX2"

for (run in runs) {
    print(run)
    dir.create(paste0(run,"/plots"), recursive = TRUE)
    dir.create(paste0(run,"/data"), recursive = TRUE)

    # Read in Cellranger count matrix and create basic seurat object
    svz.data <- Read10X(paste0("../1.CellRanger/", run, "/filtered_feature_bc_matrix"))
    cell_bcs <- colnames(svz.data)
    cell_bcs_trim <- substr(cell_bcs,1,nchar(cell_bcs)-2)
    write.csv(cell_bcs_trim, paste0("WHITE_LIST_", run, ".csv"), row.names = F)
    lmo.umi.data <- Read10X(paste0("../2.LMO/", run, "/umi_count"), gene.column = 1) # 10 12413
    dim(lmo.umi.data)

    # Remove negative control barcode and unmapped row:
    lmo.data <- lmo.umi.data[1:(nrow(lmo.umi.data) - 1), ]

    # Remove cells without barcodes and barcodes without cells
    barcode.intersect <- intersect(cell_bcs_trim, colnames(lmo.data))
    barcode.intersect.untrim <-  paste0(barcode.intersect, "-1")
    svz.data <- svz.data[, barcode.intersect.untrim]
    lmo.umi.data <- lmo.umi.data[, barcode.intersect]
    colnames(lmo.umi.data) <- paste0(barcode.intersect, "-1")

    svz <- CreateSeuratObject(counts = svz.data, project = run)
    svz <- NormalizeData(svz)
    svz <- FindVariableFeatures(svz)
    svz <- ScaleData(svz)


    #======================================================================================
    # Filter
    #======================================================================================

    # cut off of 1000 UMI instead of elbow from cell ranger
    svz[["percent.mt"]] <- PercentageFeatureSet(svz, pattern = "^mt-")
    svz <- subset(svz, subset = nCount_RNA > 1000 & percent.mt < 10 & nFeature_RNA > 500)

    # Visualize QC metrics as a violin plot
    VlnPlot(svz, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,  pt.size = .1)
    ggsave(paste0(run, "/plots/violin.mt.feature.counts.fil.pdf"), width = 10, height = 4)

    #======================================================================================
    # Dim reduction and clustering
    #======================================================================================

    # Run PCA, select 21 PCs for tSNE visualization and graph-based clustering
    svz <- RunPCA(svz, verbose = FALSE)
    svz <- FindNeighbors(svz, dims = 1:21)
    svz <- RunUMAP(svz, dims = 1:21)
    svz <- FindClusters(svz, resolution = 0.5)

    #======================================================================================
    # Add LMO sample label data to object
    #======================================================================================
    cell_bcs <- colnames(svz)
    cell_bcs_trim <- substr(cell_bcs,1,nchar(cell_bcs)-2)

    barcode.intersect <- intersect(cell_bcs_trim, colnames(lmo.data))
    lmo.data <- lmo.data[ , barcode.intersect]

    barcode.intersect.untrim <-  paste0(barcode.intersect, "-1")
    colnames(lmo.data) <- barcode.intersect.untrim

    svz[["LMO"]] <- CreateAssayObject(counts = lmo.data)
    svz <- NormalizeData(object = svz, assay = "LMO", normalization.method = "CLR")

    svz <- HTODemux(svz, assay = "LMO", positive.quantile = 0.99, nsamples = 2000)
    print(table(svz$LMO_classification.global))

    N <- length(unique(svz@meta.data$LMO_maxID))

    color_pal.8 <- tableau_color_pal(palette = "Tableau 10")(N)

    # Group cells based on the max HTO signal
    Idents(svz) <- "LMO_maxID"

    # Heatmap
    HTOHeatmap(svz, assay = "LMO", ncells = 3000) + scale_fill_viridis(option="viridis", direction=1)
    ggsave(paste0(run, "/plots/heatmap.lmo.all.png"), width = 10.9, height=2.88)


    #======================================================================================
    # Reductions based on LMO counts

    # Set three colors
    color_pal.3 <- tableau_color_pal(palette = "Tableau 10")(3)[c(2,3,1)]

    # Calculate a distance matrix using HTO
    lmo.dist.mtx <- as.matrix(dist(t(GetAssayData(object = svz, assay = "LMO"))))

    # Calculate tSNE embeddings with distance matrix
    svz <- RunTSNE(svz, distance.matrix = lmo.dist.mtx, perplexity = 30)

    # tSNE based on LMO reduction, colored by global lmo classifiction
    Idents(svz) <- "LMO_classification.global"
    DimPlot(svz, cols = color_pal.3, reduction = "tsne")
    ggsave(paste0(run, "/plots/tsne.lmo.doub.sing.neg.pdf"), width = 5.9, height = 4.6)

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

    # tSNE based on LMO reduction, colored id
    DimPlot(svz_sing, cols = color_pal.8, reduction = "tsne", group.by="hash.ID")
    ggsave(paste0(run, "/plots/tsne.lmo.sing.age.pdf"), width = 6.9, height = 4.6, useDingbats=FALSE)

    #======================================================================================
    # Group cells based on the max HTO signal
    Idents(svz_sing) <- "LMO_maxID"
    # Heatmap version of above ridgep1 (simple but less quant)
    svz_sing <- NormalizeData(object = svz_sing, assay = "LMO", normalization.method = "CLR")
    Idents(svz_sing) <- "LMO_maxID"
    HTOHeatmap(svz_sing, assay = "LMO", ncells = 2000) + scale_fill_viridis(option="viridis", direction=1)
    ggsave(paste0(run, "/plots/heatmap.lmo.select.png"), width = 6.9, height=2.88)


    #======================================================================================
    # Rerun proces on subset
    #======================================================================================

    # Normalize, Standardize
    svz_sing <- SCTransform(svz_sing)

    # Reduce
    svz_sing <- RunPCA(svz_sing, verbose = FALSE)
    svz_sing <- FindNeighbors(svz_sing, dims = 1:20)
    svz_sing <- RunUMAP(svz_sing, dims = 1:20)
    svz_sing <- FindClusters(svz_sing, resolution = 0.25)

    DimPlot(svz_sing, group.by = "hash.ID", pt.size = .5, cols = color_pal.8)
    ggsave(paste0(run, "/plots/umap.select.sct.sample.pdf"), width = 12, height = 7)

    # Don't overwrite
    #saveRDS(svz_sing, paste0(run,"/data/", run, "svz_sing_", Sys.Date(), ".rds"))

}

# Fin.
