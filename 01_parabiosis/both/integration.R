# Purpose: Merge 2019 and 2020 Parabiosis Data

library(Seurat)
library(dplyr)
library(tidyverse)
library(Matrix)
library(sctransform)
library(scales)
library(ggthemes)
library(viridis)
library(harmony)
sessionInfo()

#======================================================================================
# Load 2020 Data object

setwd("~/Dropbox/svz_singlecell_aging_clocks/01_parabiosis/")
pb2020 <- readRDS("cohort2/3.Seurat/data/pb-svz_annotatedSeuratObject_2020-08-10.rds")

#======================================================================================
# Load 2019 Data object

load("cohort1/2.Seurat/data/svz_postAssigned2019-02-14.rda")
pb2019 <- UpdateSeuratObject(svz)
rm(svz)

#======================================================================================
# Remove SCT and Merge

DefaultAssay(pb2020) <- "RNA"
pb2020[['SCT']] <- NULL
pb2020[['LMO']] <- NULL
pb2020[['Celltype']] <- NULL
pb2020[['Celltype']] <- pb2020@meta.data$Celltype.LowRes
pb2020[['Celltype.LowRes']] <- NULL
pb2020[['Experiment']] <- '2020'

pb2019[["AgeCond"]] <- pb2019@meta.data$orig.ident
pb2019[["Experiment"]] <- '2019'
pb <- merge(pb2020, y = pb2019, add.cell.ids = c("20", "19"))
pb
# An object of class Seurat
# 31071 features across 25651 samples within 2 assays
# Active assay: RNA (31053 features)

#======================================================================================
# Non corrected reduction/visualization
rm(list=c("pb2020","pb2019"))
pb <- SCTransform(pb, conserve.memory = TRUE)
pb <- RunPCA(pb, assay = "SCT", verbose = FALSE)

#======================================================================================
# Batch Correction

pb <- RunHarmony(pb, "Experiment", assay.use="SCT")

#======================================================================================

pb <- FindNeighbors(pb, reduction = "harmony", dims = 1:15)

pb <- RunUMAP(pb, reduction = "harmony",
              dims = 1:15,
              reduction.name = "umap_har",
              reduction.key = "Harmony_UMAP",
              min.dist = .4,
              spread = 1.6,
              seed.use = 5)

DimPlot(pb, reduction = "umap_har") + scale_color_tableau(palette = "Tableau 20")

pb <- FindClusters(pb, resolution = .53)
DimPlot(pb, reduction = "umap_har") + scale_color_tableau(palette = "Tableau 20")

markers <- FindAllMarkers(pb)

top <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#======================================================================================
# Annotate Cells 
levels(pb)
pb <- RenameIdents(pb,
                    '0' = "Microglia",
                     '1' = "Astrocyte_qNSC",
                     '2' = "Endothelial",
                     '3' = "Oligodendro",
                     '4' = "Oligodendro",
                     '5' = "Oligodendro",
                     '6' = "Neuroblast", 
                     '7' = "Oligodendro",
                     '8' = "aNSC_NPC",
                     '9' = "Mural",
                     '10' = "aNSC_NPC", #10
                     '11' = "Macrophage_Tcell", # 11
                     '12' = "Oligodendro", # 12
                     '13' = "Astrocyte_qNSC", #13
                     '14' = "OPC",
                     '15' = "Ependymal", #15
                     '16' = "Doublet",
                     '17' = "Doublet",
                     '18' = "Neuron", #18
                     '19' = "Doublet")
levels(pb)
pb[["Celltype2"]] <- Idents(pb)
unique(pb@meta.data$Celltype)
pb <- subset(pb, subset = Celltype != "Doublet")

# Don't overwrite
# saveRDS(pb, "data/pb_combined.rds")
