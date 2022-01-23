# Batch 1 LMO + 4 SVZ Seurat Processing
# PART TWO

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
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Batch1_4SVZ_2019_09/3.Seurat")
dir.create("plots2")

#======================================================================================
# Load data
markers <- readRDS("data/svz_select.markers_2020-04-12.rds")
svz <- readRDS("data/svz_select_2020-04-12.rds")
t20 <- tableau_color_pal(palette = "Tableau 20")(20)
t4 <- tableau_color_pal(palette = "Tableau 10")(4)[c(2,3,4,1)]

svz <- RunUMAP(svz, dims = 1:20, min.dist = .3, spread = .9, seed.use = 42)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(svz, features = top10$gene) + NoLegend() # Observe cluster 0 is noise

#======================================================================================
# Cluster based filter
#======================================================================================

# Remove Group 0
svz <- subset(svz, subset = SCT_snn_res.0.25 != 0) # 4639 samples

# Recluster
svz <- FindNeighbors(svz, dims = 1:20)
svz <- FindClusters(svz, resolution = 0.6)
svz <- RunUMAP(svz, dims = 1:20, min.dist = .3, spread = .9, seed.use = 41)
DimPlot(svz, cols = t20)

# Find markers again now that noise cluster is gone
markers2 <- FindAllMarkers(svz)
top10 <- markers2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(svz, features = top10$gene, group.colors = t20) + NoLegend() # Clean

#======================================================================================
# Redo LMO based plots with new filter
#======================================================================================

# Group cells based on the max HTO signal
svz <- NormalizeData(object = svz, assay = "LMO", normalization.method = "CLR")
svz <- HTODemux(svz, assay = "LMO", positive.quantile = 0.99)
svz <- NormalizeData(object = svz, assay = "LMO", normalization.method = "CLR")
HTOHeatmap(svz, assay = "LMO", ncells = 3000) + scale_fill_viridis(option="viridis", direction=1)
ggsave("plots2/heatmap.lmo.all.png", width = 5.9, height=2.88)

svz_sing <- subset(svz, subset = LMO_classification.global == "Singlet")
table(svz_sing[[c("LMO_classification.global")]])
HTOHeatmap(svz_sing, assay = "LMO", ncells = 3000) + scale_fill_viridis(option="viridis", direction=1)
ggsave("plots2/heatmap.lmo.sing.png", width = 5.9, height=2.88)

svz <- svz_sing
#======================================================================================
# Assign celltypes
#======================================================================================
new.cluster.ids <- c("Astrocyte",
                    "Microglia_1",
                    "Oligodendro_1",
                    "Oligodendro_2",
                    "Oligodendro_3",
                    "Oligodendro_4",
                    "Neuroblast",
                    "aNSC_NPC",
                    "Endothelial",
                    "Oligodendro_5",
                    "Microglia_2",
                    "Mural",
                    "Macrophage",
                    "Ependymal")
Idents(svz) <- "SCT_snn_res.0.6"
names(new.cluster.ids) <- levels(svz)
svz <- RenameIdents(svz, new.cluster.ids)
svz[["Celltype"]] <- Idents(svz)

# Assign lower resolution celltypes
low.res.ids <- c("Astrocyte_qNSC",
                "Microglia",
                "Oligodendro",
                "Oligodendro",
                "Oligodendro",
                "Oligodendro",
                "Neuroblast",
                "aNSC_NPC",
                "Endothelial",
                "Oligodendro",
                "Microglia",
                "Mural",
                "Macrophage",
                "Ependymal")
names(low.res.ids) <- levels(svz)
svz <- RenameIdents(svz, low.res.ids)
svz[["Celltype.LowRes"]] <- Idents(svz)

DimPlot(svz, reduction = "umap", label = TRUE, label.size = 4, pt.size = .5, cols = t20)
ggsave("plots2/umap.celltype.pdf", width=9, height=6)


#======================================================================================
# Add more metadata
#======================================================================================

head(svz[[]])

# Add Age
svz$Age <- as.character(svz$LMO_maxID)
svz$Age[svz$Age == "Sample1-TGTGATGG"] <- "4.7"
svz$Age[svz$Age == "Sample2-TCAATGGC"] <- "6.7"
svz$Age[svz$Age == "Sample3-CTCTAGAC"] <- "20.8"
svz$Age[svz$Age == "Sample4-ACCAATGC"] <- "29.0"
svz$Age <- as.numeric(svz$Age)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
svz <- CellCycleScoring(svz, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

DefaultAssay(svz) <- "RNA"


#======================================================================================
# Counts
#======================================================================================
table(svz[[c("Celltype.LowRes", "Age")]])
# Celltype.LowRes 4.7 6.7 20.8  29
#     Astrocyte   258 344  144  81
#     Microglia    48 158  224 281
#     Oligodendro 337 491  543 461
#     Neuroblast  116 111   45  23
#     aNSC_NPC     89  95   37  33
#     Endothelial  17  35   44  46

# Remove unwanted metadata columns
m.cols <- colnames(svz[[]])
good.cols <- !grepl("*_snn_*", m.cols)
svz@meta.data <- svz@meta.data[,good.cols] # includes 21 metadata columns.

# Don't overwrite
saveRDS(svz, file = "data/svz_r1_complete.rds")


#FIN
