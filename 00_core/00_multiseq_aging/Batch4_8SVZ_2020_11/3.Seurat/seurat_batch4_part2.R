# LMO + 8 SVZ Seurat Process

# PART TWO

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

#======================================================================================
# Load data
markers <- readRDS("data/svz_sing.markers.2020-12-12.rds")
svz <- readRDS(paste0("data/svz_sing_", "2020-12-12", ".rds"))
t20 <- tableau_color_pal(palette = "Tableau 20")(20)
t10 <- tableau_color_pal(palette = "Tableau 10")(10)

svz <- RunUMAP(svz, dims = 1:20, min.dist = .5, spread = 1.2, seed.use = 41)
DimPlot(svz, cols = t20)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(svz, features = top10$gene, lines.width = 10) + NoLegend()

#======================================================================================
# Assign celltypes
#======================================================================================
new.cluster.ids <- c("Microglia",
                    "Oligodendro_1",
                    "Neuroblast",
                    "Astrocyte_qNSC",
                    "aNSC_NPC_1",
                    "aNSC_NPC_2",
                    "Oligodendro_2",
                    "Endothelial",
                    "OPC",
                    "Mural",
                    "Macrophage")
names(new.cluster.ids) <- levels(svz)
svz <- RenameIdents(svz, new.cluster.ids)
svz[["Celltype"]] <- Idents(svz)
unique(svz@meta.data$Celltype)

# Assign lower resolution celltypes
new.cluster.ids <- c("Microglia",
                    "Oligodendro",
                    "Neuroblast",
                    "Astrocyte_qNSC",
                    "aNSC_NPC",
                    "aNSC_NPC",
                    "Oligodendro",
                    "Endothelial",
                    "OPC",
                    "Mural",
                    "Macrophage")
names(new.cluster.ids) <- levels(svz)
svz <- RenameIdents(svz, new.cluster.ids)
svz[["Celltype.LowRes"]] <- Idents(svz)
unique(Idents(svz))

DimPlot(svz, reduction = "umap", label = TRUE, label.size = 5, pt.size = .2, cols = t10)
ggsave("plots2/umap.celltype.wide.pdf", width=9, height=6)

#======================================================================================
# Add more metadata
#======================================================================================

head(svz[[]])

# Add Age
extractAge <- function(input) {
    x <- gsub("BC[0-9]*-", "", input)
    gsub("-[ATCG]*$", "", x)
}
Ages <- as.character(svz@meta.data$hash.ID)
Ages <- unlist(lapply(Ages, extractAge))
svz$Age <- as.numeric(Ages)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
svz <- CellCycleScoring(svz, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

DefaultAssay(svz) <- "RNA"


#======================================================================================
# Counts
#======================================================================================
table(svz[[c("Celltype.LowRes", "Age")]])
#                 Age
# Celltype.LowRes  3.6 14.77 16.83 18.87 23.9 24.9 25.93
#   Microglia       91    90   141   153  107  165   126
#   Oligodendro    126    82   132   166   75  170   148
#   Neuroblast     165    89    97   108   39   88    83
#   Astrocyte_qNSC  38    64    68    73   35   79    50
#   aNSC_NPC       167    70   106    87   25   28    53
#   Endothelial     12     8    28    28    9   19    10
# ...

# Remove unwanted metadata columns
m.cols <- colnames(svz[[]])
good.cols <- !grepl("*_snn_*", m.cols)
svz@meta.data <- svz@meta.data[,good.cols] # 3626   21

# Don't overwrite
#saveRDS(svz, file = "data/svz_r4_complete.rds") # includes 21 metadata columns.


#FIN

