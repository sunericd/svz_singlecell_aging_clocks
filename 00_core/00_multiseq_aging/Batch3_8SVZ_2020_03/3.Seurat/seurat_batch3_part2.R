# LMO + 8 SVZ Seurat Process
# Batch 3

# PART TWO

# VERSION 3.0 Seurat

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
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Batch3_8SVZ_2020_03/3.Seurat")

#======================================================================================
# Load data
markers <- readRDS("data/svz_sing.markers.2020-04-08.rds")
svz <- readRDS(paste0("data/svz_sing_", "2020-04-08", ".rds"))
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
                    "Neuroblast_2",
                    "Oligodendro_1",
                    "Neuroblast_1",
                    "aNSC_NPC_2",
                    "Astrocyte_qNSC",
                    "Oligodendro_2",
                    "aNSC_NPC_1",
                    "Endothelial",
                    "Mural",
                    "Oligodendro_3",
                    "OPC",
                    "Macrophage")
names(new.cluster.ids) <- levels(svz)
svz <- RenameIdents(svz, new.cluster.ids)
svz[["Celltype"]] <- Idents(svz)
unique(svz@meta.data$Celltype)

# Assign lower resolution celltypes
low.res.ids <- c("Microglia",
                "Neuroblast",
                "Oligodendro",
                "Neuroblast",
                "aNSC_NPC",
                "Astrocyte_qNSC",
                "Oligodendro",
                "aNSC_NPC",
                "Endothelial",
                "Mural",
                "Oligodendro",
                "OPC",
                "Macrophage")
names(low.res.ids) <- levels(svz)
svz <- RenameIdents(svz, low.res.ids)
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
    x <- gsub("BC[1-8]-", "", input)
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
# Celltype.LowRes  3.3 4.3 8.4 10.43 12.47 21.57 22.6
#   Microglia      270 336 175   276   319   151  126
#   Neuroblast     628 540 327   439   517   175  102
#   Oligodendro    437 327 293   526   369   198  219
#   aNSC_NPC       351 319 166   187   215    77   41
#   Astrocyte_qNSC  92 150  92   142   158    48   64
#   Endothelial     41  29  27    48    57    24   25
#   Mural           40  11  22    26    39     3   17
#   OPC             18  14  11    13    12     2   13
#   Macrophage       9  11  10    15    14     6    7

# Remove unwanted metadata columns
m.cols <- colnames(svz[[]])
good.cols <- !grepl("*_snn_*", m.cols)
svz@meta.data <- svz@meta.data[,good.cols]

#Don't overwrite
#saveRDS(svz, file = "data/svz_r3_complete.rds") # includes 21 metadata columns.

#FIN

