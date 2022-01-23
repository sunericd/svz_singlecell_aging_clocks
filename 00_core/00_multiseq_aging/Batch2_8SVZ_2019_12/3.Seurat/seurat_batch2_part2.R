# Batch 2 LMO + 4 SVZ Seurat Processing
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
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Batch2_8SVZ_2019_12/3.Seurat")


#======================================================================================
# Load data
markers <- readRDS("data/svz_select.markers_2020-04-12.rds")
svz <- readRDS("data/svz_select_2020-04-12.rds")
t20 <- tableau_color_pal(palette = "Tableau 20")(20)
t10 <- tableau_color_pal(palette = "Tableau 10")(10)

svz <- RunUMAP(svz, dims = 1:20, min.dist = .5, spread = 1.2, seed.use = 42)
DimPlot(svz, cols = t20)
ggsave(paste0("plots2/umaps/unfilted_", 42, ".pdf"))

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(svz, features = top10$gene, lines.width = 10) + NoLegend()

#======================================================================================
# Assign celltypes
#======================================================================================
new.cluster.ids <- c("Microglia_1",
                    "Oligodendro_1",
                    "Astrocyte",
                    "Oligodendro_2",
                    "Neuroblast_2",
                    "Endothelial",
                    "aNSC_NPC",
                    "qNSC_Primed",
                    "Neuroblast_1",
                    "Mural",
                    "Macrophage",
                    "Oligodendro_3",
                    "OPC",
                    "Ependymal",
                    "Microglia_2")

names(new.cluster.ids) <- levels(svz)
svz <- RenameIdents(svz, new.cluster.ids)
svz[["Celltype"]] <- Idents(svz)
unique(svz@meta.data$Celltype)

# Assign lower resolution celltypes
low.res.ids <- c("Microglia",
                "Oligodendro",
                "Astrocyte_qNSC",
                "Oligodendro",
                "Neuroblast",
                "Endothelial",
                "aNSC_NPC",
                "Astrocyte_qNSC",
                "Neuroblast",
                "Mural",
                "Macrophage",
                "Oligodendro",
                "OPC",
                "Ependymal",
                "Microglia")

names(low.res.ids) <- levels(svz)
svz <- RenameIdents(svz, low.res.ids)
svz[["Celltype.LowRes"]] <- Idents(svz)

DimPlot(svz, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1, cols = t20)
ggsave("plots2/umap.celltype.wide.pdf", width=11, height=6)


#======================================================================================
# Add more metadata
#======================================================================================

head(svz[[]])

# Add Age
extractAge <- function(input) {
    x <- gsub("^BC[0-9]+.-", "", input)
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
#-----------------------------------------------------------
# Celltype.LowRes  3.33 5.4 9.47 14.5 16.53 18.58 20.6 22.57
#   Microglia        52 210  230  156    24    39    9    17
#   Oligodendro      98 367  146  259    31   105   28    53
#   Astrocyte_qNSC  115 212   57  130    28    77   49    68
#   Neuroblast       84 170   81  211    19    20    9    15
#   Endothelial      38 108  149   59    10    20    7     6
#   aNSC_NPC         73 106   69   66     7    14    3     7
#   ...
#-----------------------------------------------------------

# Remove unwanted metadata columns
m.cols <- colnames(svz[[]])
good.cols <- !grepl("*_snn_*", m.cols)
svz@meta.data <- svz@meta.data[,good.cols]

# Don't overwrite
#saveRDS(svz, file = "data/svz_r2_complete.rds") # includes 21 metadata columns.


#FIN

