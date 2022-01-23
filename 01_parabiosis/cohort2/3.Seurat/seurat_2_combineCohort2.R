# Combine Demultiplexed Parabiosis Runs into one Seurat Object

library(Seurat)
library(dplyr)
library(MAST)
library(tidyverse)
library(Matrix)
library(sctransform)
library(scales)
library(ggthemes)
library(viridis)

#======================================================================================

setwd("~/Dropbox/svz_singlecell_aging_clocks/01_parabiosis/cohort2/3.Seurat")
runs <- c("GEX1", "GEX2", "GEX4")

g1 <- readRDS("GEX1/data/GEX1svz_sing_2020-08-07.rds")
g2 <- readRDS("GEX2/data/GEX2svz_sing_2020-08-08.rds")
g4 <- readRDS("GEX4/data/GEX4svz_sing_2020-08-07.rds")
obj_list <- c(g1, g2, g4)
obj <- merge(x = obj_list[[1]], y = obj_list[c(2:length(obj_list))])

obj # 48563 features across 13864 samples within 3 assays 
#======================================================================================
# Renormalize together
obj <- SCTransform(obj)

#======================================================================================
# Reduce and Plot
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- RunUMAP(obj,
               dims = 1:20,
               min.dist = 0.5,
               spread = 0.5,
               seed.use = 42)
obj <- FindClusters(obj, resolution = 0.2, assay = "SCT")
DimPlot(obj) + scale_color_tableau(palette = "Tableau 20")
ggsave(paste0("plots/umap.pdf"), width = 6.7, height = 6.5)

#======================================================================================
# Adjust metadata
hash.ID <- obj@meta.data$hash.ID
# Corrects barcode BC14 label (now corrected at source TAG_LIST_GEX2.csv)
# Harmless if run 
hash.ID[hash.ID == "BC14-Young-Iso-5.13-5B-B2-50100-AACCGAAC"] <- "BC14-Young-Het-5.13-5B-B2-50100-AACCGAAC"
obj$hash.ID <- hash.ID

sample_info <- as.data.frame(obj@meta.data$hash.ID)
colnames(sample_info) <- "hashID"
sample_df <- separate(sample_info, hashID, sep = "-",
            into = c("Barcode", "Age", "Type", "Months",
                     "Pair_ID", "Batch", "FACS", "LMO_Barcode"))
sample_df$AgeCond <- paste0(sample_df$Age, "-", sample_df$Type)

# Factor adjustments
sample_df <- tbl_df(sample_df)
sample_df$Batch <- factor(sample_df$Batch, levels = c("B1", "B2", "B4"))
sample_df$Age <- factor(sample_df$Age, levels = c("Young", "Old"))
sample_df$AgeCond <- factor(sample_df$AgeCond,
                        levels = c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"))
sample_df$Months <- as.numeric(sample_df$Months)

sample_df <- data.frame(sample_df)
rownames(sample_df) <- colnames(obj)
obj <- AddMetaData(obj, sample_df)
head(obj[[]])

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

#======================================================================================
# Marker Genes
markers <- FindAllMarkers(obj)

#======================================================================================
# Annotate Cells 
new.cluster.ids <- c("Microglia",
                     "Oligo_1",
                     "Astrocyte_qNSC",
                     "Oligo_2",
                     "Endothelial",
                     "Neuroblast", 
                     "aNSC_NPC",
                     "Mural",
                     "Oligo_3",
                     "Oligo_4",
                     "OPC", # 10
                     "Macrophage",
                     "Neutrophil",
                     "Neuron", 
                     "Doublet")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
obj[["Celltype"]] <- Idents(obj)
unique(obj@meta.data$Celltype)


# Assign lower resolution celltypes
low.res.ids <- c("Microglia",
                 "Oligo",
                 "Astrocyte_qNSC",
                 "Oligo",
                 "Endothelial",
                 "Neuroblast", 
                 "aNSC_NPC",
                 "Mural",
                 "Oligo",
                 "Oligo",
                 "OPC",
                 "Macrophage",
                 "Neutrophil",
                 "Neuron", 
                 "Doublet")

names(low.res.ids) <- levels(obj)
obj <- RenameIdents(obj, low.res.ids)
obj[["Celltype.LowRes"]] <- Idents(obj)
unique(Idents(obj))
obj <- subset(obj, Celltype != "Doublet")

# Don't overwrite
#saveRDS(obj, paste0("data/pb-svz_annotatedSeuratObject_", "2020-08-10", ".rds"))
