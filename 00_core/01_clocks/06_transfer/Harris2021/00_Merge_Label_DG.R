library(tidyverse)
library(Seurat)
library(harmony)

setwd("~/Dropbox/MULTIseq/Integrated")
svz <- readRDS("data/multi_intergrated_seurat_Dec2020.rds")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/06_transfer/Harris2021")

input_folder <- "data"
seurat_object <- Read10X(data.dir = input_folder)
seurat_object <- CreateSeuratObject(counts = seurat_object, project = "Aggregate")

# Add GEMgroup ----
seurat_object@meta.data$GEMgroup[grepl("-1", rownames(seurat_object@meta.data))] <- "1"
seurat_object@meta.data$GEMgroup[grepl("-2", rownames(seurat_object@meta.data))] <- "2"
seurat_object@meta.data$GEMgroup[grepl("-3", rownames(seurat_object@meta.data))] <- "3"
seurat_object@meta.data$GEMgroup[grepl("-4", rownames(seurat_object@meta.data))] <- "4"
seurat_object@meta.data$GEMgroup[grepl("-5", rownames(seurat_object@meta.data))] <- "5"
seurat_object@meta.data$GEMgroup[grepl("-6", rownames(seurat_object@meta.data))] <- "6"
seurat_object@meta.data$GEMgroup[grepl("-7", rownames(seurat_object@meta.data))] <- "7"
seurat_object@meta.data$GEMgroup[grepl("-8", rownames(seurat_object@meta.data))] <- "8"

# Add Batch no - got this from 10X website ----
seurat_object@meta.data$Batch[grepl("-1", rownames(seurat_object@meta.data))] <- "1"
seurat_object@meta.data$Batch[grepl("-2|-3", rownames(seurat_object@meta.data))] <- "2"
seurat_object@meta.data$Batch[grepl("-4", rownames(seurat_object@meta.data))] <- "3"
seurat_object@meta.data$Batch[grepl("-5", rownames(seurat_object@meta.data))] <- "4"
seurat_object@meta.data$Batch[grepl("-6", rownames(seurat_object@meta.data))] <- "5"
seurat_object@meta.data$Batch[grepl("-7", rownames(seurat_object@meta.data))] <- "6"
seurat_object@meta.data$Batch[grepl("-8", rownames(seurat_object@meta.data))] <- "7"

# Add Age of samples ----
seurat_object@meta.data$Age[grepl("-1", rownames(seurat_object@meta.data))] <- "1mo"
seurat_object@meta.data$Age[grepl("-2|-3|-4|-5", rownames(seurat_object@meta.data))] <- "2mo"
seurat_object@meta.data$Age[grepl("-6|-7|-8", rownames(seurat_object@meta.data))] <- "6mo"

# GFP and tdtomato or combined ----
seurat_object@meta.data$cells[grepl("-1|-4|-5|-6|-7|-8", rownames(seurat_object@meta.data))] <- "combined"
seurat_object@meta.data$cells[grepl("-2", rownames(seurat_object@meta.data))] <- "tdTomato"
seurat_object@meta.data$cells[grepl("-3", rownames(seurat_object@meta.data))] <- "Gfp"

# Add 10x chemistry ----
seurat_object@meta.data$Chemistry[grepl("-1|-6|-7|-8", rownames(seurat_object@meta.data))] <- "V3"
seurat_object@meta.data$Chemistry[grepl("-2|-3|-4|-5", rownames(seurat_object@meta.data))] <- "V2"

# Add MT content ----
seurat_object[["percentMito"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")

#==================================================================

# Removedoublets with oligos ----
seurat_object <- subset(seurat_object, subset = Mog > 4 & Aldoc > 4, slot = 'counts', invert = TRUE) %>%
        subset(subset = Mog > 4 & C1qc > 4, slot = 'counts', invert = TRUE) %>%
        subset(subset = Mog > 4 & Cldn5 > 4, slot = 'counts', invert = TRUE) %>%
        subset(subset = Mog > 4 & Stmn2 > 4, slot = 'counts', invert = TRUE) %>%
        # Removedoublets with astro
        subset(subset = Aldoc > 4 & C1qc > 4, slot = 'counts', invert = TRUE) %>%
        subset(subset = Aldoc > 4 & Cldn5 > 4, slot = 'counts', invert = TRUE) %>%
        subset(subset = Aldoc > 4 & Stmn2 > 4, slot = 'counts', invert = TRUE) %>%
        # Removedoublets with microglia
        subset(subset = C1qc > 4 & Cldn5 > 4, slot = 'counts', invert = TRUE)  %>%
        subset(subset = C1qc > 4 & Stmn2 > 4, slot = 'counts', invert = TRUE)  %>%
        subset(subset = Cldn5 > 4 & Stmn2 > 4, slot = 'counts', invert = TRUE)

# Set thresholds for gene and mt content ----
mitoHi <- 10
nGeneLo <- 500
seurat_object <- subset(seurat_object, subset = nFeature_RNA > nGeneLo & percentMito < mitoHi)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

head(seurat_object[[]], 20)
print(colnames(seurat_object[[]]))
seurat_object@meta.data$region <- "DG"
svz@meta.data$region <- "SVZ"
svz@meta.data$GEMgroup <- "9"

svzdg <- merge(svz, seurat_object)
DefaultAssay(svzdg) <- "RNA"
svzdg[["LMO"]] <- NULL
svzdg[["SCT"]] <- NULL

#svzdg <- SCTransform(svzdg)
svzdg <- NormalizeData(svzdg)
svzdg <- FindVariableFeatures(svzdg)
svzdg <- ScaleData(svzdg)
svzdg <- RunPCA(svzdg)

svzdg <- RunHarmony(svzdg, "GEMgroup", assay.use="RNA")

svzdg <- FindNeighbors(svzdg, reduction = "harmony", dims = 1:20)
svzdg <- RunUMAP(svzdg, reduction = "harmony", dims = 1:20, reduction.name = "umap_har", seed.use = 4)

DimPlot(svzdg, reduction = "umap_har", group.by = "GEMgroup")
DimPlot(svzdg, reduction = "umap_har", group.by = "Celltype.LowRes")
ggsave("plots/svz_dg_harmony_umap.pdf")

svzdg <- FindClusters(svzdg, resolution = 0.10)
DimPlot(svzdg, reduction = "umap_har", label=T)


# Create list of all numbered cluster names, replace old names for the clusters of NSCs and IPCs, check labelling ———
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12", "13", "14")
new.cluster.ids <- c("Oligodendro","Neuroblast","2","Astrocyte_qNSC","aNSC_NPC","Microglia","Endothelial","7","8","9","10","11","12","13", "14")
svzdg@active.ident <- plyr::mapvalues(x = svzdg@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(svzdg, reduction = "umap_har", label=T)
dg <- subset(svzdg, region == "DG")
dg@meta.data$Celltype.LowRes <- dg@active.ident
dg <- subset(dg, Celltype.LowRes %in% c("Oligodendro","Neuroblast", "Astrocyte_qNSC","aNSC_NPC","Microglia","Endothelial"))
dg[["SCT"]] <- NULL
DimPlot(dg, reduction = "umap_har", label=T)

saveRDS(dg, "data/dg_object.rds")










