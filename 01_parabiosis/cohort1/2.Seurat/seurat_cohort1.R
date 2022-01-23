# Parabiosis 1 Seurat Processing
# Matthew Buckley
# UPDATED: 2/12/2019

library(Seurat)
library(dplyr)
library(tidyverse)
library(Matrix)

setwd("~/Dropbox/svz_singlecell_aging_clocks/01_parabiosis/cohort1/1.CellRanger")
# Read in 10X samples
het_o <- Read10X(data.dir = "HET-O/filtered_gene_bc_matrices")
het_y <- Read10X(data.dir = "HET-Y/filtered_gene_bc_matrices")
yy <- Read10X(data.dir = "YY1/filtered_gene_bc_matrices")
oo <- Read10X(data.dir = "OO1/filtered_gene_bc_matrices")

setwd("../2.Seurat")
# Create one data object
#### Generated with an older version of Seurat
svz <- CreateSeuratObject(raw.data = yy, project = "yy", min.cells = 0, min.genes = 0)
# Add samples
svz <- AddSamples(object = svz, new.data = oo, add.cell.id = "oo")
svz <- AddSamples(object = svz, new.data = het_y, add.cell.id = "het_y")
svz <- AddSamples(object = svz, new.data = het_o, add.cell.id = "het_o")

# Inspect
svz; head(svz@meta.data); tail(svz@meta.data)

# Mitochondria Genes QC
mito.genes <- grep(pattern = "^mt-", x = rownames(x = svz@data), value = TRUE)
percent.mito <- Matrix::colSums(svz@raw.data[mito.genes, ])/Matrix::colSums(svz@raw.data)
svz <- AddMetaData(object = svz, metadata = percent.mito, col.name = "percent.mito")

# Filter. Added Hba-a1 to remove contaminating red blood cells
svz <- FilterCells(object = svz, subset.names = c("nGene", "percent.mito", "Hba-a1"), low.thresholds = c(600, -Inf, -Inf), high.thresholds = c(4500, 0.15, 3))

# Normalize
svz <- NormalizeData(object = svz, normalization.method = "LogNormalize", scale.factor = 10000)

# Find Variable Genes
svz <- FindVariableGenes(object = svz, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.75)
length(svz@var.genes) # 7000 genes

# Regress out
svz <- ScaleData(object = svz, vars.to.regress = c("nUMI", "nGene"))

# Cell cycle assignments
svz <- CellCycleScoring(object = svz, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)

# PCA
svz <- RunPCA(object = svz, pc.genes = svz@var.genes, pcs.compute = 100)

# Final cluster & tSNE choice.
svz <- FindClusters(object = svz, reduction.type = "pca", dims.use = 1:20, resolution = .2, print.output = 0, force.recalc = TRUE)

svz <- RunTSNE(object = svz, dims.use = 1:20, perplexity = 30, seed.use = 5)

svz.markers <- FindAllMarkers(object = svz, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25, type = "MAST")

top10 <- svz.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top50 <- svz.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
new.cluster.ids <- c("Oligodendrocytes", "Endothelial","Microglia","Astrocytes_qNSCs",
					 "aNSCs_NPCs", "Oligodendrocytes", "Neuroblasts", "Mural_cells",
					 "Macrophages_NK", "Neurons", "OPC", "Endothelial")
svz@ident <- plyr::mapvalues(x = svz@ident, from = current.cluster.ids, to = new.cluster.ids)
svz <- AddMetaData(object = svz, metadata = svz@ident, col.name = "Celltype")

# Don't overwrite
#save(svz, file = paste0("data/svz_postAssigned", Sys.Date(), ".rda"))