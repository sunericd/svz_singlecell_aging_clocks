# Make Seurat object from Allen Brain Atlas MTG dataset

library(tidyverse)
library(Matrix)
library(Seurat)
library(plyr)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/06_transfer/HumanMTG")


raw_counts<-read.table(file="raw_data/human_MTG_2018-06-14_exon-matrix.csv",sep=",",header=T,row.names=1)
meta<-read.table(file="raw_data/human_MTG_2018-06-14_samples-columns.csv",sep=",",header=T,row.names=1)

gene.meta<-read.table(file="raw_data/human_MTG_2018-06-14_genes-rows.csv",sep=",",header=T)
rownames(raw_counts) <- mapvalues(rownames(raw_counts), from=gene.meta$entrez_id, to=gene.meta$gene)



raw_counts <- CreateSeuratObject(raw_counts, meta.data=meta)

saveRDS(raw_counts, "human_mtg_seurat.rds")