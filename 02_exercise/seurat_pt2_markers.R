# Seurat Process Exercise Data
# Analysis Part 2
# Cluster Marker Genes

library(Seurat)
library(tidyverse)
library(future)
options(future.globals.maxSize= 2000000000)

#setwd("/labs/abrunet1/Buckley/10.Exercise/2.Seurat") # cluster
setwd("~/Dropbox/Exercise/Run2/2.Seurat") # local

# LOCAL TEST
tissue <- "SVZ"

# CLUSTER
# args <- commandArgs()
# tissue <- args[6]

print(tissue)
obj <- readRDS(paste0("data/seurat.", tissue, ".2020-04-20.rds"))
obj.subsample <- subset(obj, downsample = 100)
#plan("multiprocess", workers = 4)
#obj.markers <- FindAllMarkers(object=obj)
#saveRDS(obj.markers, paste0("data/markers.", tissue, ".", Sys.Date(), ".rds"))
obj.markers <- readRDS(paste0("data/markers.", tissue, ".rds"))

top10 <- obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(obj.subsample, features = top10$gene, lines.width = 10) + NoLegend()
ggsave(paste0("plots/", tissue, ".heatmap.markers.big.png"), height = 40, width = 40) 








