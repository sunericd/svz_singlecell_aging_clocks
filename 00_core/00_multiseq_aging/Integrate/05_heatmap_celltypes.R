library(tidyverse)
library(Seurat)
library(ggplot2)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Integrate")

svz <- readRDS("data/multi_intergrated_seurat_Dec2020.rds")

# find markers
#svz.markers <- FindAllMarkers(object=svz)
#saveRDS(svz.markers, paste0("data/svz.markers_", Sys.Date(), ".rds"))
svz.markers <- readRDS("data/svz.markers_2022-03-18.rds")

# save top 10 marker genes as a CSV table
top10 <- svz.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, 'data/top10_celltype_markers.csv')

# use top5 markers per cell type
top5 <- svz.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# make heatmap
p <- DoHeatmap(svz, features=top5$gene, label=F, group.by="Celltype.LowRes")
ggsave("plots/heatmap_celltype.pdf", p, width=7.07, height=6.78)


# OPC vs Oligodendro
VlnPlot(svz, features="Pdgfra", idents=c("Oligodendro", "OPC"))
