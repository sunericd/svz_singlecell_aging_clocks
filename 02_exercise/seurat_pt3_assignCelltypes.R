# Seurat Process Exercise Data
# Experiment 2 

library(tidyverse)
library(cerebroApp)
library(Seurat)

#setwd("/labs/abrunet1/Buckley/10.Exercise/2.Seurat")
setwd("~/Dropbox/Exercise/Run2/2.Seurat")

# args <- commandArgs()
# tissue <- args[6]
tissue <- "SVZ"
print(tissue)
obj <- readRDS(paste0("data/data/seurat.", tissue, ".2020-04-21.rds"))

#==================================================================================
# Add metadata
#==================================================================================
ids <- obj@meta.data$orig.ident
meta <- separate(obj[[]], col=orig.ident,
        into = c("ID", "Tissue", "Rep", "Day", "Condition", "Weight"),
        sep = "_", remove = F, convert = T)
meta$Age <- substr(meta$ID, 0, 1)
meta$AgeCond <- paste0(meta$Age, "_", meta$Condition)
meta$Dissection <- d
obj <- AddMetaData(obj, meta)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

#======================================================================================
# Assign celltypes
#======================================================================================
new.cluster.ids <- c("Microglia_1",#0
                    "Astrocyte",#1 
                    "Endothelial",#2
                    "Neuroblast_1",#3
                    "Oligodendro_1",#4
                    "Neuroblast_2",#5
                    "Oligodendro_2",#6
                    "Oligodendro_3",#7
                    "NSPC_G2M-Phase",#8
                    "NSPC_S-Phase",#9
                    "Pericyte",#10
                    "Microglia_2",#11
                    "OPC_1",#12
                    "Primed_qNSC",#13
                    "Smooth_Muscle",#14 (higher alpha smooth muscle actin than in pericytes)
                    "OPC_2",#15
                    "Macrophage",#16,
                    "Oligodendro_4",#17
                    "Doublet_Astrocyte-Oligodendro",#18
                    "Doublet_Microglia-Oligodendro",#19
                    "T_Cell",#20
                    "Doublet_Microglia-Astrocyte",#21
                    "Doublet_Microglia-Neuroblast",#22
                    "Inflammed_qNSC",#23
                    "Doublet_Oligodendro-Neuroblast",#24
                    "Epithelial_Choroid-Plexux",#25,
                    "Doublet_Microglia-Endothelial",#26
                    "Neuronal",#27
                    "Doublet_Endothelial-Neuroblast",#28
                    "Vascular_Leptomeningeal",#29,
                    "Doublet_Pericyte_Microglia",#30
                    "Doublet_Microglia-NSPC",#31
                    "Ependymal",#32
                    "Doublet_Astrocyte-Pericyte"#33
)
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
obj[["Celltype"]] <- Idents(obj)
unique(obj@meta.data$Celltype)

# Assign lower resolution celltypes
low.res.ids <- c("Microglia",#0
                "Astrocyte_qNSC",#1 
                "Endothelial",#2
                "Neuroblast",#3
                "Oligodendro",#4
                "Neuroblast",#5
                "Oligodendro",#6
                "Oligodendro",#7
                "aNSC_NPC",#8
                "aNSC_NPC",#9
                "Mural",#10
                "Microglia",#11
                "OPC",#12
                "Astrocyte_qNSC",#13
                "Mural",#14
                "OPC",#15
                "Macrophage",#16,
                "Oligodendro",#17
                "Doublet",#18
                "Doublet",#19
                "T_Cell",#20
                "Doublet",#21
                "Doublet",#22
                "Astrocyte_qNSC",#23
                "Doublet",#24
                "Epithelial",#25,
                "Doublet",#26
                "Neuronal",#27
                "Doublet",#28
                "Vascular_Leptomeningeal",#29,
                "Doublet",#30
                "Doublet",#31
                "Ependymal",#32
                "Doublet"#33
)
names(low.res.ids) <- levels(obj)
obj <- RenameIdents(obj, low.res.ids)
obj[["Celltype.LowRes"]] <- Idents(obj)
unique(Idents(obj))

DimPlot(obj, reduction = "umap", label = TRUE, label.size = 5, pt.size = .1)

# table(obj[["Celltype.LowRes"]])
#               Microglia          Astrocyte_qNSC             Endothelial 
#                   18267                   13144                    7294 
#              Neuroblast             Oligodendro                aNSC_NPC 
#                   10777                   14961                    7294 
#                   Mural                     OPC              Macrophage 
#                    3566                    2911                     667 
#                 Doublet                  T_Cell              Epithelial 
#                    1550                     226                     132 
#                Neuronal Vascular_Leptomeningeal               Ependymal 
#                     123                      87                      39


#obj <- subset(obj, subset = Celltype.LowRes != "Doublet")
DimPlot(obj, reduction = "umap", label = TRUE, label.size = 5, pt.size = .1)
ggsave(paste0("plots/umap.celltype.", tissue, ".pdf"))

#saveRDS(obj, paste0("data/data/seurat.", tissue, ".annotated.2020-04-27.rds"))

setwd("~/Dropbox/Exercise2020_Shared/SeuratObjects")
saveRDS(obj, paste0("seurat.", tissue, ".annotated.2020-04-27.rds"))





