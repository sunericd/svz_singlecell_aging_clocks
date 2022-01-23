# Integrate 4 Multi-seq batches

library(Seurat)
library(dplyr)
library(MAST)
library(tidyverse)
library(Matrix)
library(sctransform)
library(scales)
library(ggthemes)
library(viridis)
library(harmony)
library(cowplot)

# All subsequent paths will be relative.
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Integrate")

#======================================================================================
# Load svz seurat objects
#======================================================================================

r1 <- readRDS("../Batch1_4SVZ_2019_09/3.Seurat/data/svz_r1_complete.rds") # 4245 cells
r2 <- readRDS("../Batch2_8SVZ_2019_12/3.Seurat/data/svz_r2_complete.rds") # 4218 cells
r3 <- readRDS("../Batch3_8SVZ_2020_03/3.Seurat/data/svz_r3_complete.rds") # 9416 cells
r4 <- readRDS("../Batch4_8SVZ_2020_11/3.Seurat/data/svz_r4_complete.rds") # 3626

#======================================================================================
# Merge
#======================================================================================
r1$orig.ident <- "Batch-1"
r2$orig.ident <- "Batch-2"
r3$orig.ident <- "Batch-3"
r4$orig.ident <- "Batch-4"
r1[["SCT"]] <- NULL
r2[["SCT"]] <- NULL
r3[["SCT"]] <- NULL

svz <- merge(x = r1, y = c(r2, r3, r4)) # 21505 cells

#======================================================================================
# Harmonize
#======================================================================================
svz <- SCTransform(svz)
svz <- RunPCA(svz, verbose = FALSE)
svz <- RunHarmony(svz, "orig.ident", assay.use="SCT")

svz <- FindNeighbors(svz, reduction = "harmony", dims = 1:20)
svz <- RunUMAP(svz, reduction = "harmony", dims = 1:20, reduction.name = "umap_har", seed.use = 4)
svz <- FindClusters(svz, resolution = 0.30)
DimPlot(svz, reduction = "umap_har", shuffle = TRUE, group.by = "Phase")
ggsave(paste0("plots/umap_harmony_", Sys.Date(), "_phase.pdf"), useDingbats = F)

#======================================================================================
# Joint clustering
#======================================================================================
svz.markers <- FindAllMarkers(object=svz)
svz.markers <- readRDS("data/markers.rds")
#saveRDS(svz.markers, "data/markers.rds")

top10 <- svz.markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_logFC)
top <- rev(unique(top10$gene))
DotPlot(svz, features = top) + coord_flip() + scale_color_viridis()
ggsave("plots/dotplot_markergenes.pdf")

DimPlot(svz, reduction = "umap_har") + scale_color_tableau(palette = "Tableau 20")


#======================================================================================
# Assign celltypes
#======================================================================================
new.cluster.ids <- c("Microglia",
                    "Oligodendro_1",
                    "Neuroblast_2",
                    "Astrocyte_qNSC",
                    "aNSC_NPC_2",
                    "Neuroblast_1",
                    "Oligodendro_3",
                    "Oligodendro_2",
                    "Endothelial",
                    "aNSC_NPC_1",
                    "Mural",
                    "Macrophage",
                    "OPC",
                    "Neuron",
                    "Ependymal",
                    "Doublet")
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
                    "Neuroblast",
                    "Oligodendro",
                    "Oligodendro",
                    "Endothelial",
                    "aNSC_NPC",
                    "Mural",
                    "Macrophage",
                    "OPC",
                    "Neuron",
                    "Ependymal",
                    "Doublet")
names(new.cluster.ids) <- levels(svz)
svz <- RenameIdents(svz, new.cluster.ids)
svz[["Celltype.LowRes"]] <- Idents(svz)
unique(Idents(svz))
DimPlot(svz, reduction = "umap_har") + scale_color_tableau(palette = "Tableau 20")
#======================================================================================
# Proliferative Fraction
#======================================================================================

e1 <- svz[[]] %>% select(Celltype.LowRes, Phase, hash.ID)
e1$Lineage <- as.character(e1$Celltype.LowRes)
e1$Lineage[e1$Lineage %in% c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast") &
            e1$Phase %in% c("G2M", "S")] <- "Prolif_Lineage"
e2 <- e1 %>% group_by(hash.ID, Lineage) %>%
             summarise(n = n()) %>%
             mutate(freq = n / sum(n))
e3 <- dplyr::filter(e2, Lineage == "Prolif_Lineage") %>%
        select(hash.ID, freq)
colnames(e3)[2] <- "Prolif_Lineage_Fraction_of_SVZ"
meta <- svz[[]] %>% rownames_to_column("Cell")
e4 <- plyr::join(meta, e3)
rownames(e4) <- e4$Cell
e4 <- e4 %>% select(Prolif_Lineage_Fraction_of_SVZ)
svz <- AddMetaData(svz, e4)
e4 <- e4 %>% select(Prolif_Lineage_Fraction_of_SVZ, everything())

meta2 <- select(meta, hash.ID, Age, Prolif_Lineage_Fraction_of_SVZ) %>% unique()
ggplot(meta2, aes(x=Age, y=Prolif_Lineage_Fraction_of_SVZ*100)) +
     geom_point() +
     geom_smooth(method = "lm") +
     theme_classic() +
     annotate(x=10, y=10, 
         label=paste("R = ", round(cor(meta2$Age, meta2$Prolif_Lineage_Fraction_of_SVZ),2)), 
         geom="text", size=5)


#======================================================================================
# Visualize with ggplot
#======================================================================================
# There's some slight stochasticity in the projection orientations.
# Read this file to make the same plots:
d <- readRDS("data/metadata_Dec2020.rds")

set.seed(42)
d <- d[sample(nrow(d)), ]
d$Age <- factor(d$Age, levels = c("3.3","3.33","3.6","4.3","4.7","5.4","6.7","8.4",  "9.47","10.43","12.47","14.5","14.77","16.53", "16.83", "18.58", "18.87", "20.6", "20.8", "21.57", "22.57", "22.6", "23.9", "24.9","25.93","29"))
manual_pal <- c("#3d405b", # dark bluish
                "#aa4465", # Purple
                "#f8961e", # spectrum orange, Neurob
                "#277da1", # Spectrum Bight blue, Astro
                "#90be6d", # Spectrum yellow green, NSC
                "#f3722c", # spectrum bright red
                "#f9c74f", # Spectrum yellow
                "#ff87ab", # Pink, OPC
                "#f94144", # Spectrum bright red
                "#9db4c0", # blue-grey
                "#606c38") # green

p1 <- ggplot(d, aes(x = umap_har_1, y = umap_har_2, color = Celltype.LowRes)) +
          geom_point(size = .35, stroke=0) +
          scale_color_manual(values = manual_pal) +
          theme_classic() +
          guides(colour = guide_legend(override.aes = list(size=3)))
p1
p2 <- ggplot(d, aes(x = umap_har_1, y = umap_har_2, color = Age)) +
          geom_point(size = .55, stroke=0) +
          scale_color_viridis(option="viridis", discrete = T, begin = 0, end = 1) +
          theme_classic() +
          guides(colour = guide_legend(override.aes = list(size=3)))
p2
pc12 <- plot_grid(p1, p2, ncol = 2, align = "hv")
ggsave("plots/Figure1bc_UMAPs.png", pc12, dpi = 600, width = 11, height = 4)

p3 <- ggplot(d, aes(x = umap_har_1, y = umap_har_2, color = Phase)) +
          geom_point(size = .55, stroke=0) +
          scale_color_tableau(palette = "Superfishel Stone") +
          theme_classic() +
          guides(colour = guide_legend(override.aes = list(size=3)))
ggsave("plots/Figure1c.png", p3,
     dpi = 600, width = 4.6, height = 4)
#======================================================================================
# Save
#======================================================================================

# Don't overwrite
#saveRDS(svz, file = "data/multi_intergrated_seurat_Dec2020.rds")

meta <- svz[[]] # Metadata. Age, Celltype, and Celltype.LowRes most important
harmony <- Reductions(svz, slot="umap_har")@cell.embeddings # 2D batch effect correction reduction
colnames(harmony) <- c("umap_har_1", "umap_har_2") 
d <- tbl_df(cbind(harmony, meta))
# Don't overwrite
#saveRDS(d, file = "data/metadata_Dec2020.rds")












