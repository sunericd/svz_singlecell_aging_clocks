library(Seurat)
library(tidyverse)
library(ggthemes)
library(cowplot)


setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Integrate")

Celltypes <- c("Oligodendro", "Microglia", "Endothelial",
               "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

#======================================== Shared Genes

svz <- readRDS("data/multi_intergrated_seurat_Dec2020.rds")
meta <- svz[[]] # Metadata. Age, Celltype, and Celltype.LowRes most important
harmony <- Reductions(svz, slot="umap_har")@cell.embeddings # 2D batch effect correction reduction
colnames(harmony) <- c("umap_har_1", "umap_har_2") 
umap <- Reductions(svz, slot="umap")@cell.embeddings # Non corrected UMAP reduction
d <- tbl_df(cbind(umap, harmony, meta))

DefaultAssay(svz) <- "RNA"
exps <- svz[["RNA"]]@data


#========================================
# Bootstrap Cell Chrono Celltype Specific Specific Genes
celltypeDF <- readRDS("../../01_clocks/05_visualize/data/models_upset_data.rds")
generalG <- rownames(celltypeDF[rowSums(celltypeDF) >  1, ])


genes1 <- c("Map2", "Asap1")#aNSC
genes2 <- c("Rida", "Tspan7")#qnsc
genes3 <- c("Pkm","Itga1") # Endo
genes4 <- c("Crlf2", "Fcrls") # Micro
genes5 <- c("Kansl2", "Ap1s1") # Neuro
genes6 <- c("Gstm1", "Cd59a") # Oligo

genes <- c(genes1, genes2, genes3, genes4, genes5, genes6)


exps2 <- exps[rownames(exps) %in% genes, ]
exprs_t <- t(as.matrix(exps2))
meta_counts <- cbind(d, exprs_t)
mclong <- meta_counts %>% pivot_longer(c(36:length(colnames(meta_counts))), names_to = "gene", values_to = "expression")
mclong$Celltype.LowRes <- factor(mclong$Celltype.LowRes, levels = Celltypes)

mclong2 <- mclong %>%
            select(hash.ID, Age, orig.ident, Celltype.LowRes, gene, expression) %>%
            group_by(hash.ID, Age, orig.ident, Celltype.LowRes, gene) %>%
            filter(Celltype.LowRes %in% Celltypes) %>%
            summarize(Mean = mean(expression))


mclong3 <- mclong2 %>%
            filter((gene %in% genes1 & Celltype.LowRes == "aNSC_NPC") |
                    (gene %in% genes2 & Celltype.LowRes == "Astrocyte_qNSC") |
                    (gene %in% genes3 & Celltype.LowRes == "Endothelial") |
                    (gene %in% genes4 & Celltype.LowRes == "Microglia") |
                    (gene %in% genes5 & Celltype.LowRes == "Neuroblast") |
                    (gene %in% genes6 & Celltype.LowRes == "Oligodendro")
                    )

# Fig 2d
ggplot(mclong3, aes(x = Age, y = Mean, color = gene)) +
    facet_wrap(Celltype.LowRes~., scales= "free_y") +
    geom_point() +
    geom_smooth(method = "loess", aes(fill = gene)) +
    theme_cowplot() +
    scale_color_tableau(palette = "Tableau 20") +
    scale_fill_tableau(palette = "Tableau 20")

ggsave("plots/specificgene_trajectories.pdf", width=7.14, height=3.53)



#===================================================
# Shared across several bootstrap chrono clocks

exps2 <- exps[rownames(exps) %in% c("Ifi27", "AC149090.1", "Tma7"), ]
exprs_t <- t(as.matrix(exps2))
meta_counts <- cbind(d, exprs_t)
mclong <- meta_counts %>% pivot_longer(c(36:length(colnames(meta_counts))), names_to = "gene", values_to = "expression")
mclong$Celltype.LowRes <- factor(mclong$Celltype.LowRes, levels = Celltypes)

mclong2 <- mclong %>%
            select(hash.ID, Age, orig.ident, Celltype.LowRes, gene, expression) %>%
            group_by(hash.ID, Age, orig.ident, Celltype.LowRes, gene) %>%
            filter(Celltype.LowRes %in% Celltypes) %>%
            #filter(gene %in% c("Ifi27", "AC149090.1", "Mdk", "Pax6")) %>%
            summarize(Mean = mean(expression))
# Fig 2e
ggplot(mclong2, aes(x = Age, y = Mean, color = gene)) +
    facet_wrap(Celltype.LowRes~., scales= "free_y") +
    geom_point() +
    geom_smooth(method = "loess", aes(fill = gene)) +
    theme_cowplot() +
    scale_color_few() +
    scale_fill_few()
    
ggsave("plots/threegene_trajectories.pdf", width=7.14, height=3.53)




