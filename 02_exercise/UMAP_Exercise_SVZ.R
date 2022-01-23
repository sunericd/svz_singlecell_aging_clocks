library(Seurat)
library(tidyverse)
library(viridis)
library(ggthemes)
library(cowplot)
library(scales)


setwd("~/Dropbox/Exercise/Run2/2.Seurat/")

svz <- readRDS("~/Dropbox/Exercise2020_Shared/SeuratObjects/seurat.SVZ.annotated.2020-04-27.rds")


meta <- svz[[]]
umap <- Reductions(svz, slot="umap")@cell.embeddings
d <- tbl_df(cbind(umap, meta))

d$AgeCond <- d$AgeCond  %>% as_factor() %>%
    fct_recode(YoungControl="Y_Control", YoungExercise="Y_Exercise",
               OldControl="O_Control", OldExercise="O_Exercise")

d$AgeCond <- fct_relevel(d$AgeCond, "YoungControl", "YoungExercise", "OldControl", "OldExercise")
d$Age <- d$Age %>% as_factor() #%>% as.character() %>% as.numeric()
d <- d[sample(nrow(d)),]
d <- filter(d, Celltype.LowRes != "Doublet")

COLORS <- c("#2E3049", # Microglia
			"#1E6A92", #Astro
			"#F3732C", #Endo
			"#FBB040", #Neuroblast
			"#993254", #oligo
			"#85B860", #aNSC
			"#FFDE17", #Mural
			"#F16D9A", #OPC
			"#EE2E39", #Macrophage
			"#C574AF", #Tcell
			"#008E04", #Epithelial
			"#8CA5B1", #Neuron
			"#37F8FF", #Ependymal/Vascular-Lepto
			"black") #Ependymal/Vascular-Lepto

p <- ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = Celltype.LowRes)) +
    geom_point(size = 0.3, stroke = 0, shape = 16) + 
    scale_color_manual(values = COLORS) +
    scale_fill_tableau() +
    theme_classic() +
    guides(color=guide_legend(ncol=1, override.aes = list(size=4)))
p
ggsave("plots/Pretty_SVZ_October_Celltype.png", p, dpi = 600, width = 7, height = 4)