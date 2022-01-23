library(tidyverse)
library(viridis)
library(RColorBrewer)
library(cowplot)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/05_visualize")
d <- tbl_df(read.csv("data/comparison_data_bioage.csv"))
colnames(d)[1] <- "Celltype"
d$R <- sqrt(d$R2)
d$MAE <- d$MAE * 100

CELLTYPES <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")
d$Celltype <- factor(d$Celltype, levels = CELLTYPES)

d$RNA <- factor(d$RNA, levels = c("Pseudobulk", "BootstrapPseudocell_Median", "EnsemblePseudoCell_Median", "SingleCell"))

#==============================================================================

p1 <- ggplot(data = d, mapping = aes(y = RNA, x = Celltype, size = R, color = MAE)) +
    geom_point() +
    #scale_color_viridis(option = "inferno", begin = 0.1, end = .9, direction = -1) +
    scale_color_distiller(palette="Blues", limits = c(1, 8)) +
    geom_point(shape = 1, color = "black") +
    scale_size(range = c(2,10), limits = c(0, 1), breaks = c(.1,.3,.5,.7,.9)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 30, hjust=1)) +
    theme(legend.position = "bottom")

#ggsave(paste0("plots/dotplot_", Sys.Date(), "_leg_bioage", ".pdf"), p1, useDingbats=F, width=4.02, height=2.79) 

#==============================================================================

d <- tbl_df(read.csv("data/comparison_data.csv"))
colnames(d)[1] <- "Celltype"
d$R <- sqrt(d$R2)
d <- filter(d, RNA != "BootstrapPseudocell_Distrib", RNA != "EnsemblePseudoCell_Distrib")
CELLTYPES <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")
d$Celltype <- factor(d$Celltype, levels = CELLTYPES)
d$RNA <- factor(d$RNA, levels = c("Pseudobulk", "BootstrapPseudocell_Median",
                                  "EnsemblePseudoCell_Median", "SingleCell"))
p2 <- ggplot(data = d, mapping = aes(y = RNA, x = Celltype, size = R, color = MAE)) +
    geom_point() +
    scale_color_distiller(palette="Blues", limits = c(1, 8)) +
    geom_point(shape = 1, color = "black") +
    scale_size(range = c(2,10), limits = c(0, 1), breaks = c(.1,.3,.5,.7,.9)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 30, hjust=1)) +
    theme(legend.position = "bottom")

plot_grid(p2, p1)
ggsave(paste0("plots/dotplot_", Sys.Date(), "_bio_leg_combo", ".pdf"), useDingbats=F, width=7.63, height=2.76)


# For vertical legend only
p2 <- ggplot(data = d, mapping = aes(y = RNA, x = Celltype, size = R, color = MAE)) +
    geom_point() +
    scale_color_distiller(palette="Blues", limits = c(1, 8)) +
    geom_point(shape = 1, color = "black") +
    scale_size(range = c(2,10), limits = c(0, 1), breaks = c(.1,.3,.5,.7,.9)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 30, hjust=1)) +
    theme(legend.position = "right")
ggsave(paste0("plots/dotplot_", Sys.Date(), "_vertical", ".pdf"))
