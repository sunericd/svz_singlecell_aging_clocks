library(ggthemes)
library(cowplot)
library(tidyverse)
library(scales)
library(ggridges)
library(ggpubr)

CELLTYPES <- c("Oligodendro", "Microglia",  "Endothelial",
               "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

#=====================================================================
#=====================================================================
## BootstrapCell Predictions- Parabiosis

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")
d <- readRDS("data/parabiosis_predictions.rds")
d$Batch <- d$Mouse
d$Batch[!grepl("2019", d$Batch)] <- "2"
d$Batch[grepl("2019", d$Batch)] <- "1"
d$Batch <- factor(d$Batch, levels = c("1", "2"))
d$Sample <- plyr::mapvalues(x = d$Sample, from = c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"), to = c("Young Isochronic", "Young Heterochronic", "Old Isochronic", "Old Heterochronic"))
d$Sample <- factor(d$Sample, levels = c("Young Isochronic", "Young Heterochronic", "Old Isochronic", "Old Heterochronic"))
d$Celltype <- factor(d$Celltype, levels = CELLTYPES)


# get median for each mouse
df_med <- d %>% filter(Batch=="2") %>% group_by(Mouse, Celltype, Sample) %>% summarise(median_Pred=median(Pred)) %>% ungroup()
#df_med <- df_med %>% filter(Sample %in% c("Old Heterochronic", "Old Isochronic"))



##############################
# Parabiosis
##############################

p <- ggplot(df_med, aes(x=Sample, y=median_Pred, fill=Sample)) +
  geom_violin() +
  geom_point(color="black") +
  ylab("Median Chronological Age Prediction") + xlab("Sample") +
  stat_compare_means(method="wilcox.test", comparisons=list(c("Young Isochronic", "Young Heterochronic"), c("Old Isochronic", "Old Heterochronic"))) +
  facet_wrap(~ Celltype) +
  scale_fill_tableau() +
  scale_color_tableau() +
  theme_classic()
p
ggsave("plots/parabiosis2_oldhet_oldiso_scatter.pdf", p, width=10 , height=6)







d <- readRDS("data/exercise_predictions.rds")

d$Batch <- as.character(d$Year)
d$Batch <- factor(d$Batch, levels = c("R1", "R2"))
d$Sample <- plyr::mapvalues(x = d$Sample, from = c("Y_Control", "Y_Exercise", "O_Control", "O_Exercise"), to = c("Young Sedentary", "Young Exercise", "Old Sedentary", "Old Exercise"))
d$Sample <- factor(d$Sample, levels = c("Young Sedentary", "Young Exercise", "Old Sedentary", "Old Exercise"))
d$Celltype <- factor(d$Celltype, levels = CELLTYPES)

# get median for each mouse
df_med <- d %>% group_by(Mouse, Celltype, Sample) %>% summarise(median_Pred=median(Pred)) %>% ungroup()
#df_med <- df_med %>% filter(Sample %in% c("Old Exercise", "Old Sedentary"))


##############################
# Exercise
##############################

p <- ggplot(df_med, aes(x=Sample, y=median_Pred, fill=Sample)) +
  geom_violin() +
  geom_point(color="black") +
  ylab("Median Chronological Age Prediction") + xlab("Sample") +
  stat_compare_means(method="wilcox.test", comparisons=list(c("Young Exercise", "Young Sedentary"), c("Old Exercise", "Old Sedentary"))) +
  facet_wrap(~ Celltype) +
  scale_fill_tableau() +
  scale_color_tableau() +
  theme_classic()
p
ggsave("plots/exercise_oldhet_oldiso_scatter.pdf", p, width=10 , height=6)

