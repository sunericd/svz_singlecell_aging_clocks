library(ggthemes)
library(cowplot)
library(ggthemes)
library(tidyverse)
library(scales)


setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/01_bootstrap_bioage")

d <- readRDS("data/parabiosis_bootstrapCell_predictions_bio.rds")
cell_order <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

d$Batch <- d$Mouse
d$Batch[!grepl("2019", d$Batch)] <- "2"
d$Batch[grepl("2019", d$Batch)] <- "1"
d$Batch <- factor(d$Batch, levels = c("1", "2"))
d$Sample <- factor(d$Sample, levels = c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"))
d$Celltype <- factor(d$Celltype, levels = cell_order)
df <- d

df$Pred <- 35 - 100*df$Pred

#=============================================================================
MouseMedians <- df %>%
    group_by(Celltype, Sample, Batch, Mouse) %>%
    summarise(MedianMouseAges = median(Pred))

SampleMeans <- MouseMedians %>% 
    group_by(Celltype, Sample, Batch) %>%
    summarise(MedianPseudoCellAge = mean(MedianMouseAges))

effect_table <- c()
for (celltype in unique(SampleMeans$Celltype)) {
    print(celltype)
    for (batch in c("1", "2")) {
        print(batch)
        # Young blood effect:
        d <- filter(SampleMeans, Batch == batch, Celltype == celltype)
        eff <- d[d$Sample == "Old-Iso", "MedianPseudoCellAge"] - d[d$Sample == "Old-Het", "MedianPseudoCellAge"]
        effect_table <- rbind(effect_table, c(celltype, batch, "YoungBlood", eff))
        # Old Blood effects
        eff2 <- d[d$Sample == "Young-Iso", "MedianPseudoCellAge"] - d[d$Sample == "Young-Het", "MedianPseudoCellAge"]
        effect_table <- rbind(effect_table, c(celltype,  batch, "OldBlood", eff2))
    }
}
colnames(effect_table) <- c("Celltype", "Batch", "Condition", "Effect")
effect_table

et <- data.frame(effect_table) %>% unnest(cols = c(Celltype, Batch, Condition, Effect))
et$Celltype <- factor(et$Celltype, levels = cell_order)
et$Clock <- "Bootstrap"
et$CondBatch <- paste0(et$Condition, "_", et$Batch)
et$CondBatch <- factor(et$CondBatch, levels = c("YoungBlood_2", "YoungBlood_1",
                                                "OldBlood_2", "OldBlood_1"))

et_bootstrap <- et


#=============================================================================
#=============================================================================
rm(et)
setwd("../02_ensemble")
d <- read.csv("data/parabiosis_predicted_fraction.csv")
d <- as_tibble(d)

d$Batch <- d$batch
d$Batch[!grepl("2019", d$Batch)] <- "2"
d$Batch[grepl("2019", d$Batch)] <- "1"
d$Batch <- factor(d$Batch, levels = c("1", "2"))
d$Sample <- factor(d$condition, levels = c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"))
d$Celltype <- factor(d$celltype, levels = cell_order)
d$Pred <- 35 - 100*d$predicted_fraction
d$Mouse <- d$hash_id
df <- d

#=============================================================================
MouseMedians <- df %>%
    group_by(Celltype, Sample, Batch, Mouse) %>%
    summarise(MedianMouseAges = median(Pred))

SampleMeans <- MouseMedians %>% 
    group_by(Celltype, Sample, Batch) %>%
    summarise(MedianPseudoCellAge = mean(MedianMouseAges))

effect_table <- c()
for (celltype in unique(SampleMeans$Celltype)) {
    print(celltype)
    for (batch in c("1", "2")) {
        print(batch)
        # Young blood effect:
        d <- filter(SampleMeans, Batch == batch, Celltype == celltype)
        eff <- d[d$Sample == "Old-Iso", "MedianPseudoCellAge"] - d[d$Sample == "Old-Het", "MedianPseudoCellAge"]
        effect_table <- rbind(effect_table, c(celltype, batch, "YoungBlood", eff))
        # Old Blood effects
        eff2 <- d[d$Sample == "Young-Iso", "MedianPseudoCellAge"] - d[d$Sample == "Young-Het", "MedianPseudoCellAge"]
        effect_table <- rbind(effect_table, c(celltype,  batch, "OldBlood", eff2))
    }
}
colnames(effect_table) <- c("Celltype", "Batch", "Condition", "Effect")
effect_table

et <- data.frame(effect_table) %>% unnest(cols = c(Celltype, Batch, Condition, Effect))
et$Celltype <- factor(et$Celltype, levels = cell_order)

et$Clock <- "Ensemble"
et$CondBatch <- paste0(et$Condition, "_", et$Batch)
et$CondBatch <- factor(et$CondBatch, levels = c("YoungBlood_2", "YoungBlood_1",
                                                "OldBlood_2", "OldBlood_1"))
data.frame(et)

#          Celltype Batch  Condition      Effect    Clock    CondBatch
# 1     Oligodendro     1 YoungBlood  1.68105848 Ensemble YoungBlood_1
# 2     Oligodendro     1   OldBlood  0.34128147 Ensemble   OldBlood_1
# 3     Oligodendro     2 YoungBlood -2.20907287 Ensemble YoungBlood_2
# 4     Oligodendro     2   OldBlood -0.72875601 Ensemble   OldBlood_2
# 5       Microglia     1 YoungBlood  2.80197206 Ensemble YoungBlood_1
# 6       Microglia     1   OldBlood -0.71631569 Ensemble   OldBlood_1
# 7       Microglia     2 YoungBlood -1.04692554 Ensemble YoungBlood_2
# 8       Microglia     2   OldBlood  0.02380623 Ensemble   OldBlood_2
# 9     Endothelial     1 YoungBlood  0.36236099 Ensemble YoungBlood_1
# 10    Endothelial     1   OldBlood -0.43043169 Ensemble   OldBlood_1
# 11    Endothelial     2 YoungBlood -0.89641507 Ensemble YoungBlood_2
# 12    Endothelial     2   OldBlood -0.17391599 Ensemble   OldBlood_2
# 13 Astrocyte_qNSC     1 YoungBlood  3.33263575 Ensemble YoungBlood_1
# 14 Astrocyte_qNSC     1   OldBlood -1.58901770 Ensemble   OldBlood_1
# 15 Astrocyte_qNSC     2 YoungBlood -1.03453265 Ensemble YoungBlood_2
# 16 Astrocyte_qNSC     2   OldBlood -1.09808814 Ensemble   OldBlood_2
# 17       aNSC_NPC     1 YoungBlood  3.56938675 Ensemble YoungBlood_1
# 18       aNSC_NPC     1   OldBlood -0.83908132 Ensemble   OldBlood_1
# 19       aNSC_NPC     2 YoungBlood  1.44102023 Ensemble YoungBlood_2
# 20       aNSC_NPC     2   OldBlood -0.40330812 Ensemble   OldBlood_2
# 21     Neuroblast     1 YoungBlood  2.58821284 Ensemble YoungBlood_1
# 22     Neuroblast     1   OldBlood -1.05730834 Ensemble   OldBlood_1
# 23     Neuroblast     2 YoungBlood  1.77020788 Ensemble YoungBlood_2
# 24     Neuroblast     2   OldBlood -0.81565049 Ensemble   OldBlood_2


#=============================================================================
#=============================================================================

et2 <- rbind(et, et_bootstrap)


ggplot(data = et2, mapping = aes(x = Clock, y = CondBatch, size = abs(Effect), color = -Effect)) +
        geom_point() +
        scale_color_distiller(palette="RdBu", limits = c(-3,3), oob = scales::squish) +
        geom_point(shape = 1, color = "black") +
        facet_grid(.~Celltype) +
        theme_cowplot() +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7),   strip.text.x = element_text(size = 6)) +
        theme(legend.position = "bottom") +
        theme(axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1.15, hjust=1)) +
        scale_size_area(
            max_size = 6,
            breaks = c(0,1,2,3,5),
            labels = c("0","1","2","3","5+"),
            guide = "legend",
            limits = c(0, 5),
            oob = scales::squish)

ggsave("plots/dotplot_combo_parabiosis_bio.pdf", useDingbats=F, width = 4, height = 2.35)
