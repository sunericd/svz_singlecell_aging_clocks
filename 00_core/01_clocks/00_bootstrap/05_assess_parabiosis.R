library(ggthemes)
library(cowplot)
library(tidyverse)
library(scales)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")
d <- readRDS("data/parabiosis_predictions.rds")

d$Batch <- d$Mouse
d$Batch[!grepl("2019", d$Batch)] <- "2"
d$Batch[grepl("2019", d$Batch)] <- "1"
d$Batch <- factor(d$Batch, levels = c("1", "2"))
d$Sample <- factor(d$Sample, levels = c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"))
d$Celltype <- factor(d$Celltype, levels = c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
                                            "Microglia", "Oligodendro", "Endothelial"))
df <- d
#=============================================================================
MouseMedians <- df %>%
    group_by(Celltype, Sample, Batch, Mouse) %>%
    summarise(MedianMouseAges = median(Pred))

SampleMedians <- MouseMedians %>% 
    group_by(Celltype, Sample, Batch) %>%
    summarise(MedianPseudoCellAge = median(MedianMouseAges))

SampleMeans <- MouseMedians %>% 
    group_by(Celltype, Sample, Batch) %>%
    summarise(MedianPseudoCellAge = mean(MedianMouseAges))

# Sample Medians
effect_table <- c()
for (celltype in unique(SampleMedians$Celltype)) {
    print(celltype)
    for (batch in c("1", "2")) {
        print(batch)
        # Young blood effect:
        d <- filter(SampleMedians, Batch == batch, Celltype == celltype)
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
et$Celltype <- factor(et$Celltype, levels = c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
                                            "Microglia", "Oligodendro", "Endothelial"))
data.frame(et)

#          Celltype Batch  Condition      Effect
# 1  Astrocyte_qNSC     1 YoungBlood  3.78794604
# 2  Astrocyte_qNSC     1   OldBlood -1.89734904
# 3  Astrocyte_qNSC     2 YoungBlood -0.64647650
# 4  Astrocyte_qNSC     2   OldBlood -1.99972394
# 5        aNSC_NPC     1 YoungBlood  5.37724286 <<<
# 6        aNSC_NPC     1   OldBlood -0.98442196
# 7        aNSC_NPC     2 YoungBlood  4.18418116 <<<
# 8        aNSC_NPC     2   OldBlood -0.53492940
# 9      Neuroblast     1 YoungBlood -1.22330041
# 10     Neuroblast     1   OldBlood -1.60123542
# 11     Neuroblast     2 YoungBlood  0.89523543
# 12     Neuroblast     2   OldBlood -0.04422853
# 13      Microglia     1 YoungBlood  3.37415032
# 14      Microglia     1   OldBlood -1.53136914
# 15      Microglia     2 YoungBlood  1.48453302
# 16      Microglia     2   OldBlood  0.38165705
# 17    Oligodendro     1 YoungBlood  2.01489730
# 18    Oligodendro     1   OldBlood -0.35142001
# 19    Oligodendro     2 YoungBlood -1.22156522
# 20    Oligodendro     2   OldBlood -1.41409357
# 21    Endothelial     1 YoungBlood  2.08316634
# 22    Endothelial     1   OldBlood -1.26547066
# 23    Endothelial     2 YoungBlood -1.29646614
# 24    Endothelial     2   OldBlood -0.49157901


# Sample Means Version
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
et$Celltype <- factor(et$Celltype, levels = c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
                                            "Microglia", "Oligodendro", "Endothelial"))
data.frame(et)

#          Celltype Batch  Condition      Effect
# 1  Astrocyte_qNSC     1 YoungBlood  3.78794604
# 2  Astrocyte_qNSC     1   OldBlood -1.89734904
# 3  Astrocyte_qNSC     2 YoungBlood -0.44299253
# 4  Astrocyte_qNSC     2   OldBlood -1.51717851
# 5        aNSC_NPC     1 YoungBlood  5.37724286 <<<
# 6        aNSC_NPC     1   OldBlood -0.98442196
# 7        aNSC_NPC     2 YoungBlood  3.65949844 <<<
# 8        aNSC_NPC     2   OldBlood -1.53022649
# 9      Neuroblast     1 YoungBlood -1.22330041
# 10     Neuroblast     1   OldBlood -1.60123542
# 11     Neuroblast     2 YoungBlood  0.92876414
# 12     Neuroblast     2   OldBlood  0.40590629
# 13      Microglia     1 YoungBlood  3.37415032
# 14      Microglia     1   OldBlood -1.53136914
# 15      Microglia     2 YoungBlood -0.33292093
# 16      Microglia     2   OldBlood -0.08398407
# 17    Oligodendro     1 YoungBlood  2.01489730
# 18    Oligodendro     1   OldBlood -0.35142001
# 19    Oligodendro     2 YoungBlood -0.68202183
# 20    Oligodendro     2   OldBlood -1.15986406
# 21    Endothelial     1 YoungBlood  2.08316634
# 22    Endothelial     1   OldBlood -1.26547066
# 23    Endothelial     2 YoungBlood -1.63741403
# 24    Endothelial     2   OldBlood -0.50128428

# Mean of batches for aNSC
# > (5.37724286 + 3.65949844) / 2
# [1] 4.518371

# Avg acceleration effect
mean(filter(et, Condition == "OldBlood")$Effect) # -1.001491
# Avg reversal effect
mean(filter(et, Condition == "YoungBlood")$Effect) # 1.408918