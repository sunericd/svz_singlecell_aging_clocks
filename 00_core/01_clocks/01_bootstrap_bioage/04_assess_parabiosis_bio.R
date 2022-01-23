library(tidyverse)
library(ggthemes)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/01_bootstrap_bioage")
d <- readRDS("data/parabiosis_bootstrapCell_predictions_bio.rds")


d$Batch <- d$Mouse
d$Batch[!grepl("2019", d$Batch)] <- "2"
d$Batch[grepl("2019", d$Batch)] <- "1"
d$Batch <- factor(d$Batch, levels = c("1", "2"))
d$Sample <- factor(d$Sample, levels = c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"))
d$Celltype <- factor(d$Celltype, levels = c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
                                            "Microglia", "Oligodendro", "Endothelial"))
df <- d

# Convert fraction to score
df$Pred <- 35 - 100*df$Pred

#=============================================================================
MouseMedians <- df %>%
    group_by(Celltype, Sample, Batch, Mouse) %>%
    summarise(MedianMouseAges = median(Pred))

SampleMeans <- MouseMedians %>% 
    group_by(Celltype, Sample, Batch) %>%
    summarise(MedianPseudoCellAge = mean(MedianMouseAges))

#=============================================================================

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
#          Celltype Batch  Condition       Effect
# 1  Astrocyte_qNSC     1 YoungBlood  2.383736919
# 2  Astrocyte_qNSC     1   OldBlood -2.431726594
# 3  Astrocyte_qNSC     2 YoungBlood -1.416349192
# 4  Astrocyte_qNSC     2   OldBlood -1.211008221
# 5        aNSC_NPC     1 YoungBlood  2.579259143
# 6        aNSC_NPC     1   OldBlood -0.841252512
# 7        aNSC_NPC     2 YoungBlood  0.670768494
# 8        aNSC_NPC     2   OldBlood -1.131673968
# 9      Neuroblast     1 YoungBlood  2.165708999 <<<
# 10     Neuroblast     1   OldBlood -2.240388839
# 11     Neuroblast     2 YoungBlood  1.845594149 <<<
# 12     Neuroblast     2   OldBlood -0.085228863
# 13      Microglia     1 YoungBlood  2.508383068
# 14      Microglia     1   OldBlood -0.180959984
# 15      Microglia     2 YoungBlood -1.419327663
# 16      Microglia     2   OldBlood -0.008436267
# 17    Oligodendro     1 YoungBlood  0.300493250
# 18    Oligodendro     1   OldBlood  0.005350468
# 19    Oligodendro     2 YoungBlood -2.265821659
# 20    Oligodendro     2   OldBlood -0.658590264
# 21    Endothelial     1 YoungBlood  0.816709254
# 22    Endothelial     1   OldBlood -1.831576790
# 23    Endothelial     2 YoungBlood -1.435661188
# 24    Endothelial     2   OldBlood  0.137522805

# aNSC
# (2.579259143 + 0.670768494) / 2
# [1] 1.625014

# Mean of batches for Neuroblast
# (2.1657+ 1.84559) / 2
# [1] 2.00

mean(filter(et, Condition == "OldBlood")$Effect) # -0.8731641
mean(filter(et, Condition == "YoungBlood")$Effect) # 0.5611245


