library(tidyverse)
library(glmnet)
library(ggthemes)
library(ggpubr)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/01_bootstrap_bioage")
d <- readRDS("data/exercise_predictions_bio.rds")
cell_order <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

d$Batch <- as.character(d$Year)
d$Batch <- factor(d$Batch, levels = c("R1", "R2"))
d$Sample <- factor(d$Sample, levels = c("Y_Control", "Y_Exercise", "O_Control", "O_Exercise"))
d$Celltype <- factor(d$Celltype, levels = cell_order)
df <- d
df$Year <- NULL
df <- filter(df, Mouse != "O5_R1")
df$Pred <- 35 - 100*df$Pred
df <- filter(df, Batch == "R2")

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
        # Young Ex effect:
    d <- filter(SampleMeans, Celltype == celltype)
    eff <- d[d$Sample == "O_Control", "MedianPseudoCellAge"] - d[d$Sample == "O_Exercise", "MedianPseudoCellAge"]
    effect_table <- rbind(effect_table, c(celltype, "OldExerciseEffect", eff))
    # Old Ex effects
    eff2 <- d[d$Sample == "Y_Control", "MedianPseudoCellAge"] - d[d$Sample == "Y_Exercise", "MedianPseudoCellAge"]
    effect_table <- rbind(effect_table, c(celltype, "YoungExerciseEffect", eff2))
}
colnames(effect_table) <- c("Celltype","Condition", "Effect")
et <- data.frame(effect_table) %>% unnest(cols = c(Celltype,Condition, Effect))
et$Celltype <- factor(et$Celltype, levels = cell_order)
data.frame(et)
#          Celltype           Condition     Effect
# 1     Oligodendro   OldExerciseEffect  0.7736392 <<<
# 2     Oligodendro YoungExerciseEffect  0.5577027 <<<
# 3       Microglia   OldExerciseEffect -1.1608973
# 4       Microglia YoungExerciseEffect -0.4094212
# 5     Endothelial   OldExerciseEffect -0.3133071
# 6     Endothelial YoungExerciseEffect -0.3992197
# 7  Astrocyte_qNSC   OldExerciseEffect -1.0131277
# 8  Astrocyte_qNSC YoungExerciseEffect  0.1537159
# 9        aNSC_NPC   OldExerciseEffect  1.1563850 <<<
# 10       aNSC_NPC YoungExerciseEffect  1.3645865 <<<
# 11     Neuroblast   OldExerciseEffect -0.4150609
# 12     Neuroblast YoungExerciseEffect  0.7959546

