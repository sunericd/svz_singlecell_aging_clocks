library(ggthemes)
library(cowplot)
library(tidyverse)
library(scales)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")
d <- readRDS("data/exercise_predictions.rds")

d$Batch <- as.character(d$Year)
d$Batch <- factor(d$Batch, levels = c("R1", "R2"))
d$Sample <- factor(d$Sample, levels = c("Y_Control", "Y_Exercise", "O_Control", "O_Exercise"))
d$Celltype <- factor(d$Celltype, levels = c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
                                            "Microglia", "Oligodendro", "Endothelial"))
df <- d
df$Year <- NULL
df <- filter(df, Batch == "R2")


#=============================================================================
MouseMedians <- df %>%
    group_by(Celltype, Sample, Batch, Mouse) %>%
    summarise(MedianMouseAges = median(Pred))

SampleMeans <- MouseMedians %>% 
    group_by(Celltype, Sample, Batch) %>%
    summarise(MedianPseudoCellAge = mean(MedianMouseAges))

#==============================================================================
# Means
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
colnames(effect_table) <- c("Celltype", "Condition", "Effect")
effect_table

et <- data.frame(effect_table) %>% unnest(cols = c(Celltype, Condition, Effect))
et$Celltype <- factor(et$Celltype, levels = c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
                                            "Microglia", "Oligodendro", "Endothelial"))

data.frame(et)
#          Celltype           Condition      Effect
# 1  Astrocyte_qNSC   OldExerciseEffect -1.45533382
# 2  Astrocyte_qNSC YoungExerciseEffect  0.15588236
# 3        aNSC_NPC   OldExerciseEffect  0.31358232
# 4        aNSC_NPC YoungExerciseEffect  1.87492641
# 5      Neuroblast   OldExerciseEffect -1.06053999
# 6      Neuroblast YoungExerciseEffect  1.01475809
# 7       Microglia   OldExerciseEffect -0.83890869
# 8       Microglia YoungExerciseEffect  0.03227218
# 9     Oligodendro   OldExerciseEffect  1.99692100 <<<
# 10    Oligodendro YoungExerciseEffect  1.37495884 <<<
# 11    Endothelial   OldExerciseEffect  0.36484821
# 12    Endothelial YoungExerciseEffect  1.22334163
