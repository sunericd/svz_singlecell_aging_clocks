# Dotplots both methods Exercise

library(ggthemes)
library(cowplot)
library(tidyverse)
library(scales)


setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/01_bootstrap_bioage")
d <- readRDS("data/exercise_predictions_bio.rds")
cell_order <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")
d$Batch <- as.character(d$Year)
d$Batch <- factor(d$Batch, levels = c("R1", "R2"))
d$Sample <- factor(d$Sample, levels = c("Y_Control", "Y_Exercise", "O_Control", "O_Exercise"))
d$Celltype <- factor(d$Celltype, levels = cell_order)
df <- d
df$Year <- NULL
df <- filter(df, Batch == "R2")
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
# 1  Astrocyte_qNSC   OldExerciseEffect -1.0131277
# 2  Astrocyte_qNSC YoungExerciseEffect  0.1537159
# 3        aNSC_NPC   OldExerciseEffect  1.1563850 <<<
# 4        aNSC_NPC YoungExerciseEffect  1.3645865 <<<
# 5      Neuroblast   OldExerciseEffect -0.4150609
# 6      Neuroblast YoungExerciseEffect  0.7959546
# 7       Microglia   OldExerciseEffect -1.1608973
# 8       Microglia YoungExerciseEffect -0.4094212
# 9     Oligodendro   OldExerciseEffect  0.7736392 <<<
# 10    Oligodendro YoungExerciseEffect  0.5577027 <<<
# 11    Endothelial   OldExerciseEffect -0.3133071
# 12    Endothelial YoungExerciseEffect -0.3992197
et$Clock <- "Bootstrap"
et_bootstrap <- et

#=============================================================================
setwd("../02_ensemble")
d <- read.csv("data/exercise_predicted_fraction.csv")
d <- as_tibble(d)

d$Batch <- as.character(d$hash_id)
d$Batch[!grepl("_R1_", d$Batch)] <- "R2"
d$Batch[grepl("_R1_", d$Batch)] <- "R1"
d$Batch <- factor(d$Batch, levels = c("R1", "R2"))
d$Sample <- factor(d$condition, levels = c("Y_Control", "Y_Exercise", "O_Control", "O_Exercise"))
d$Celltype <- factor(d$celltype, levels = cell_order)
d$Pred <- d$predicted_age
d$Mouse <- d$hash_id
df <- d
df$Year <- NULL
df <- filter(df, Batch == "R2")
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
        # Young Ex effect:
        d <- filter(SampleMeans,Celltype == celltype)
        eff <- d[d$Sample == "O_Control", "MedianPseudoCellAge"] - d[d$Sample == "O_Exercise", "MedianPseudoCellAge"]
        effect_table <- rbind(effect_table, c(celltype,"OldExerciseEffect", eff))
        # Old Ex effects
        eff2 <- d[d$Sample == "Y_Control", "MedianPseudoCellAge"] - d[d$Sample == "Y_Exercise", "MedianPseudoCellAge"]
        effect_table <- rbind(effect_table, c(celltype,"YoungExerciseEffect", eff2))
    }

colnames(effect_table) <- c("Celltype", "Condition", "Effect")
effect_table
#  [1,] "Oligodendro"    "OldExerciseEffect"   -0.1403197
#  [2,] "Oligodendro"    "YoungExerciseEffect" 0.7127449 
#  [3,] "Microglia"      "OldExerciseEffect"   -1.217294 
#  [4,] "Microglia"      "YoungExerciseEffect" -0.7554955
#  [5,] "Endothelial"    "OldExerciseEffect"   -0.2790609
#  [6,] "Endothelial"    "YoungExerciseEffect" 0.09314584
#  [7,] "Astrocyte_qNSC" "OldExerciseEffect"   -1.153985 
#  [8,] "Astrocyte_qNSC" "YoungExerciseEffect" 0.5066315 
#  [9,] "aNSC_NPC"       "OldExerciseEffect"   0.5002919 
# [10,] "aNSC_NPC"       "YoungExerciseEffect" 1.099566  
# [11,] "Neuroblast"     "OldExerciseEffect"   -0.4284901
# [12,] "Neuroblast"     "YoungExerciseEffect" 0.9236344 

et <- data.frame(effect_table) %>% unnest(cols = c(Celltype, Condition, Effect))
et$Celltype <- factor(et$Celltype, levels = cell_order)
data.frame(et)

#          Celltype           Condition      Effect
# 1     Oligodendro   OldExerciseEffect -0.14031966
# 2     Oligodendro YoungExerciseEffect  0.71274492
# 3       Microglia   OldExerciseEffect -1.21729395
# 4       Microglia YoungExerciseEffect -0.75549551
# 5     Endothelial   OldExerciseEffect -0.27906092
# 6     Endothelial YoungExerciseEffect  0.09314584
# 7  Astrocyte_qNSC   OldExerciseEffect -1.15398472
# 8  Astrocyte_qNSC YoungExerciseEffect  0.50663151
# 9        aNSC_NPC   OldExerciseEffect  0.50029195
# 10       aNSC_NPC YoungExerciseEffect  1.09956561
# 11     Neuroblast   OldExerciseEffect -0.42849012
# 12     Neuroblast YoungExerciseEffect  0.92363442

et$Clock <- "Ensemble"


#=============================================================================
#=============================================================================

et2 <- rbind(et, et_bootstrap)

et2$Condition <- factor(et2$Condition, levels = c("YoungExerciseEffect", "OldExerciseEffect"))

ggplot(data = et2, mapping = aes(x = Clock, y = Condition, size = abs(Effect), color = -Effect)) +
        geom_point() +
        scale_color_distiller(palette="RdBu", limits = c(-2.5,2.5), oob = scales::squish) +
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
ggsave("plots/exercise_dotplot_combo_bio.pdf", useDingbats=F, width = 4, height = 1.8)
