# Dotplot Overview of Effects (Chronological-trained)
# Ensemble and Bootstrap methods on same plot
# Parabiosis then Exercise

library(ggthemes)
library(cowplot)
library(tidyverse)
library(scales)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")
d <- readRDS("data/parabiosis_predictions.rds")
cell_order <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

d$Batch <- d$Mouse
d$Batch[!grepl("2019", d$Batch)] <- "2"
d$Batch[grepl("2019", d$Batch)] <- "1"
d$Batch <- factor(d$Batch, levels = c("1", "2"))
d$Sample <- factor(d$Sample, levels = c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"))
d$Celltype <- factor(d$Celltype, levels = cell_order)
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
data.frame(et)

et$CondBatch <- paste0(et$Condition, "_", et$Batch)
et$CondBatch <- factor(et$CondBatch, levels = c("YoungBlood_2", "YoungBlood_1",
                                                "OldBlood_2", "OldBlood_1"))
et$Clock <- "bootstrap"

et_bootstrap <- et
data.frame(et)

#====================================================================
# Ensemble Data

setwd("../02_ensemble")
d <- read.csv("data/parabiosis_predicted_age.csv")
d <- as_tibble(d)

d$Batch <- d$batch
d$Batch[!grepl("2019", d$Batch)] <- "2"
d$Batch[grepl("2019", d$Batch)] <- "1"
d$Batch <- factor(d$Batch, levels = c("1", "2"))
d$Sample <- factor(d$condition, levels = c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"))
d$Celltype <- factor(d$celltype, levels = cell_order)
d$Pred <- d$predicted_age
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
data.frame(et)


et$CondBatch <- paste0(et$Condition, "_", et$Batch)
et$CondBatch <- factor(et$CondBatch, levels = c("YoungBlood_2", "YoungBlood_1",
                                                "OldBlood_2", "OldBlood_1"))

et$Clock <- "ensemble"

#=============================================================================

et2 <- rbind(et, et_bootstrap)

ggplot(data = et2, mapping = aes(x = Clock, y = CondBatch, size = abs(Effect), color = -Effect)) +
        geom_point() +
        scale_color_distiller(palette="RdBu", limits = c(-1,1)*max(abs(3)), oob = scales::squish) +
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
            labels = c("0","1","2","3","5"),
            guide = "legend",
            limits = c(0, 5),
            oob = scales::squish)
ggsave("plots/dotplot_parabiosis_bootstrap_ensemble_avg.pdf", useDingbats=F, width = 4, height = 2.35)




#=============================================================================
################## Exercise
#=============================================================================

pb <- et2
setwd("../00_bootstrap")
d <- readRDS("data/exercise_predictions.rds")

d$Batch <- as.character(d$Year)
d$Batch <- factor(d$Batch, levels = c("R1", "R2"))
d$Sample <- factor(d$Sample, levels = c("Y_Control", "Y_Exercise", "O_Control", "O_Exercise"))
d$Celltype <- factor(d$Celltype, levels = cell_order)
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
effect_table

et <- data.frame(effect_table) %>% unnest(cols = c(Celltype,Condition, Effect))
et$Celltype <- factor(et$Celltype, levels = cell_order)
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
# 9     Oligodendro   OldExerciseEffect  1.99692100
# 10    Oligodendro YoungExerciseEffect  1.37495884
# 11    Endothelial   OldExerciseEffect  0.36484821
# 12    Endothelial YoungExerciseEffect  1.22334163

et$Clock <- "Bootstrap"
et_bootstrap <- et



#=============================================================================
setwd("../02_ensemble")
d <- read.csv("data/exercise_predicted_age.csv")
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
et <- data.frame(effect_table) %>% unnest(cols = c(Celltype, Condition, Effect))
et$Celltype <- factor(et$Celltype, levels = cell_order)
data.frame(et)
#          Celltype           Condition      Effect
# 1  Astrocyte_qNSC   OldExerciseEffect -1.54759711
# 2  Astrocyte_qNSC YoungExerciseEffect  0.47872655
# 3        aNSC_NPC   OldExerciseEffect -0.03601384
# 4        aNSC_NPC YoungExerciseEffect  1.39384384
# 5      Neuroblast   OldExerciseEffect -1.41773407
# 6      Neuroblast YoungExerciseEffect  0.38159798
# 7       Microglia   OldExerciseEffect -0.83904288
# 8       Microglia YoungExerciseEffect -0.14718335
# 9     Oligodendro   OldExerciseEffect  1.13291453
# 10    Oligodendro YoungExerciseEffect  1.20853807
# 11    Endothelial   OldExerciseEffect -0.22244969
# 12    Endothelial YoungExerciseEffect  0.85626383

et$Clock <- "Ensemble"

#=============================================================================

et2 <- rbind(et, et_bootstrap)

et2$Condition <- factor(et2$Condition, levels = c("YoungExerciseEffect", "OldExerciseEffect"))

ggplot(data = et2, mapping = aes(x = Clock, y = Condition, size = abs(Effect), color = -Effect)) +
        geom_point() +
        scale_color_distiller(palette="RdBu", limits = c(-1,1)*max(abs(2.5)), oob = scales::squish) +
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
ggsave("plots/dotplot_exercise_combo.pdf", useDingbats=F, width = 4, height = 1.8)
