# Figure 5 Simple Comparision between Old Exercise Effect, Old Parabiosis Effect
# Bootstrap predictions only
# Average results from cohort 1 and 2 parabiosis


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
# PARABIOSIS
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

et <- et %>% 
  group_by(Celltype, Condition) %>% 
  summarize(Effect = mean(Effect))

et_pb <- et
data.frame(et_pb)

#=============================================================================
# EXERCISE

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
et_ex <- et

#=============================================================================
# COMBINE 
epb <- bind_rows(et_pb, et_ex)

epb$Condition <- factor(epb$Condition, levels = c("YoungExerciseEffect", "OldBlood",
                                            "OldExerciseEffect", "YoungBlood"))

ggplot(epb, aes(x = " ", y = -Effect, fill = Condition)) +
	geom_bar(stat = "identity", position = "dodge") +
	facet_grid(.~Celltype) +
	scale_fill_manual(values = c("#F7BB7F", "#F28E2B", "#acd3d0", "#76B7B2")) +
	theme_cowplot() +
	theme(axis.text.x = element_blank(),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		strip.text.x = element_text(size = 8)) +
	geom_hline(yintercept=0, size = .4)
ggsave("plots/barplot_simpleCombo.pdf", width = 11, height = 2.3)

