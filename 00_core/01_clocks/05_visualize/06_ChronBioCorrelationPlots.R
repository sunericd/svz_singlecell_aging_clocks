library(tidyverse)
library(Seurat)
library(ggpubr)
library(cowplot)


setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/05_visualize")


resultsChron <- readRDS("../00_bootstrap/data/leaveOneBatchPredictions_bootstrapPseudocell.rds")
resultsBio <- readRDS("../01_bootstrap_bioage/data/leaveOneBatchPredictions_bootstrapPC_bio.rds")

#======================================================================================================================
chron.medians <- resultsChron %>% 
  group_by(Celltype, Mouse) %>% 
  mutate(median = median(Prediction, na.rm=TRUE))  %>% select(-Prediction) %>% distinct()

chron.meanmed <- chron.medians %>% 
  group_by(Mouse) %>% 
  mutate(predchron = mean(median, na.rm=TRUE)) %>% select(-median, -Celltype) %>% distinct()



bio.medians <- resultsBio %>% 
  group_by(Celltype, Mouse) %>% 
  mutate(median = median(Prediction, na.rm=TRUE))  %>% select(-Prediction) %>% distinct()

bio.meanmed <- bio.medians %>% 
  group_by(Mouse) %>% 
  mutate(predbio = mean(median, na.rm=TRUE)) %>% select(-median, -Celltype) %>% distinct()


full.meanmed <- left_join(chron.meanmed, bio.meanmed, by = c("Mouse", "Batch"))

# convert prolf frac to bio age
full.meanmed$predbio <- 35 - 100*full.meanmed$predbio
full.meanmed$BioAge <- 35 - 100*full.meanmed$Prolif_Lineage_Fraction_of_SVZ


# make plots

p2 <- ggplot(full.meanmed, aes(x = Age, y = Prolif_Lineage_Fraction_of_SVZ)) +
  geom_point() +
  geom_smooth(method = 'lm', span = 10, formula = 'y ~ x', color="#27285C") +
  stat_cor( aes(label = ..r.label..), method = "spearman", label.x = 10, label.y = .3) +
  theme_pubr() +
  xlab("Chronological Age (Months)") +
  ylab("Proliferative Fraction (%)")
ggsave("plots/Age_vs_Fraction.pdf", width = 4, height = 3, useDingbats=F)

p2 <- ggplot(full.meanmed, aes(x = Age, y = BioAge)) +
  geom_point() +
  geom_smooth(method = 'lm', span = 10, formula = 'y ~ x', color="#27285C") +
  stat_cor( aes(label = ..r.label..), method = "spearman", label.x = 10, label.y = 25) +
  theme_pubr() +
  xlab("Chronological Age (Months)") +
  ylab("Biological Age (Months)")
ggsave("plots/Age_vs_BioAge.pdf", width = 4, height = 3, useDingbats=F)

p2 <- ggplot(full.meanmed, aes(x = predchron, y = predbio)) +
  geom_point() +
  geom_smooth(method = 'lm', span = 10, formula = 'y ~ x', color="#27285C") +
  stat_cor( aes(label = ..r.label..), method = "spearman", label.x = 10, label.y = 25) +
  theme_pubr() +
  xlab("Predicted Chronological Age (Months)") +
  ylab("Predicted Biological Age (Months)")
ggsave("plots/Predicted_Age_vs_BioAge.pdf", width = 4, height = 3, useDingbats=F)

