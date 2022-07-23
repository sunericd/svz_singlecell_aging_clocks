library(tidyverse)
library(glmnet)
library(viridis)
library(ggpubr)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/06_transfer/HumanMTG")


d <- readRDS("data/humans_predictions.rds")


celltype_order <- c("Oligodendro", "Astrocyte_qNSC")
d <- d %>% filter(Celltype %in% celltype_order)


# Oligodendro rescale
a1 <- 0.1
b1 <- -0.55
c1 <- 120

# Astrocyte rescale
a2 <- 0.2
b2 <- -0.55
c2 <- 30

# linear rescale
d$Pred <- ifelse(d$Celltype=="Oligodendro",(d$Pred-b1)/a1 + c1, (d$Pred-b2)/a2 + c2)

d$Age <- as.double(d$Age) / 365
d2 <- d %>% group_by(Celltype, Donor) %>% mutate(med = median(Pred)) %>% select(-Pred) %>% distinct()

d$Celltype <- factor(d$Celltype, levels = celltype_order)
d2$Celltype <- factor(d2$Celltype, levels = celltype_order)


ggplot(d, aes(x = Age, y = Pred)) +
  geom_bin2d(binwidth = 1.5) +
  scale_fill_distiller(palette="PuBuGn", direction = 1, limits = c(0, 70), oob = scales::squish) +
  facet_wrap(Celltype~.) +
  theme_classic() +
  xlab("Chronological Age") +
  ylab("Rescaled Predicted Age") +
  geom_point(data = d2, aes(x=Age, y=med), size = 1.5, color = "black") +
  geom_smooth(data = d2, aes(x=Age, y=med), method = "lm", color = "black", size = 1) +
  stat_cor(data = d2, aes(x=Age, y=med), y.label=70, size = 2.7) # pearson r not rho
ggsave("predicted_ages_human.pdf", width = 6.36, height = 3.29, useDingbats = F)
