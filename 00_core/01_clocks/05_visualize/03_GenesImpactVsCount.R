# Sum coeffecients of shared genes vs specific genes for each celltype

library(tidyverse)
library(ggthemes)
library(cowplot)

# Load Coefficients

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/05_visualize")
geneBinaryMatrix <- readRDS("data/models_upset_data.rds")
sharedGenes <- rownames(geneBinaryMatrix[rowSums(geneBinaryMatrix) >  1, ])
all <- readRDS("data/clock_genes_all.rds")


celltypes <- c("Microglia", "Oligodendro", "Endothelial",
                "Astrocyte_qNSC", "aNSC_NPC",  "Neuroblast")
# Non shared
summary <- c()
for (celltype in celltypes) {
    print(celltype)
    all_ct <- filter(all, Celltype == celltype)
    all_specific <- filter(all_ct, ! gene %in% sharedGenes)
    all_shared <- filter(all_ct, gene %in% sharedGenes)
    sum_abs_specific <- sum(abs(all_specific$coef))
    sum_abs_shared <- sum(abs(all_shared$coef))
    count_specific <- length(all_specific$coef)
    count_shared <- length(all_shared$coef)
    avg_specific <- sum_abs_specific/count_specific
    avg_shared <- sum_abs_shared/count_shared
    dataDf <- data.frame("SumSpecific" = sum_abs_specific,
                         "SumShared" = sum_abs_shared,
                         "CountSpecific" = count_specific,
                         "CountShared" = count_shared,
                         "AvgSpecific" = avg_specific,
                         "AvgShared" = avg_shared, 
                         "PercentShared" = count_shared/(count_shared+count_specific),
                         "ImpactShared" = sum_abs_shared/(sum_abs_shared+sum_abs_specific),
                         "Celltype" = celltype)
    summary <- rbind(summary, dataDf)
}

summary2  <- summary %>% select(Celltype, PercentShared)
summary2$PercentNotShared <- 1- summary2$PercentShared
summary3 <- summary2 %>% pivot_longer(cols = starts_with("Perc"), names_to = "Shared", values_to = "Percent")
summary3$Shared <- gsub('Percent', '', summary3$Shared)
summary3$Metric <- "Count"

summary4  <- summary %>% select(Celltype, ImpactShared)
summary4$ImpactNotShared <- 1- summary4$ImpactShared
summary5 <- summary4 %>% pivot_longer(cols = starts_with("Imp"), names_to = "Shared", values_to = "Percent")
summary5$Shared <- gsub('Impact', '', summary3$Shared)
summary5$Metric <- "Impact"

summarylong <- rbind(summary3, summary5)

#summary2$ImpactNotShared <- 1- summary2$ImpactShared

ggplot(summarylong, aes(x = Metric, y = Percent, fill = Shared)) +
    geom_bar(stat = "identity") +
    facet_grid(. ~ Celltype ) +
    scale_fill_tableau() +
    theme_cowplot() +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    coord_flip()
ggsave("plots/impactSharedSpecific.pdf", width = 10.4, height= 1.16)


############################
# Biological Age
############################

geneBinaryMatrix <- readRDS("data/models_upset_data_bio.rds")
sharedGenes <- rownames(geneBinaryMatrix[rowSums(geneBinaryMatrix) >  1, ])
all <- readRDS("data/clock_genes_all_bio.rds")


celltypes <- c("Microglia", "Oligodendro", "Endothelial",
               "Astrocyte_qNSC", "aNSC_NPC",  "Neuroblast")
# Non shared
summary <- c()
for (celltype in celltypes) {
  print(celltype)
  all_ct <- filter(all, Celltype == celltype)
  all_specific <- filter(all_ct, ! gene %in% sharedGenes)
  all_shared <- filter(all_ct, gene %in% sharedGenes)
  sum_abs_specific <- sum(abs(all_specific$coef))
  sum_abs_shared <- sum(abs(all_shared$coef))
  count_specific <- length(all_specific$coef)
  count_shared <- length(all_shared$coef)
  avg_specific <- sum_abs_specific/count_specific
  avg_shared <- sum_abs_shared/count_shared
  dataDf <- data.frame("SumSpecific" = sum_abs_specific,
                       "SumShared" = sum_abs_shared,
                       "CountSpecific" = count_specific,
                       "CountShared" = count_shared,
                       "AvgSpecific" = avg_specific,
                       "AvgShared" = avg_shared, 
                       "PercentShared" = count_shared/(count_shared+count_specific),
                       "ImpactShared" = sum_abs_shared/(sum_abs_shared+sum_abs_specific),
                       "Celltype" = celltype)
  summary <- rbind(summary, dataDf)
}

summary2  <- summary %>% select(Celltype, PercentShared)
summary2$PercentNotShared <- 1- summary2$PercentShared
summary3 <- summary2 %>% pivot_longer(cols = starts_with("Perc"), names_to = "Shared", values_to = "Percent")
summary3$Shared <- gsub('Percent', '', summary3$Shared)
summary3$Metric <- "Count"

summary4  <- summary %>% select(Celltype, ImpactShared)
summary4$ImpactNotShared <- 1- summary4$ImpactShared
summary5 <- summary4 %>% pivot_longer(cols = starts_with("Imp"), names_to = "Shared", values_to = "Percent")
summary5$Shared <- gsub('Impact', '', summary3$Shared)
summary5$Metric <- "Impact"

summarylong <- rbind(summary3, summary5)

#summary2$ImpactNotShared <- 1- summary2$ImpactShared

ggplot(summarylong, aes(x = Metric, y = Percent, fill = Shared)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Celltype ) +
  scale_fill_tableau() +
  theme_cowplot() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  coord_flip()
ggsave("plots/impactSharedSpecific_bio.pdf", width = 10.4, height= 1.16)


