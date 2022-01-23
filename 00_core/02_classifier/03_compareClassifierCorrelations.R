library(tidyverse)
library(glmnet)
library(ggthemes)
library(viridis)
library(ggridges)
library(ggpubr)
library(broom)
library(ggpubr)

#==================================================================================================
# Load Data
#==================================================================================================
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/02_classifier")

output <- readRDS("data/exercise_classifier_on_multi.rds")
output2 <- readRDS("data/parabiosis_classifier_on_multi.rds")
output2$odds <- NULL
output2$prob <- NULL

output$Intervention <- "Exercise"
output2$Intervention <- "Parabiosis"

combo <- rbind(output, output2)
combo$Intervention <- factor(combo$Intervention, levels = c("Parabiosis", "Exercise"))
Celltypes <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")
combo$Celltype <- factor(combo$Celltype, levels = Celltypes)

ggplot(combo, aes(Age, Pred)) +
    geom_bin2d(bins = 20, alpha = .75) +
    facet_grid(Intervention~Celltype) +
    geom_smooth(method = "lm") +
    geom_hline(yintercept = 0) +
    scale_fill_distiller(palette="PuBuGn", direction = 1, limits = c(0, 90), oob = scales::squish) +
    theme_classic() + ylim(c(-15,15)) +
    ylab("log( Pr(Control)/Pr(Intervention))") +
    stat_cor(size = 2, aes(label = ..r.label..))
ggsave("plots/logOdds_correlations.pdf", width = 6.25, height = 2.73)

#===============================
# Dotplot Version

cor_fun <- function(df) cor.test(df$Age, df$Pred, method = "pearson") %>% tidy()

combo2 <- combo %>% group_by(Celltype, Intervention) %>% nest()
combo2 <- combo2 %>% mutate(model = map(data, cor_fun))
combo2 <- combo2 %>% select(-data) %>% unnest(cols = c(model))
combo2$Intervention <- factor(combo2$Intervention, levels = c("Exercise", "Parabiosis"))

# Dotplots
ggplot(data = combo2, mapping = aes(x = Celltype, y = Intervention, size = -log10(p.value), color = estimate)) +
        geom_point() +
        scale_color_distiller(palette="RdBu", limits = c(-.5,.5), oob = scales::squish) +
        geom_point(shape = 1, color = "black") +
        #facet_grid(.~Celltype) +
        theme_classic() +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), strip.text.x = element_text(size = 7)) +
        theme(legend.position = "right") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(legend.direction = "vertical", legend.box = "horizontal")
ggsave("plots/dotplot_correlations.pdf", width = 4.14 , height = 1.22)


#===============================
# aNSC Highlight

d <- filter(combo, Intervention == "Parabiosis", Celltype == "aNSC_NPC")
ggplot(d, aes(Age, Pred)) +
    geom_bin2d(bins = 20, alpha = .75) +
    geom_smooth(method = "lm") +
    geom_hline(yintercept = 0) +
    scale_fill_distiller(palette="PuBuGn", direction = 1, limits = c(0, 60), oob = scales::squish) +
    theme_classic() +
    ylab("log(odds) Sedentary old")
ggsave("plots/aNSC_parabiosis.pdf", width=2.67, height=1.88)