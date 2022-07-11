# Histograms of predictions for 
# Bootstrap, parabiosis then exercise

library(ggthemes)
library(cowplot)
library(tidyverse)
library(scales)
library(ggridges)

CELLTYPES <- c("Oligodendro", "Microglia",  "Endothelial",
			   "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

#=====================================================================
#=====================================================================
## BootstrapCell Predictions- Parabiosis

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")
d <- readRDS("data/parabiosis_predictions.rds")
d$Batch <- d$Mouse
d$Batch[!grepl("2019", d$Batch)] <- "2"
d$Batch[grepl("2019", d$Batch)] <- "1"
d$Batch <- factor(d$Batch, levels = c("1", "2"))
d$Sample <- factor(d$Sample, levels = c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"))
d$Celltype <- factor(d$Celltype, levels = CELLTYPES)
df <- filter(d, Celltype == "aNSC_NPC")

##################
# aNSC Highlight
##################
ggplot(df, aes(x=Pred, fill=Sample, color=Sample)) +
 geom_density(alpha = 0.7) +
 facet_grid(rows=vars(Batch)) +
 scale_fill_tableau() +
 scale_color_tableau() +
 theme_classic() +
 xlim(c(0, 30)) + ylim(c(0, 0.3)) +
 ylab("Density") + xlab("Predicted Age")
ggsave("plots/parabiosis_aNSC_histo.pdf", width=4.81, height=3.11)


mice <- rev(c("2019_Young-Iso", "BC12-Young-Iso-5.67-10B-B2-70100-AAGGCTAG",
        "BC11-Young-Iso-5.67-10A-B2-70700-GAGTCGAT", 
        "BC23-Young-Iso-4.60-15A-B4-28000-ATCTACGG", "BC24-Young-Iso-4.60-15B-B4-22800-TGTACCAG",
        "2019_Young-Het", "BC2-Young-Het-5.10-1B-B1-14300-TCAATGGC",
        "BC14-Young-Het-5.13-5B-B2-50100-AACCGAAC", "BC8-Young-Het-5.13-2B-B2-30500-GAAGCTTG",
        "BC22-Young-Het-5.40-14B-B4-38000-CGATTAGC", "2019_Old-Het",
        "BC1-Old-Het-20.50-1A-B1-25000-TGTGATGG", "BC13-Old-Het-20.53-5A-B2-50000-CAGTTAGG",
        "BC7-Old-Het-20.53-2A-B2-51000-GTACCTGT", "BC21-Old-Het-20.80-14A-B4-24000-GAGAGACT",
        "2019_Old-Iso", "BC6-Old-Iso-20.50-9B-B1-29863-CGAACAAG",
        "BC5-Old-Iso-20.50-9A-B1-24388-AGTTGCGT", "BC3-Old-Iso-20.50-6A-B1-35260-CTCTAGAC",
        "BC4-Old-Iso-20.50-6B-B1-16224-ACCAATGC", "BC9-Old-Iso-20.53-7A-B2-33700-AAGTACGC",
        "BC10-Old-Iso-20.53-7B-B2-21000-ATTCGCAC"))

df$Mouse <- factor(df$Mouse, levels = mice)

df$Cohort <- df$Batch
ggplot(df, aes(x = Pred, y = Mouse, color=Sample, fill=Sample))+
        geom_density_ridges(alpha = .95, size = .3) +
        scale_fill_tableau() +
        scale_color_tableau() +
        theme_classic()
ggsave("plots/parabiosis_histo_facs_mouse_bothCohorts_aNSC.pdf", width=6.76, height=3.4)


##################
# All
##################

ggplot(filter(d, Batch == 1), aes(x=Pred, fill=Sample, color=Sample)) +
 geom_density(alpha = 0.7) +
 facet_wrap(.~Celltype) +
 scale_fill_tableau() +
 scale_color_tableau() +
 theme_classic() +
 xlim(c(0, 30)) + ylim(c(0, 0.3)) +
 ylab("Density") + xlab("Predicted Age")
ggsave("plots/parabiosis1_histo_all.pdf", width=6.74 , height=3.09)

ggplot(filter(d, Batch == 2), aes(x=Pred, fill=Sample, color=Sample)) +
 geom_density(alpha = 0.7) +
 facet_wrap(.~Celltype) +
 scale_fill_tableau() +
 scale_color_tableau() +
 theme_classic() +
 xlim(c(0, 30)) + ylim(c(0, 0.3)) +
 ylab("Density") + xlab("Predicted Age")
ggsave("plots/parabiosis2_histo_all.pdf", width=6.74 , height=3.09)


#=====================================================================
#=====================================================================
## BootstrapCell Predictions- Exercise

d <- readRDS("data/exercise_predictions.rds")

d$Batch <- as.character(d$Year)
d$Batch <- factor(d$Batch, levels = c("R1", "R2"))
d$Sample <- factor(d$Sample, levels = c("Y_Control", "Y_Exercise", "O_Control", "O_Exercise"))
d$Celltype <- factor(d$Celltype, levels = CELLTYPES)


##################
# Oligodendrocyte Highlight
##################
ggplot(filter(d, Batch =="R2", Celltype == "Oligodendro") , aes(x=Pred, fill=Sample, color=Sample)) +
 geom_density(alpha = 0.7) +
 scale_fill_tableau() +
 scale_color_tableau() +
 theme_classic() +
 xlim(c(0, 30)) + ylim(c(0, 0.25)) +
 ylab("Density") + xlab("Predicted Age")
ggsave("plots/exercise_histo_oligo.pdf", width=5, height=2.1)

##################
# All
##################
ggplot(filter(d, Batch =="R2"), aes(x=Pred, fill=Sample, color=Sample)) +
 geom_density(alpha = 0.7) +
 facet_wrap(.~Celltype) +
 scale_fill_tableau() +
 scale_color_tableau() +
 theme_classic() +
 xlim(c(0, 30)) + ylim(c(0, 0.3)) +
 ylab("Density") + xlab("Predicted Age")
ggsave("plots/exercise_histo_all.pdf", width=6.74 , height=3.09)










