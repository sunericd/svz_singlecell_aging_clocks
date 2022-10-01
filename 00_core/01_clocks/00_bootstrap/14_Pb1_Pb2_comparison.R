# Compare PB 1 and PB 2 clock performances and genes
library(tidyverse)
library(resample)
library(glmnet)
library(viridis)
library(ggpubr)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")

# read in clocks
models <- readRDS("data/models_all_bootstrap.rds")
models$lognormalized <- NULL

# read in sc data
dfp <- readRDS("data/bootstrap_pseudocell_15_parabiosis.rds")


# Lognormalize counts, including a pseudocount
meta_p <- as.data.frame(dfp[, c(1:3)])
umi <- dfp[, -c(1:3)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)
dfp[, -c(1:3)] <- logged

# Separate cohort 1 and 2
dfp1 <- dfp[grepl("2019", dfp$Mouse), ]
dfp2 <- dfp[!grepl("2019", dfp$Mouse), ]

Celltypes <- c("Endothelial", "Astrocyte_qNSC", "aNSC_NPC",
               "Neuroblast", "Oligodendro", "Microglia")

for (i in c(1:6)) {
  ct <- as.character(models[i,1][[1]][[1]])
  dfp1_CT <- filter(dfp1, Celltype == ct)[, -c(1:3)]
  dfp2_CT <- filter(dfp2, Celltype == ct)[, -c(1:3)]
  
  # compute mean and variance expressions
  mean_exp_CT1 <- colMeans(dfp1_CT)
  mean_exp_CT2 <- colMeans(dfp2_CT)

  # take log
  log2_mean_exp_CT1 <- log2(mean_exp_CT1)
  log2_mean_exp_CT2 <- log2(mean_exp_CT2)
  
  
  ##### CLOCK GENES #####
  
  # select for clock genes only
  lasso <- models[i,2][[1]][[1]]
  selected_genes <- rownames(coef(lasso, s = 'lambda.min'))[coef(lasso, s = 'lambda.min')[,1]!= 0]
  colors <- ifelse(colnames(dfp1_CT) %in% selected_genes, 'Used in Clock', 'Not in Clock')
  alphas <- ifelse(colnames(dfp1_CT) %in% selected_genes, 1.0, 0.01)
  
  
  # Make dataframe
  res.df<- data.frame(colnames(dfp1_CT), log2_mean_exp_CT1, log2_mean_exp_CT2)
  
  # Make scatter plot of clock genes
  p <- ggplot(res.df, aes(x=log2_mean_exp_CT2, y=log2_mean_exp_CT1, color=colors)) +
    theme_classic() +
    geom_point(size=2, alpha=alphas) +
    xlab("Parabiosis Cohort 2, log2 mean gene expression") + 
    ylab("Parabiosis Cohort 1, log2 mean gene expression") + 
    labs(colour = "Gene") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    scale_colour_manual(values = c("black","red"))
  p <- p + theme(text = element_text(size = 20))   
  ggsave(paste0("plots/PB12_mean_expression_Clock_",ct,".png"), p, width=9.07, height=6.78, dpi=500)
}





# Violin Plots of Predicted Ages Compared
library(ggthemes)
library(scales)
library(dplyr)
library(tidyr)
library(ggpubr)
library(gridExtra)

d <- readRDS("data/parabiosis_predictions.rds")
d$Batch <- d$Mouse
d$Batch[!grepl("2019", d$Batch)] <- "2"
d$Batch[grepl("2019", d$Batch)] <- "1"
d$Batch <- factor(d$Batch, levels = c("1", "2"))
d$Sample <- factor(d$Sample, levels = c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"))
d$Celltype <- factor(d$Celltype, levels = c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
                                            "Microglia", "Oligodendro", "Endothelial"))
d <- d[d$Sample == "Young-Iso", ]

d <- d %>% group_by(Mouse, Celltype, Sample, Batch) %>% summarise(Pred=median(Pred)) %>% ungroup()



p <-  ggplot(d, aes_string(x="Batch", y="Pred", fill="Batch")) +
  geom_violin(alpha=0.7) +
  geom_point(color="black") +
  geom_boxplot(width=0.1, fill="lightgray") +
  scale_fill_tableau() +
  scale_color_tableau() +
  #scale_fill_manual(values = c("skyblue","orange")) +
  stat_compare_means(aes_string(group = "Sample"), comparisons=list(c("1", "2")), method="wilcox.test") +
  xlab("Condition") +
  facet_wrap(~Celltype, scales="free_y") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
ggsave(paste0("plots/violinplot_PB12.pdf"),p, height=5, width=8)
  

