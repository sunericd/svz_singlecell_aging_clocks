# Predict Bootstrap Pseudocells using leave-one-batch-out cross validation
library(tidyverse)
library(resample)
library(glmnet)
library(viridis)
library(ggpubr)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")

fdr_cut <- 0.1
fold_cut <- log(1.1)

# read in clocks
models <- readRDS("data/models_all_bootstrap.rds")
models$lognormalized <- NULL

# DGE results
dgep <- readRDS("../../00_multiseq_aging/Integrate/data/pb_both_de_df_Sep2021.rds")
dgee <- readRDS("../../00_multiseq_aging/Integrate/data/ex_de_df_April2021.rds")


# read in sc data
dfe <- readRDS("data/bootstrap_pseudocell_15_exercise.rds")
dfp <- readRDS("data/bootstrap_pseudocell_15_parabiosis.rds")


# Lognormalize counts, including a pseudocount


meta_e <- as.data.frame(dfe[, c(1:4)])
umi <- dfe[, -c(1:4)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)
dfe[, -c(1:4)] <- logged

meta_p <- as.data.frame(dfp[, c(1:3)])
umi <- dfp[, -c(1:3)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)
dfp[, -c(1:3)] <- logged


Celltypes <- c("Endothelial", "Astrocyte_qNSC", "aNSC_NPC",
               "Neuroblast", "Oligodendro", "Microglia")


for (i in c(1:6)) {
  ct <- as.character(models[i,1][[1]][[1]])
  dfe_CT <- filter(dfe, Celltype == ct)[, -c(1:4)]
  dfp_CT <- filter(dfp, Celltype == ct)[, -c(1:3)]
  
  # compute mean and variance expressions
  mean_exp_CTe <- colMeans(dfe_CT)
  mean_exp_CTp <- colMeans(dfp_CT)

  # take log
  log2_mean_exp_CTe <- log2(mean_exp_CTe)
  log2_mean_exp_CTp <- log2(mean_exp_CTp)
  
  
  
  ##### CLOCK GENES #####
  
  # select for clock genes only
  lasso <- models[i,2][[1]][[1]]
  selected_genes <- rownames(coef(lasso, s = 'lambda.min'))[coef(lasso, s = 'lambda.min')[,1]!= 0]
  colors <- ifelse(colnames(dfe_CT) %in% selected_genes, 'Used in Clock', 'Not in Clock')
  alphas <- ifelse(colnames(dfe_CT) %in% selected_genes, 1.0, 0.01)
  
  
  # Make dataframe
  res.df<- data.frame(colnames(dfe_CT), log2_mean_exp_CTe, log2_mean_exp_CTp)
  
  # Make scatter plot of clock genes
  p <- ggplot(res.df, aes(x=log2_mean_exp_CTp, y=log2_mean_exp_CTe, color=colors)) +
    theme_classic() +
    geom_point(size=2, alpha=alphas) +
    xlab("Parabiosis, log2 mean gene expression") + 
    ylab("Exercise, log2 mean gene expression") + 
    labs(colour = "Gene") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    scale_colour_manual(values = c("black","red"))
  p <- p + theme(text = element_text(size = 20))   
  ggsave(paste0("plots/PB_EX_mean_expression_Clock_",ct,".png"), p, width=9.07, height=6.78, dpi=500)
  
  
  
  ##### DEGS #####
  
  # select for clock genes only
  exercise_genes <- dgee %>% filter(celltype == ct & abs(avg_logFC) > fold_cut & fdr < fdr_cut)
  parabiosis_genes <- dgep %>% filter(celltype == ct & abs(avg_logFC) > fold_cut & fdr < fdr_cut)
  colors <- ifelse(colnames(dfe_CT) %in% exercise_genes$gene, "Exercise", ifelse(colnames(dfe_CT) %in% parabiosis_genes$gene, "Parabiosis", ""))
  alphas <- ifelse(colnames(dfe_CT) %in% exercise_genes$gene, 1.0, ifelse(colnames(dfe_CT) %in% parabiosis_genes$gene, 1.0, 0.01))
  
  # Make dataframe
  res.df<- data.frame(colnames(dfe_CT), log2_mean_exp_CTe, log2_mean_exp_CTp)
  
  # Make scatter plot of clock genes
  p <- ggplot(res.df, aes(x=log2_mean_exp_CTp, y=log2_mean_exp_CTe, color=colors)) +
    theme_classic() +
    geom_point(size=2, alpha=alphas) +
    xlab("Parabiosis, log2 mean gene expression") + 
    ylab("Exercise, log2 mean gene expression") + 
    labs(colour = "Gene") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    scale_colour_manual(values = c("black","blue","green"))
  p <- p + theme(text = element_text(size = 20))   
  ggsave(paste0("plots/PB_EX_mean_expression_DEG_",ct,".png"), p, width=9.07, height=6.78, dpi=500)
  
}

