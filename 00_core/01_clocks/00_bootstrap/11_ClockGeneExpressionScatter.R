# Predict Bootstrap Pseudocells using leave-one-batch-out cross validation
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
df <- readRDS("data/bootstrap_pseudocell_15.rds")
df <- df %>% ungroup %>% select(-c(data))


# Lognormalize counts, including a pseudocount
lognorm <- function(input) {
  norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
  log1p(norm * 10000)
}

df2 <- df %>% mutate(lognorm = map(pseudocell_all, lognorm)) 
df2$pseudocell_all <- NULL
df2 <- unnest(df2, lognorm)
dim(df2) # 16800 20952
df <- NULL # Clear memory


Celltypes <- c("Endothelial", "Astrocyte_qNSC", "aNSC_NPC",
               "Neuroblast", "Oligodendro", "Microglia")


for (i in c(1:6)) {
  ct <- as.character(models[i,1][[1]][[1]])
  df_CT <- filter(df2, Celltype.LowRes == ct)[, -c(1:4)]
  
  # compute mean and variance expressions
  mean_exp_CT <- colMeans(df_CT)
  var_exp_CT <- colStdevs(df_CT) / mean_exp_CT
  perc_nonzero_CT <- colSums(df_CT != 0) / dim(df_CT)[1]
  
  # take log
  log2_mean_exp_CT <- log2(mean_exp_CT)
  log2_var_exp_CT <- log2(var_exp_CT)
  
  
  # select for clock genes only
  lasso <- models[i,2][[1]][[1]]
  selected_genes <- rownames(coef(lasso, s = 'lambda.min'))[coef(lasso, s = 'lambda.min')[,1]!= 0]
  colors <- ifelse(colnames(df_CT) %in% selected_genes, 'Used in Clock', 'Not in Clock')
  alphas <- ifelse(colnames(df_CT) %in% selected_genes, 1.0, 0.01)
  
  # Make dataframe
  res.df<- data.frame(colnames(df_CT), log2_mean_exp_CT, log2_var_exp_CT, perc_nonzero_CT)
  
  
  # Make scatter plot of clock genes
  p <- ggplot(res.df, aes(x=log2_mean_exp_CT, y=log2_var_exp_CT, color=colors)) +
    theme_classic() +
    geom_point(size=2, alpha=alphas) +
    xlab("log2 mean gene expression") + 
    ylab("log2 CV gene expression") + 
    labs(colour = "Gene") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    scale_colour_manual(values = c("black","red"))
  p <- p + theme(text = element_text(size = 20))   
  ggsave(paste0("plots/clock_mean_cv_",ct,".png"), p, width=9.07, height=6.78, dpi=500)
  p
  
  p <- ggplot(res.df, aes(x=log2_mean_exp_CT, y=perc_nonzero_CT, color=colors)) +
    theme_classic() +
    geom_point(size=2, alpha=alphas) +
    xlab("log2 mean gene expression") + 
    ylab("Fraction of cells with nonzero counts") + 
    labs(colour = "Gene") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    scale_colour_manual(values = c("black","red"))
  p <- p + theme(text = element_text(size = 20))   
  ggsave(paste0("plots/clock_mean_percnz_",ct,".png"), p, width=9.07, height=6.78, dpi=500)
  p
  
  p <- ggplot(res.df, aes(x=log2_var_exp_CT, y=perc_nonzero_CT, color=colors)) +
    theme_classic() +
    geom_point(size=2, alpha=alphas) +
    xlab("log2 CV gene expression") + 
    ylab("Fraction of cells with nonzero counts") + 
    labs(colour = "Gene") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    scale_colour_manual(values = c("black","red"))
  p <- p + theme(text = element_text(size = 20))   
  ggsave(paste0("plots/clock_cv_percnz_",ct,".png"), p, width=9.07, height=6.78, dpi=500)
  p
}

