library(vip)
library(caret)
library(glmnet)
library(tidyverse)
library(tidyr)
library(ggthemes)
library(ggrepel)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/05_visualize")


#############################
# Chronological Age
#############################


models <- readRDS("../00_bootstrap/data/models_all_bootstrap.rds")
models$lognormalized <- NULL
dge <- readRDS("../../00_multiseq_aging/Integrate/data/age_de_df_April2021.rds")


celltypes = c("Oligodendro", "Microglia", "Endothelial",
              "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

# Apply cut offs for * marker
fdr_cut <- 0.1
fold_cut <- log(1.1)


# Plots
for (i in c(1:6)) {
  ct <- as.character(models[i,1][[1]][[1]])
  print(ct)
  dge_ct <- dge[dge$celltype==ct,]
  lasso <- models[i,2][[1]][[1]]
  print(lasso)
  c1 = vip::vi_model(lasso$glmnet.fit, s = lasso$lambda.min, method = "shap")
  d <- c1 %>% filter(Importance != 0)
  
  
  d_POS <- d %>% filter(Sign == "POS")
  d_NEG <- d %>% filter(Sign == "NEG")
  
  # Make clock gene labels -- separate into up down clock genes
  dge_ct <- dge_ct %>% mutate(clockgene = ifelse(gene %in% d_POS$Variable, "Positive coefficient in clock", ifelse(gene %in% d_NEG$Variable, "Negative coefficient in clock", "Not used by clock"))) 
  dge_ct$delabel <- NA
  dge_ct$delabel[dge_ct$clockgene != "Not used by clock"] <- dge_ct$gene[dge_ct$clockgene != "Not used by clock"]
  
  
  # Do the volcano plot
  p <- ggplot(data=dge_ct, aes(x=avg_logFC, y=-log10(fdr), col=dge_ct$clockgene, label=delabel)) +
    geom_point() + 
    theme_classic() +
    geom_text_repel(fontface=3) +
    scale_color_manual(values=c("#4E79A7", "gray", "#F28E2B")) +
    geom_vline(xintercept=c(-fold_cut, fold_cut), col="black") +
    geom_hline(yintercept=-log10(0.05), col="black") +
    labs(y= "-log10(FDR)", x = "Average Log Fold Change")
  ggsave(paste0("plots/volcano_", ct, ".png"), p, width=7.07, height=6.78, dpi=500)
}



#############################
# Biological Age
#############################


models <- readRDS("../01_bootstrap_bioage/data/models_all_bootstrapCell_bio.rds")
models$lognormalized <- NULL
dge <- readRDS("../../00_multiseq_aging/Integrate/data/bioage_de_df_March2022.rds")


celltypes = c("Oligodendro", "Microglia", "Endothelial",
              "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

# Apply cut offs for * marker
fdr_cut <- 0.1
fold_cut <- log(1.1)


# Plots
for (i in c(1:6)) {
  ct <- as.character(models[i,1][[1]][[1]])
  print(ct)
  dge_ct <- dge[dge$celltype==ct,]
  lasso <- models[i,2][[1]][[1]]
  print(lasso)
  c1 = vip::vi_model(lasso$glmnet.fit, s = lasso$lambda.min, method = "shap")
  d <- c1 %>% filter(Importance != 0)
  
  d$Sign[d$Sign == 'POS'] <- 'placeholder'
  d$Sign[d$Sign == 'NEG'] <- 'POS'
  d$Sign[d$Sign == 'placeholder'] <- 'NEG'
  
  
  d_POS <- d %>% filter(Sign == "POS")
  d_NEG <- d %>% filter(Sign == "NEG")
  
  # Make clock gene labels -- separate into up down clock genes
  dge_ct <- dge_ct %>% mutate(clockgene = ifelse(gene %in% d_POS$Variable, "Positive coefficient in clock", ifelse(gene %in% d_NEG$Variable, "Negative coefficient in clock", "Not used by clock"))) 
  dge_ct$delabel <- NA
  dge_ct$delabel[dge_ct$clockgene != "Not used by clock"] <- dge_ct$gene[dge_ct$clockgene != "Not used by clock"]
  
  dge_ct$avg_logFC <- log(2^(dge_ct$avg_log2FC))
  
  
  # Do the volcano plot
  p <- ggplot(data=dge_ct, aes(x=avg_logFC, y=-log10(fdr), col=dge_ct$clockgene, label=delabel)) +
    geom_point() + 
    theme_classic() +
    geom_text_repel(fontface=3) +
    scale_color_manual(values=c("#4E79A7", "gray", "#F28E2B")) +
    geom_vline(xintercept=c(-fold_cut, fold_cut), col="black") +
    geom_hline(yintercept=-log10(fdr_cut), col="black") +
    labs(y= "-log10(FDR)", x = "Average Log Fold Change")
  ggsave(paste0("plots/volcano_", ct, "_bio.png"), p, width=7.07, height=6.78, dpi=500)
}



