# Make gene set score validation plots

library(Seurat)
library(ggplot2)
library(ggthemes)
library(scales)
library(dplyr)
library(tidyr)
library(ggpubr)
library(gridExtra)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Integrate")


# gene sets and labels
filenames <- c("data/ifn_gamma_response_genes_msigdb", "data/GO_term_negative_regulation_of_neurogenesis")
filename_labels <- c("Interferon", "Negative Regulation of Neurogenesis")
single_genes <- c()

all_labels <- c(filename_labels, single_genes)

# Load Data
svz <- readRDS("data/pb_combined_2.rds")
meta <- svz[[]]
d <- t(as.matrix(svz[['RNA']]@counts))

rm(svz)
gc()

# Normalize counts to expression values
d <- sweep(d, MARGIN = 1, FUN = "/", STATS = rowSums(d))
d <- log1p(d * 10000)
d <- as.data.frame(d)
colnames(d) <- tolower(colnames(d))

# Subset data to signature genes
counter <- 1
for (fn in filenames){
  gene_list_GO <- readLines(paste0(fn,".csv"))
  genelist <- tolower(c(gene_list_GO))
  sub_data <- d[, colnames(d) %in% genelist]
  d[filename_labels[counter]] <- rowSums(sub_data, na.rm=TRUE)
  counter <- counter + 1
  
  print(paste0("Identified ", dim(sub_data)[2], " genes out of ", length(genelist)))
}


meta <- meta[, colnames(meta) %in% c("S.Score", "G2M.Score", "Phase","AgeCond", "Replicate", "Celltype")]
meta$SScore <- meta$S.Score
meta$G2MScore <- meta$G2M.Score
d <- cbind(meta, d)

# Reorder factors
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Neuron",
           "OPC", "Oligodendro", "Endothelial", "Ependymal", "Mural",
           "Microglia", "Macrophage")
d <- d %>% filter(Celltype %in% CELLS)
d$celltype <- factor(d$Celltype,  levels=CELLS, ordered=T)
d$AgeCond <- plyr::mapvalues(d$AgeCond, 
                             from=c("yy", "het_y", "oo", "het_o"),
                             to=c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"))
d$AgeCond <- factor(d$AgeCond,  levels=c("Young-Iso", "Young-Het", "Old-Iso", "Old-Het"), ordered=T)



####### Make plot
comparisons <- list(c("Young-Iso", "Old-Iso"), c("Old-Iso", "Old-Het"))


for (ct in c("aNSC_NPC", "Oligodendro")){
  p2 <- list()
  dct <- d %>% filter(celltype == ct) %>% pivot_longer(cols = all_labels, names_to = "Variable", values_to = "Value")
  p <-  ggplot(dct, aes_string(x="AgeCond", y="Value", fill="AgeCond")) +
    geom_violin() +
    geom_boxplot(width=0.1, fill="lightgray") +
    scale_fill_tableau() +
    scale_color_tableau() +
    stat_compare_means(aes_string(group = "AgeCond"), comparisons=comparisons, method="wilcox.test") +
    xlab("Condition") +
    facet_wrap(~Variable, scales="free_y") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  ggsave(paste0("plots/gene_set_validation/violinplot_Parabiosis_",ct,".pdf"),p, height=5, width=8)
}