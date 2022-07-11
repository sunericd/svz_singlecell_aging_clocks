# Make gene set score validation plots

library(Seurat)
library(ggplot2)
library(ggthemes)
library(scales)
library(dplyr)
library(tidyr)
library(ggpubr)
library(gridExtra)


# gene sets and labels
filenames <- c("data/ifn_gamma_response_genes_msigdb", "data/neurogenesis_regulation")
filename_labels <- c("Interferon", "Regulation of Neurogenesis")
single_genes <- c()

# filenames <- c("data/up_senescence_genes", "data/casella_up_senescence", "data/nagano_up_senescence",
#                "data/sasp_geneset", "data/GO_Cell_Adhesion_GO0007155", "data/ifn_gamma_response_genes_msigdb",
#                "data/generation_neurons", "data/neurogenesis", "data/neurogenesis_regulation")
# filename_labels <- c("CellAge_UP_Senescence", "Casella_UP_Senescence", "Nagano_UP_Senescence", "Reactome_SASP", "GO_Cell_Adhesion", "mSigDB_IFN_gamma_response", "generation_neurons", "neurogenesis", "neurogenesis_regulation")
# single_genes <- c("cdkn2a", "mki67", "SScore", "G2MScore")

all_labels <- c(filename_labels, single_genes)

# Load Data
svz <- readRDS("data/ex_seurat.SVZ.annotated.2020-04-27.rds")
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


meta <- meta[, colnames(meta) %in% c("S.Score", "G2M.Score", "Phase","AgeCond", "Replicate", "Celltype.LowRes")]
meta$SScore <- meta$S.Score
meta$G2MScore <- meta$G2M.Score
d <- cbind(meta, d)

# Reorder factors
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "OPC", "Oligodendro", "Endothelial")
d <- d %>% filter(Celltype.LowRes %in% CELLS)
d$celltype <- factor(d$Celltype.LowRes,  levels=CELLS, ordered=T)
d$AgeCond <- factor(d$AgeCond,  levels=c("Y_Control", "Y_Exercise", "O_Control", "O_Exercise"), ordered=T)


####### Make plot
comparisons <- list(c("Y_Control", "O_Control"), c("O_Control", "O_Exercise"))
#d <- d %>% filter(AgeCond %in% c("Y_Control", "O_Control"))
#comparisons <- list(c("Y_Control", "O_Control"))


for (ct in c("aNSC_NPC", "Oligodendro")){
  p2 <- list()
  dct <- d %>% filter(celltype == ct) %>% pivot_longer(cols = all_labels, names_to = "Variable", values_to = "Value")
  p <-  ggplot(dct, aes_string(x="AgeCond", y="Value", fill="AgeCond")) +
    geom_violin(alpha=0.7) +
    geom_boxplot(width=0.1, fill="lightgray") +
    scale_fill_tableau() +
    scale_color_tableau() +
    #scale_fill_manual(values = c("skyblue","orange")) +
    stat_compare_means(aes_string(group = "AgeCond"), comparisons=comparisons, method="wilcox.test") +
    xlab("Condition") +
    facet_wrap(~Variable, scales="free_y") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  ggsave(paste0("plots/gene_set_validation/violinplot_Exercise_",ct,".pdf"),p, height=5, width=8)
}