# Make SASP score vs G1S score plot

library(Seurat)
library(ggplot2)
library(dplyr)

filename <- "data/sasp_geneset"
gene_list_GO <- readLines(paste0(filename,".txt"))
genelist <- tolower(c(gene_list_GO))


# Load Data
svz <- readRDS("data/multi_intergrated_seurat_Dec2020.rds")
meta <- svz[[]]
d <- t(as.matrix(svz[['RNA']]@counts))

rm(svz)
gc()

# Normalize counts to expression values
d <- sweep(d, MARGIN = 1, FUN = "/", STATS = rowSums(d))
d <- log1p(d * 10000)
d <- as.data.frame(d)


# Subset data to signature genes
colnames(d) <- tolower(colnames(d))
adhesion_data <- d[, colnames(d) %in% genelist]
meta <- meta[, colnames(meta) %in% c("S.Score", "G2m.Score", "Phase","Age", "Replicate", "Celltype.LowRes")]
adhesion_data$adhesion_response <- rowSums(adhesion_data)
adhesion_data <- cbind(meta, adhesion_data)

# Reorder factors
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Neuron",
           "OPC", "Oligodendro", "Endothelial", "Ependymal", "Mural",
           "Microglia", "Macrophage")
adhesion_data$celltype <- factor(adhesion_data$Celltype.LowRes,  levels=CELLS, ordered=T)
adhesion_data <- adhesion_data[(adhesion_data$Age>20 & adhesion_data$Age<24)  | (adhesion_data$Age>3 & adhesion_data$Age<5),]
adhesion_data$age[adhesion_data$Age>20 & adhesion_data$Age<24] = "o"
adhesion_data$age[adhesion_data$Age>3 & adhesion_data$Age<5] = "y"
adhesion_data$age <- factor(adhesion_data$age, levels=c("y", "o"), ordered=T)


qNSCs <- dplyr::filter(adhesion_data, celltype == "Oligodendro") #1597 qNSCs/Astrocytes
Y_qNSCs <- dplyr::filter(qNSCs, age == "y") #1117
O_qNSCs <- dplyr::filter(qNSCs, age == "o") #480

aNSCs <- dplyr::filter(adhesion_data, celltype == "aNSC_NPC") #866 aNSCs/NPCs
Y_aNSCs <- dplyr::filter(aNSCs, age == "y") #784
O_aNSCs <- dplyr::filter(aNSCs, age == "o") #82

Neuroblasts <- dplyr::filter(adhesion_data, celltype == "Neuroblast") #1014 Neuroblasts
Y_neuro <- dplyr::filter(Neuroblasts, age == "y") #868
O_neuro <- dplyr::filter(Neuroblasts, age == "o") #146


young_NSC_lineage <- rbind(Y_qNSCs, Y_aNSCs, Y_neuro)
old_NSC_lineage <- rbind(O_qNSCs, O_aNSCs, O_neuro)
whole_NSC_lineage <- rbind(young_NSC_lineage, old_NSC_lineage)


theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

g <- ggplot(data=old_NSC_lineage, aes(x=S.Score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4")) + ylim(0,55) + xlim(-0.2,1) + theme
g <- g + ylab("SASP Score") + xlab("S.Score")
ggsave("plots/scatterplot_SASP_CellCycle_Buckley_Old.pdf",g, height=6, width=6)
g

g <- ggplot(data=young_NSC_lineage, aes(x=S.Score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4"))  + ylim(0,55) + theme
g <- g + ylab("SASP Score") + xlab("S.Score")
ggsave("plots/scatterplot_SASP_CellCycle_Buckley_Young.pdf",g, height=6, width=6)
g
