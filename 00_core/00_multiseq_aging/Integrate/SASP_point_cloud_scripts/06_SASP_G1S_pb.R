# Make SASP score vs G1S score plot

library(Seurat)
library(ggplot2)
library(dplyr)

filename <- "data/sasp_geneset"
gene_list_GO <- readLines(paste0(filename,".txt"))
genelist <- tolower(c(gene_list_GO))


# Load Data
#svz <- readRDS("../../../01_parabiosis/both/data/pb_combined.rds")
svz <- readRDS("data/pb_combined_2.rds")
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
meta <- meta[, colnames(meta) %in% c("S.Score", "G2M.Score", "Phase","AgeCond", "Replicate", "Celltype")]
adhesion_data$adhesion_response <- rowSums(adhesion_data)
adhesion_data <- cbind(meta, adhesion_data)

# Reorder factors
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Neuron",
           "OPC", "Oligodendro", "Endothelial", "Ependymal", "Mural",
           "Microglia", "Macrophage")
adhesion_data$celltype <- factor(adhesion_data$Celltype,  levels=CELLS, ordered=T)


qNSCs <- dplyr::filter(adhesion_data, celltype == "Oligodendro") 
YI_qNSCs <- dplyr::filter(qNSCs, AgeCond == "Young-Iso") 
OI_qNSCs <- dplyr::filter(qNSCs, AgeCond == "Old-Iso")
YH_qNSCs <- dplyr::filter(qNSCs, AgeCond == "Young-Het") 
OH_qNSCs <- dplyr::filter(qNSCs, AgeCond == "Old-Het") 

aNSCs <- dplyr::filter(adhesion_data, celltype == "aNSC_NPC") 
YI_aNSCs <- dplyr::filter(aNSCs, AgeCond == "Young-Iso") 
OI_aNSCs <- dplyr::filter(aNSCs, AgeCond == "Old-Iso")
YH_aNSCs <- dplyr::filter(aNSCs, AgeCond == "Young-Het") 
OH_aNSCs <- dplyr::filter(aNSCs, AgeCond == "Old-Het") 

Neuroblasts <- dplyr::filter(adhesion_data, celltype == "Neuroblast") 
YI_neuro <- dplyr::filter(Neuroblasts, AgeCond == "Young-Iso") 
OI_neuro <- dplyr::filter(Neuroblasts, AgeCond == "Old-Iso")
YH_neuro <- dplyr::filter(Neuroblasts, AgeCond == "Young-Het") 
OH_neuro <- dplyr::filter(Neuroblasts, AgeCond == "Old-Het") 


YI_NSC_lineage <- rbind(YI_qNSCs, YI_aNSCs, YI_neuro)
YH_NSC_lineage <- rbind(YH_qNSCs, YH_aNSCs, YH_neuro)
OI_NSC_lineage <- rbind(OI_qNSCs, OI_aNSCs, OI_neuro)
OH_NSC_lineage <- rbind(OH_qNSCs, OH_aNSCs, OH_neuro)


theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

g <- ggplot(data=OH_NSC_lineage, aes(x=S.Score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4")) + ylim(0,55) + xlim(-0.2,1) + theme
g <- g + ylab("SASP Score") + xlab("S.Score")
ggsave("plots/scatterplot_SASP_CellCycle_Pb_OldHet.pdf",g, height=6, width=6)
g

g <- ggplot(data=OI_NSC_lineage, aes(x=S.Score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4"))  + ylim(0,55) + xlim(-0.2,1) + theme
g <- g + ylab("SASP Score") + xlab("S.Score")
ggsave("plots/scatterplot_SASP_CellCycle_Pb_OldIso.pdf",g, height=6, width=6)
g

g <- ggplot(data=YI_NSC_lineage, aes(x=S.Score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4"))  + ylim(0,55) + xlim(-0.2,1) + theme
g <- g + ylab("SASP Score") + xlab("S.Score")
ggsave("plots/scatterplot_SASP_CellCycle_Pb_YoungIso.pdf",g, height=6, width=6)
g

g <- ggplot(data=YH_NSC_lineage, aes(x=S.Score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4"))  + ylim(0,55) + xlim(-0.2,1) + theme
g <- g + ylab("SASP Score") + xlab("S.Score")
ggsave("plots/scatterplot_SASP_CellCycle_Pb_YoungHet.pdf",g, height=6, width=6)
g
