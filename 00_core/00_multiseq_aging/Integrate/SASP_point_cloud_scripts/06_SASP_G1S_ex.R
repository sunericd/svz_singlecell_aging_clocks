# Make SASP score vs G1S score plot

library(Seurat)
library(ggplot2)
library(dplyr)

filename <- "data/sasp_geneset"
gene_list_GO <- readLines(paste0(filename,".txt"))
genelist <- tolower(c(gene_list_GO))


# Load Data
#svz <- readRDS("../../../02_exercise/data/seurat.SVZ.annotated.2020-04-27.rds")
svz <- readRDS("data/ex_seurat.SVZ.annotated.2020-04-27.rds")
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
meta <- meta[, colnames(meta) %in% c("S.Score", "G2M.Score", "Phase","AgeCond", "Replicate", "Celltype.LowRes")]
adhesion_data$adhesion_response <- rowSums(adhesion_data)
adhesion_data <- cbind(meta, adhesion_data)

# Reorder factors
print(unique(adhesion_data$Celltype.LowRes))
print(unique(adhesion_data$AgeCond))
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "OPC", "Oligodendro", "Endothelial")
adhesion_data <- adhesion_data %>% filter(Celltype.LowRes %in% CELLS)
adhesion_data$celltype <- factor(adhesion_data$Celltype.LowRes,  levels=CELLS, ordered=T)


qNSCs <- dplyr::filter(adhesion_data, celltype == "Oligodendro") 
YC_qNSCs <- dplyr::filter(qNSCs, AgeCond == "Y_Control") 
OC_qNSCs <- dplyr::filter(qNSCs, AgeCond == "O_Control")
YE_qNSCs <- dplyr::filter(qNSCs, AgeCond == "Y_Exercise") 
OE_qNSCs <- dplyr::filter(qNSCs, AgeCond == "O_Exercise") 

aNSCs <- dplyr::filter(adhesion_data, celltype == "aNSC_NPC") 
YC_aNSCs <- dplyr::filter(aNSCs, AgeCond == "Y_Control") 
OC_aNSCs <- dplyr::filter(aNSCs, AgeCond == "O_Control")
YE_aNSCs <- dplyr::filter(aNSCs, AgeCond == "Y_Exercise") 
OE_aNSCs <- dplyr::filter(aNSCs, AgeCond == "O_Exercise") 

Neuroblasts <- dplyr::filter(adhesion_data, celltype == "Neuroblast") 
YC_neuro <- dplyr::filter(Neuroblasts, AgeCond == "Y_Control") 
OC_neuro <- dplyr::filter(Neuroblasts, AgeCond == "O_Control")
YE_neuro <- dplyr::filter(Neuroblasts, AgeCond == "Y_Exercise") 
OE_neuro <- dplyr::filter(Neuroblasts, AgeCond == "O_Exercise") 


YC_NSC_lineage <- rbind(YC_qNSCs, YC_aNSCs, YC_neuro)
YE_NSC_lineage <- rbind(YE_qNSCs, YE_aNSCs, YE_neuro)
OC_NSC_lineage <- rbind(OC_qNSCs, OC_aNSCs, OC_neuro)
OE_NSC_lineage <- rbind(OE_qNSCs, OE_aNSCs, OE_neuro)


theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

g <- ggplot(data=OE_NSC_lineage, aes(x=S.Score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4")) + ylim(0,55) + xlim(-0.2,1) + theme
g <- g + ylab("SASP Score") + xlab("S.Score")
ggsave("plots/scatterplot_SASP_CellCycle_Ex_OldExercise.pdf",g, height=6, width=6)
g

g <- ggplot(data=OC_NSC_lineage, aes(x=S.Score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4"))  + ylim(0,55) + xlim(-0.2,1) + theme
g <- g + ylab("SASP Score") + xlab("S.Score")
ggsave("plots/scatterplot_SASP_CellCycle_Ex_OldControl.pdf",g, height=6, width=6)
g

g <- ggplot(data=YC_NSC_lineage, aes(x=S.Score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4"))  + ylim(0,55) + xlim(-0.2,1) + theme
g <- g + ylab("SASP Score") + xlab("S.Score")
ggsave("plots/scatterplot_SASP_CellCycle_Ex_YoungControl.pdf",g, height=6, width=6)
g

g <- ggplot(data=YE_NSC_lineage, aes(x=S.Score, y=adhesion_response, color=celltype)) + geom_point() + scale_color_manual(values=c("dodgerblue","red2","purple4"))  + ylim(0,55) + xlim(-0.2,1) + theme
g <- g + ylab("SASP Score") + xlab("S.Score")
ggsave("plots/scatterplot_SASP_CellCycle_Ex_YoungExercise.pdf",g, height=6, width=6)
g
