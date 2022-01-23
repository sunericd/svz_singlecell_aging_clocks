# Visualize overlap of genes selected by l1 regression models to predict age in 
# main SVZ celltypes. 

library(tidyverse)
library(UpSetR)
library(glmnet)
library(viridis)
library(ComplexHeatmap)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/05_visualize")
models <- readRDS("../00_bootstrap/data/models_all_bootstrap.rds")
models$lognormalized <- NULL

celltypes = c("Oligodendro", "Microglia", "Endothelial",
            "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

#======================================================================================================================
# Save genes
all <- c()
for (celltype in celltypes) {
    print(celltype)
    clock <- filter(models, Celltype.LowRes == celltype)[1,2][[1]][[1]]
    # Coefficients
    betas <- coef(clock, s = clock$lambda.min) # includes intercept, length 13194
    genes <- unlist(betas@Dimnames[1]) # includes intercept, length 13194
    clock.genes <- genes[which(betas !=0)]
    print(paste0("Genes in model: ", length(clock.genes)))
    clock.betas <- betas[which(betas!=0)]
    clock.df <- data.frame("gene"=clock.genes[-1], "coef"=clock.betas[-1], "Celltype" = celltype)
    clock.df <- clock.df[order(clock.df$coef, decreasing=T),]
    write.csv(clock.df, file = paste0("data/clock_genes_", celltype, ".csv"), quote = F)
    all <- rbind(all, clock.df)
}
saveRDS(all, "data/clock_genes_all.rds")

#======================================================================================================================
# Read and combine genes into named list
celltypeGeneList <- vector(mode="list", length=length(celltypes))
names(celltypeGeneList) <- celltypes
for (celltype in celltypes) {
    print(celltype)
    genes.csv <- read_csv(file=paste0("data/clock_genes_", celltype, ".csv"))
    genes <- genes.csv$gene
    celltypeGeneList[[celltype]] <- genes
}

# Find union of genes
geneUnion = c()
for (celltype in celltypes) {
    geneUnion = union(geneUnion, celltypeGeneList[[celltype]])
    print(length(geneUnion))
}

celltypeGeneBinary <- vector(mode="list", length=length(celltypes))
names(celltypeGeneBinary) <- celltypes
for (celltype in celltypes) {
    gBinary <- vector(length = length(geneUnion))
    for (index in c(1:length(gBinary))) {
        gBinary[index] <- ifelse((geneUnion[index] %in% celltypeGeneList[[celltype]]), 1, 0) 
    }
    gBinary_df <- data.frame(gBinary, row.names=geneUnion)
    colnames(gBinary_df) <- celltype
    celltypeGeneBinary[[celltype]] <- gBinary_df
}


# Combine list of dataframe
celltypeDF <- dplyr::bind_cols(celltypeGeneBinary)
rownames(celltypeDF) <- geneUnion
celltypeDF[1:5,1:5]
#       Oligodendro Microglia Endothelial Astrocyte_qNSC aNSC_NPC
# C4b             1         0           1              1        0
# Apoe            1         1           0              0        0
# Gstm1           1         0           0              0        0
# Eif3f           1         0           0              1        0
# Pmp22           1         0           0              0        0

saveRDS(celltypeDF, "data/models_upset_data.rds")


#======================================================================================================================
# Upset Plot using complexheatmap
#========================================================
celltypeDF <- readRDS("data/models_upset_data.rds")

pdf("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/05_visualize/plots/Complex_Upset.pdf", width = 8, height = 3)
m = make_comb_mat(celltypeDF)
ht = draw(UpSet(m))
od = column_order(ht)
cs = comb_size(m)

decorate_annotation("intersection_size", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
        default.units = "native", just = c("left", "bottom"), gp = gpar(fontsize = 8), rot = 45)
})
dev.off()

#======================================================================================================================
# For Manual Annotation of upset plot
celltypeDF[rowSums(celltypeDF) == 3, ]
# H2-D1, Tmem176b, Rsrp1, Il33, C4b, Tma7, Rbbp7, Ctnnbip1, Scd2, Pcdhb7, Ckb
celltypeDF[rowSums(celltypeDF) == 4, ] # H1f0, Gm42418, Hist1h1c, Sparc
celltypeDF[rowSums(celltypeDF) == 5, ] # Ifi27
celltypeDF[rowSums(celltypeDF) == 6, ] # AC149090.1

generalG <- rownames(celltypeDF[rowSums(celltypeDF) >  1, ])
cat(generalG, sep = '\n')

v_generalG <- rownames(celltypeDF[rowSums(celltypeDF) >  2, ])
cat(v_generalG, sep = '\n')

