# Purpose:
# Take a seurat object and return a bootstrapCell dataframe
# Bootstrap sampling rather than random partitions
    # Equally weights all samples

library(Seurat)
library(tidyverse)
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/06_transfer/DBN2019")

load("data/svz_celltypes_2019-01-31.rda") #V2 seurat object useding in DBN2019 paper
svz <- UpdateSeuratObject(svz)

Convert_to_Dataframe <- function(svz) {
    DefaultAssay(svz) <- "RNA"
    Celltypes <- c("Oligodendrocytes", "Microglia", "aNSCs_NPCs",
                   "Astrocytes_qNSCs", "Neuroblasts", "Endothelial")
    svz <- subset(svz, subset = Celltype %in% Celltypes)
    meta <- svz@meta.data
    meta <- meta[, c("Age", "Celltype", "orig.ident")]
    raw_counts <- t(as.matrix(svz[["RNA"]]@counts))
    raw_counts <- raw_counts[, colSums(raw_counts) > 0]
    df <- as_tibble(cbind(meta, raw_counts))
    return(df)
}

df <- Convert_to_Dataframe(svz)
df <- df %>% group_by(Celltype, Age, orig.ident) %>% nest()

#===================================================================================================

bootstrap.pseudocells <- function(df, size=15, n=100, replace="dynamic") {
    pseudocells <- c()
    # If dynamic then only sample with replacement if required due to shortage of cells.
    if (replace == "dynamic") {
        if (nrow(df) <= size) {replace <- TRUE} else {replace <- FALSE}
    }
    for (i in c(1:n)) {
        batch <- df[sample(1:nrow(df), size = size, replace = replace), ]
        pseudocells <- rbind(pseudocells, colSums(batch))
    }
    colnames(pseudocells) <- colnames(df)
    return(as_tibble(pseudocells))
}
    
# Apply boostrap.pseudocells using map()
df2 <- df %>% mutate(pseudocell_all = map(data, bootstrap.pseudocells))

#==================================================================================================
# Remove single cell data; keep just key metadata and pseudocells
df2$data <- NULL
saveRDS(df2, "data/bootstrap_pseudocell_15.rds")
