# Purpose:
# Take a seurat object and return a pseudobulked dataframe
# Technique of pseudobulking
# Bootstrap sampling rather than random partitions
    # Equally weights all samples

library(Seurat)
library(tidyverse)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/01_bootstrap_bioage")
svz <- readRDS("../../00_multiseq_aging/Integrate/data/multi_intergrated_seurat_Dec2020.rds")

Convert_to_Dataframe <- function(svz) {
    DefaultAssay(svz) <- "RNA"
    svz[["SCT"]] <- NULL
    svz[["LMO"]] <- NULL
    Celltypes <- c("Oligodendro", "Microglia", "aNSC_NPC",
                   "Astrocyte_qNSC", "Neuroblast", "Endothelial")
    svz <- subset(svz, subset = Celltype.LowRes %in% Celltypes)
    meta <- svz@meta.data
    meta <- meta[, c("hash.ID", "Prolif_Lineage_Fraction_of_SVZ", "Celltype.LowRes", "orig.ident")]
    raw_counts <- t(as.matrix(svz[["RNA"]]@counts))
    raw_counts <- raw_counts[, colSums(raw_counts) > 0]
    df <- as_tibble(cbind(meta, raw_counts))
    return(df)
}

df <- Convert_to_Dataframe(svz) %>%
        group_by(Celltype.LowRes, Prolif_Lineage_Fraction_of_SVZ, orig.ident, hash.ID) %>%
        nest()

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

# Apply bootstrap.pseudocells using map()
set.seed(42)
df2 <- df %>% mutate(pseudocell_all = map(data, bootstrap.pseudocells))

#==================================================================================================
# Remove single cell data; keep just key metadata and pseudocells
df2$data <- NULL
# Original did not specify set.seed, don't overwrite
#saveRDS(df2, "data/bootstrap_pseudocell_15_bio.rds")
