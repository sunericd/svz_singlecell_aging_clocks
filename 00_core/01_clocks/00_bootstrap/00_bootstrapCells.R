# Purpose:
# Take a seurat object and return a bootstrapcell dataframe
# Technique of cell type specific pseudo bulking to optimize
# trade off between cell count and cell complexity
# Bootstrap sampling rather than random partitions
# Equally weights samples rather than cells

library(Seurat)
library(tidyverse)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")
svz <- readRDS("../../00_multiseq_aging/Integrate/data/multi_intergrated_seurat_Dec2020.rds")

Convert_to_Dataframe <- function(svz) {
    DefaultAssay(svz) <- "RNA"
    svz[["SCT"]] <- NULL
    svz[["LMO"]] <- NULL
    Celltypes <- c("Oligodendro", "Microglia", "aNSC_NPC",
                   "Astrocyte_qNSC", "Neuroblast", "Endothelial")
    svz <- subset(svz, subset = Celltype.LowRes %in% Celltypes)
    meta <- svz@meta.data
    meta <- meta[, c("hash.ID", "Age", "Celltype.LowRes", "orig.ident")]
    raw_counts <- t(as.matrix(svz[["RNA"]]@counts))
    raw_counts <- raw_counts[, colSums(raw_counts) > 0]
    df <- as_tibble(cbind(meta, raw_counts))
    return(df)
}

df <- Convert_to_Dataframe(svz) %>%
        group_by(Celltype.LowRes, Age, orig.ident, hash.ID) %>%
        nest()

head(df)
#   <chr>            <dbl> <fct>           <chr>      <list>                 
# 1 Sample4-ACCAATGC  29   Oligodendro     Batch-1    <tibble [452 x 20,948]>
# 2 Sample4-ACCAATGC  29   Endothelial     Batch-1    <tibble [46 x 20,948]> 
# 3 Sample4-ACCAATGC  29   Microglia       Batch-1    <tibble [284 x 20,948]>
# 4 Sample1-TGTGATGG   4.7 Oligodendro     Batch-1    <tibble [295 x 20,948]>
# 5 Sample3-CTCTAGAC  20.8 Microglia       Batch-1    <tibble [226 x 20,948]>
# 6 Sample2-TCAATGGC   6.7 Oligodendro     Batch-1    <tibble [468 x 20,948]>


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
set.seed(42)
df2 <- df %>% mutate(pseudocell_all = map(data, bootstrap.pseudocells)) # ~15 minutes on a laptop

#==================================================================================================
# Remove single cell data; keep just key metadata and pseudocells
df2$data <- NULL
saveRDS(df2, "data/bootstrap_pseudocell_15_seed42.rds")

# Original did not specify set.seed, don't overwrite
#saveRDS(df2, "data/bootstrap_pseudocell_15.rds")

