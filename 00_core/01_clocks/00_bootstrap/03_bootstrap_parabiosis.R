# Use Clocks Trained on 28 Mice (v3 10x chem) to predict
# ages of parabiosis cells (v2 and v3 10x chem).

library(tidyverse)
library(Seurat)

#==================================================================================================
# Load Parabiosis 10x Data (v3 and v2 Combined)
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")
pb <- readRDS("../../../01_parabiosis/both/data/pb_combined.rds")
pb.umi <- as_tibble(t(as.matrix((pb[["RNA"]]@counts)))) # 25531 cells 31053 genes
pb.genes <- colnames(pb.umi)

#==================================================================================================
# Load Clock Training Object to match column structure
svz <- readRDS("data/bootstrap_pseudocell_15.rds")
genes <- colnames(svz[1,6][[1]][[1]]) # use these genes to subset parabiosis data. 20948 genes.

#==================================================================================================
# Format parabiosis data
pb.missing <- setdiff(genes, pb.genes)
if (length(as.numeric(pb.missing)) > 0) {
  print("There are some genes in the training data
        that are not in test data. Filling with Zeros.")
  missing_df <- matrix(0, nrow(pb.umi), length(pb.missing))
  colnames(missing_df) <- pb.missing
  rownames(missing_df) <- rownames(pb.umi)
  pb.umi <- cbind(pb.umi, missing_df)
}
pb.umi <- pb.umi[, genes] # Reorder and cut to size to match clock training data format

# Fix metadata in seurat object
pb$AgeCond <- plyr::mapvalues(pb$AgeCond, 
          from=c("yy", "oo", "het_y", "het_o"),
          to=c("Young-Iso", "Old-Iso", "Young-Het", "Old-Het"))
Mouse <- pb@meta.data$hash.ID
Mouse[is.na(Mouse)] <- paste0(pb@meta.data$Experiment[is.na(Mouse)], "_",
                              pb@meta.data$AgeCond[is.na(Mouse)])
pb@meta.data$hash.ID <- Mouse

# Change Cell type labels
df <- tibble("Celltype"=as.character(pb@meta.data$Celltype),
             "Sample"=as.character(pb@meta.data$AgeCond),
             "Mouse"=as.character(pb@meta.data$hash.ID),
             pb.umi)

# Clean memory
rm(pb); rm(pb.umi); rm(svz)


# Convert to nested tidy form & nest
Celltypes <- c("Microglia", "Endothelial", "Oligodendro",
               "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")
by_celltype <- df %>%
    group_by(Celltype, Mouse, Sample) %>%
    nest() %>%
    dplyr::filter(Celltype %in% Celltypes)
colnames(by_celltype)[4] <- "parabiosisData"

#==================================================================================================
# Bootstrap and pseudobulk

bootstrap.pseudocells <- function(df, size=15, n=100, replace="dynamic") {
    pseudocells <- c()
    # If dynamic then only sample with replacement if required due to shortage of cells
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

# Apply boostrap.pseudocells using map(). Takes ~20 minutes
df2 <- by_celltype %>% mutate(pseudocell_parabiosis = map(parabiosisData, bootstrap.pseudocells))
               
# Remove single cell data; keep just key metadata and pseudocells
df2$parabiosisData <- NULL

df2 <- df2 %>% unnest(pseudocell_parabiosis) # A tibble: 13,200 x 20,951

# Don't overwrite
# saveRDS(df2, "data/bootstrap_pseudocell_15_parabiosis.rds")








