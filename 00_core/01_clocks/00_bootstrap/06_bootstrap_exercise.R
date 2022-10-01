# REQUIRES COMPUTE CLUSTER TO RUN

library(tidyverse)
library(Seurat)
library(Matrix)

setwd("/labs/abrunet1/Buckley/17.Clock/SvzClockV5/bootstrap/Exercise")
#==================================================================================================
print("Load Exercise 10x Data (v3 and v2 Combined")
pb <- readRDS("data/joint_harmony_SVZ_2020-11-30.rds") # 26 samples, 92326 cells
### NOTE: The above Seurat object contains all cells in "data/ex_seurat.SVZ.annotated.2020-04-27.rds" along with
###       cells from an earlier pilot experiment. In all downstream analyses, we subset to just those cells that
###       are in "data/ex_seurat.SVZ.annotated.2020-04-27.rds" (i.e. Batch=="R2"). In practice, this is identical
###       to running all analyses with "data/ex_seurat.SVZ.annotated.2020-04-27.rds".
pb.umi <- t(pb[["RNA"]]@counts)
pb.genes <- colnames(pb.umi)

#==================================================================================================
print("Load Clock Training Object to match column structure")
svz <- readRDS("../data/bootstrap_pseudocell_15.rds")
genes <- colnames(svz[1,6][[1]][[1]]) # use these genes to subset exercise data. 20948 genes.

#==================================================================================================
print("Format exercise data")
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



# Change Cell type labels
df <- tibble("Celltype"=as.character(pb@meta.data$Celltype),
             "Sample"=as.character(pb@meta.data$AgeCond),
             "Mouse"=as.character(pb@meta.data$Mouse),
             "Year"=as.character(pb@meta.data$Year),
             pb.umi)

# Clean memory
rm(pb); rm(pb.umi); rm(svz)


# Convert to nested tidy form & nest
Celltypes <- c("Microglia", "Endothelial", "Oligodendro",
               "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")
by_celltype <- df %>%
    group_by(Celltype, Mouse, Sample, Year) %>%
    nest() %>%
    dplyr::filter(Celltype %in% Celltypes)
colnames(by_celltype)[4] <- "exData"

#==================================================================================================
print("Bootstrap and pseudobulk")

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

# Apply boostrap.pseudocells using map()
df2 <- by_celltype %>% mutate(pseudocell_ex = map(exData, bootstrap.pseudocells))
print(dim(df2))            
# Remove single cell data; keep just key metadata and pseudocells
df2$exData <- NULL

df2 <- df2 %>% unnest(pseudocell_ex)
print(dim(df2))

saveRDS(df2, "data/bootstrap_pseudocell_15_exercise.rds")
print("Done")







