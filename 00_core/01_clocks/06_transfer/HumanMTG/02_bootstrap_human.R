library(tidyverse)
library(Matrix)
library(Seurat)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/06_transfer/HumanMTG")


#==================================================================================================
# Load Human MTG data
pb <- readRDS("human_mtg_seurat.rds")
celltypes <- c("Astro L1-2 FGFR3 GFAP", "Astro L1-6 FGFR3 SLC14A1", "Endo L2-6 NOSTRIN", "Micro L1-3 TYROBP",
               "Oligo L1-6 OPALIN")
pb <- subset(x=pb, subset = cluster %in% celltypes)
saveRDS(pb, "human_mtg_seurat_subset.rds")

pb.umi <- as_tibble(t(as.matrix((pb[["RNA"]]@counts))))
colnames(pb.umi) <- tolower(colnames(pb.umi))
pb.genes <- colnames(pb.umi)


#==================================================================================================
# Load Clock Training Object to match column structure
svz <- readRDS("data/bootstrap_pseudocell_15.rds")
genes <- tolower(colnames(svz[1,6][[1]][[1]])) # use these genes to subset data. 20948 genes.

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
pb@meta.data$cluster <- plyr::mapvalues(pb@meta.data$cluster, 
                              from=c("Astro L1-2 FGFR3 GFAP", "Astro L1-6 FGFR3 SLC14A1",
                                     "Endo L2-6 NOSTRIN", "Micro L1-3 TYROBP",
                                     "Oligo L1-6 OPALIN"),
                              to=c("Astrocyte_qNSC", "Astrocyte_qNSC",
                                   "Endothelial", "Microglia",
                                   "Oligodendro"))

# Change Cell type labels
df <- tibble("Celltype"=as.character(pb@meta.data$cluster),
             "Donor"=as.character(pb@meta.data$donor),
             "Age"=pb@meta.data$age_days,
             pb.umi)

# Clean memory
rm(pb); rm(pb.umi); rm(svz)
gc()

# Convert to nested tidy form & nest
Celltypes <- c("Microglia", "Endothelial", "Oligodendro",
               "Astrocyte_qNSC")
by_celltype <- df %>%
  group_by(Celltype, Donor, Age) %>%
  nest() %>%
  dplyr::filter(Celltype %in% Celltypes)
colnames(by_celltype)[4] <- "humanData"

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
df2 <- by_celltype %>% mutate(pseudocell_human = map(humanData, bootstrap.pseudocells))
      
# Remove single cell data; keep just key metadata and pseudocells
df2$humanData <- NULL

df2 <- df2 %>% unnest(pseudocell_human) 

# Don't overwrite
saveRDS(df2, "data/bootstrap_pseudocell_15_human.rds")








