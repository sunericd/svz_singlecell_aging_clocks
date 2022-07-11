# Differential Expression 
# Multiseq, Parabiosis, Exercise

library(tidyverse)
library(Seurat)
library(MAST)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Integrate")

celltypes = c("Microglia", "Endothelial", "Oligodendro",
            "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")


#########################################################################
#===========================================================
# Change between young and old
#===========================================================

svz <- readRDS("data/multi_intergrated_seurat_Dec2020.rds")
svz2 <- subset(svz, subset = Age < 7 | Age > 20)
AgeClass <- svz2@meta.data$Age
AgeClass[AgeClass > 20] <- "Old"
AgeClass[AgeClass != "Old"] <- "Young"
svz2$AgeClass <- AgeClass
DefaultAssay(svz2) <- "RNA"
svz2[["LMO"]] <- NULL
svz2[["SCT"]] <- NULL

# Assign as main identity
svz2$Celltype_AgeClass <- paste0(svz2@meta.data$Celltype.LowRes, "_", svz2@meta.data$AgeClass)
Idents(svz2) <- svz2[["Celltype_AgeClass"]]

de_list <- vector(mode="list", length = length(celltypes))
names(de_list) = celltypes

#
for (CELLTYPE in celltypes) {
    print(CELLTYPE)

    # Find  cluster marker genes if possible

        # Change arguments to Find Markers as desired.
    obj_de <- FindMarkers(object = svz2,
                          test.use = "MAST",
                          ident.1 = paste0(CELLTYPE, "_Old"),
                          ident.2 = paste0(CELLTYPE, "_Young"),
                          #max.cells.per.ident = 1000,
                          min.pct = 0,
                          logfc.threshold = 0,
                          random.seed = 3)

    obj_de <- as.data.frame(obj_de)
    obj_de <- rownames_to_column(obj_de, var = "gene")
    obj_de$fdr <- p.adjust(obj_de$p_val, method = "fdr", n = length(rownames(svz2[["RNA"]]@data)))
    obj_de$celltype <- CELLTYPE
    obj_de$comparision <- "Old/Young"
    obj_de$test <- "mast"
    de_list[[CELLTYPE]] <- obj_de
}

# Combine all matrices into one dataframe
de_df <- data.frame()
for (CELLTYPE in celltypes) {
    print(CELLTYPE)
    de_df  <- rbind(de_df, de_list[[CELLTYPE]])
}
de_df <- tbl_df(de_df)

saveRDS(de_df, paste0("data/age_de_df_April2021.rds"))




#########################################################################
#===========================================================
# Change between young (bio) and old (bio)
#===========================================================

svz <- readRDS("data/multi_intergrated_seurat_Dec2020.rds")
svz@meta.data$BioAge <- 35-100*svz[[]]$Prolif_Lineage_Fraction_of_SVZ
svz2 <- subset(svz, subset = BioAge < 7 | Age > 20)
AgeClass <- svz2@meta.data$BioAge
AgeClass[AgeClass > 20] <- "Old"
AgeClass[AgeClass != "Old"] <- "Young"
svz2$AgeClass <- AgeClass
DefaultAssay(svz2) <- "RNA"
svz2[["LMO"]] <- NULL
svz2[["SCT"]] <- NULL

# Assign as main identity
svz2$Celltype_AgeClass <- paste0(svz2@meta.data$Celltype.LowRes, "_", svz2@meta.data$AgeClass)
Idents(svz2) <- svz2[["Celltype_AgeClass"]]

de_list <- vector(mode="list", length = length(celltypes))
names(de_list) = celltypes

#
for (CELLTYPE in celltypes) {
  print(CELLTYPE)
  
  # Find  cluster marker genes if possible
  
  # Change arguments to Find Markers as desired.
  obj_de <- FindMarkers(object = svz2,
                        test.use = "MAST",
                        ident.1 = paste0(CELLTYPE, "_Old"),
                        ident.2 = paste0(CELLTYPE, "_Young"),
                        #max.cells.per.ident = 1000,
                        min.pct = 0,
                        logfc.threshold = 0,
                        random.seed = 3)
  
  obj_de <- as.data.frame(obj_de)
  obj_de <- rownames_to_column(obj_de, var = "gene")
  obj_de$fdr <- p.adjust(obj_de$p_val, method = "fdr", n = length(rownames(svz2[["RNA"]]@data)))
  obj_de$celltype <- CELLTYPE
  obj_de$comparision <- "Old/Young"
  obj_de$test <- "mast"
  de_list[[CELLTYPE]] <- obj_de
}

# Combine all matrices into one dataframe
de_df <- data.frame()
for (CELLTYPE in celltypes) {
  print(CELLTYPE)
  de_df  <- rbind(de_df, de_list[[CELLTYPE]])
}
de_df <- tbl_df(de_df)

saveRDS(de_df, paste0("data/bioage_de_df_March2022.rds"))












#########################################################################
#===========================================================
# Change between Old Iso to Old Het # Different ages across cohorts!
#===========================================================

obj <- readRDS("../../../01_parabiosis/both/data/pb_combined.rds")

# Fix metadata in object seurat object
obj$AgeCond <- plyr::mapvalues(obj$AgeCond, 
      from=c("yy", "oo", "het_y", "het_o"),
      to=c("Young-Iso", "Old-Iso", "Young-Het", "Old-Het"))

# Append age to celltype label and add to metadata
obj[["CelltypeAgeCond"]] <- paste0(obj@meta.data$Celltype, "_", obj@meta.data$AgeCond)
Idents(obj) <- obj[["CelltypeAgeCond"]]

#===========================
# Combined Parabiosis Cohorts

de_list <- vector(mode="list", length = length(celltypes))
names(de_list) = celltypes
# Compare Exercise vs Control for each celltype, save each matrix in a list. 
for (CELLTYPE in celltypes) {
    print(CELLTYPE)

    # Find  cluster marker genes if possible

        # Change arguments to Find Markers as desired.
    obj_de <- FindMarkers(object = obj,
                          test.use = "MAST",
                          ident.1 = paste0(CELLTYPE, "_Old-Iso"),
                          ident.2 = paste0(CELLTYPE, "_Old-Het"),
                          min.pct = 0,
                          logfc.threshold = 0,
                          random.seed = 3)

    obj_de <- as.data.frame(obj_de)
    obj_de <- rownames_to_column(obj_de, var = "gene")
    obj_de$fdr <- p.adjust(obj_de$p_val, method = "fdr", n = length(rownames(obj2[["RNA"]]@data)))
    obj_de$celltype <- CELLTYPE
    obj_de$comparision <- "OldIso/OldHet"
    obj_de$test <- "mast"
    de_list[[CELLTYPE]] <- obj_de
}

# Combine all matrices into one dataframe
de_df <- data.frame()
for (CELLTYPE in celltypes) {
    print(CELLTYPE)
    de_df  <- rbind(de_df, de_list[[CELLTYPE]])
}

saveRDS(de_df, paste0("data/pb_both_de_df_Sep2021.rds"))

#########################################################################
#===========================================================
# Change between Old Sedentary and Old Exercise
#===========================================================

obj <- readRDS("../../../02_exercise/data/seurat.SVZ.annotated.2020-04-27.rds")

obj$Celltype_AgeCond <- paste0(obj@meta.data$Celltype.LowRes, "_", obj@meta.data$AgeCond)

Idents(obj) <- obj[["Celltype_AgeCond"]]

de_list <- vector(mode="list", length = length(celltypes))
names(de_list) = celltypes


for (CELLTYPE in celltypes) {
    print(CELLTYPE)

    # Find  cluster marker genes if possible

        # Change arguments to Find Markers as desired.
    obj_de <- FindMarkers(object = obj,
                          test.use = "MAST",
                          ident.1 = paste0(CELLTYPE, "_O_Control"),
                          ident.2 = paste0(CELLTYPE, "_O_Exercise"),
                          min.pct = 0,
                          logfc.threshold = 0,
                          random.seed = 3)

    obj_de <- as.data.frame(obj_de)
    obj_de <- rownames_to_column(obj_de, var = "gene")
    obj_de$fdr <- p.adjust(obj_de$p_val, method = "fdr", n = length(rownames(obj[["RNA"]]@data)))
    obj_de$celltype <- CELLTYPE
    obj_de$comparision <- "OldSed/OldEx"
    obj_de$test <- "mast"
    de_list[[CELLTYPE]] <- obj_de
}

# Combine all matrices into one dataframe
de_df <- data.frame()
for (CELLTYPE in celltypes) {
    print(CELLTYPE)
    de_df  <- rbind(de_df, de_list[[CELLTYPE]])
}
de_df <- tbl_df(de_df)

saveRDS(de_df, "data/ex_de_df_April2021.rds")


