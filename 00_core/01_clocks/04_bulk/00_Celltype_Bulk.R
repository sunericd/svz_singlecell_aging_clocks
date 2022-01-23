# Cell type specific models after pseudobulking by sample-celltype

library(tidyverse)
library(Seurat)
library(glmnet)
library(glmnetUtils)
library(cowplot)
library(ggpubr)
set.seed(7)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/04_bulk")
d <- readRDS("../../MULTIseq/Integrated/data/multi_intergrated_seurat_Dec2020.rds")


#============================================================================
# Convert to Tidy form and group by cell type
rna.data <- as_tibble(t(as.matrix((d[["RNA"]]@counts))))

df <- tibble("Celltype"=as.factor(d@meta.data$Celltype.LowRes),
            "Mouse"=as.factor(d@meta.data$hash.ID),
            "Batch"=as.factor(d@meta.data$orig.ident),
             "Age"=as.numeric(d@meta.data$Age),
             rna.data)


Celltypes <- c("Microglia", "Endothelial", "Oligodendro", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

by_celltype <- df %>%
    group_by(Celltype, Mouse, Batch, Age) %>%
    nest() %>%
    filter(Celltype %in% Celltypes)

head(by_celltype)


#============================================================================
# Pseudobulk
pseudobulk_sample_celltype <- function(x) {
                            x2 <- t(as.matrix(colSums(x)))
                            x3 <- log1p((x2 / rowSums(x2)) * 10000)
                            as.tibble(x3)
                        }

# Apply pseudobulk function using map()
by_celltype <- by_celltype %>%
    mutate(pseudobulked = map(data, pseudobulk_sample_celltype))

head(by_celltype)


a <- by_celltype[1,4][[1]][[1]]
b <- t(as.matrix(colSums(a)))
c <- log1p((b / rowSums(b)) * 10000)

range(by_celltype[1,5][[1]][[1]])

#============================================================================
# Fit the model, add predictions

# Initialize main df
full_df <- c()


Mice <- unique(by_celltype$Mouse)
Batches <- unique(by_celltype$Batch)
for (CT in Celltypes) {
    print(CT)
    svz_CT <- filter(by_celltype, Celltype == CT)
    for (batch in Batches) {
        print(batch)

        # split
        svz_train <- filter(svz_CT, Batch != batch)
        print(svz_train)
        svz_predict <- filter(svz_CT, Batch == batch)

        # Train
        model <- cv.glmnet(x = as.matrix(unnest(svz_train[,'pseudobulked']), cols = c(pseudobulked)),
                         y = as.matrix(svz_train[,'Age']),
                         type.measure="mae",
                         standardize=F,
                         relax=F,
                         nfolds = 5)
        print(model)

        # Predict
        testPredictions <- predict(model,
                                  newx=as.matrix(unnest(svz_predict[,'pseudobulked']),
                                                        cols = c(pseudobulked)),
                                  s="lambda.min")
        test_df <- data.frame("Age" = svz_predict$Age, "Mouse" = svz_predict$Mouse,
                              "Batch" = svz_predict$Batch, "Prediction" = testPredictions[,1])
        test_df$Celltype <- CT

        # Update the outer world
        full_df <- rbind(full_df, test_df)
    }
}

#============================================================================
saveRDS(full_df, "data/celltype_specific_bulk_leave1batch_results.rds")

ggplot(data = full_df, aes(x = Age, y = Prediction)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(.~Celltype) +
    theme_classic()
ggsave("plots/celltype-specific.pdf", width=6.76, height=4.31)










