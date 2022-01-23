# Train with full set of data. Tidy form.

library(tidyverse)
library(glmnet)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/01_bootstrap_bioage")
df <- readRDS("data/bootstrap_cell_15_bio.rds")
df <- df %>% ungroup %>%
            select(-c(orig.ident))
head(df)

# # A tibble: 6 x 4
#   hash.ID          Prolif_Lineage_Fraction_of_SVZ Celltype.LowRes pseudocell_all
#   <chr>                                     <dbl> <fct>           <list>        
# 1 Sample4-ACCAATGC                         0.0587 Oligodendro     <tibble [100 ~
# 2 Sample4-ACCAATGC                         0.0587 Endothelial     <tibble [100 ~
# 3 Sample4-ACCAATGC                         0.0587 Microglia       <tibble [100 ~
# 4 Sample1-TGTGATGG                         0.208  Oligodendro     <tibble [100 ~
# 5 Sample3-CTCTAGAC                         0.0710 Microglia       <tibble [100 ~
# 6 Sample2-TCAATGGC                         0.157  Oligodendro     <tibble [100 ~


lognorm <- function(input) {
    norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
    log1p(norm * 10000)
}

df2 <- df %>% mutate(lognorm = map(pseudocell_all, lognorm))
df2$pseudocell_all <- NULL
df2 <- unnest(df2, lognorm)
df <- NULL # Clear memory


by_celltype <- df2 %>% group_by(Celltype.LowRes) %>% nest()
colnames(by_celltype)[2] <- "lognormalized"

# Fit the model, add predictions
# Model function
celltype_model <- function(input) {
    cv.glmnet(x = as.matrix(input[, -c(1,2)]), y = as.matrix(input[, 2]),
              type.measure = "mae", standardize = F, relax = F, nfolds = 5)
}

# Apply model function using map()
models <- by_celltype %>% mutate(model = map(lognormalized, celltype_model))
# Don't overwrite
#saveRDS(models, "data/models_all_bootstrapCell_bio.rds")
