
library(tidyverse)
library(glmnet)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/01_bootstrap_bioage")

models <- readRDS("data/models_all_bootstrapCell_bio.rds")
colnames(models)[1] <- "Celltype"
models$lognormalized <- NULL

#==================================================================================================
# Load parabiosis bootstrapCells

parabiosis_data <- readRDS("../00_bootstrap/data/bootstrap_pseudocell_15_parabiosis.rds")

meta <- as.data.frame(parabiosis_data[, c(1:3)])
umi <- parabiosis_data[, -c(1:3)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)

by_celltype <- as_tibble(cbind(meta, logged)) %>%
            group_by(Celltype) %>%
            nest()
by_celltype <-  dplyr::inner_join(by_celltype, models, by = "Celltype")
by_celltype
# # A tibble: 6 x 3
# # Groups:   Celltype [6]
#   Celltype       data                      model     
#   <chr>          <list>                    <list>    
# 1 Astrocyte_qNSC <tibble [2,200 x 20,950]> <cv.glmnt>
# 2 Oligodendro    <tibble [2,200 x 20,950]> <cv.glmnt>
# 3 Microglia      <tibble [2,200 x 20,950]> <cv.glmnt>
# 4 aNSC_NPC       <tibble [2,200 x 20,950]> <cv.glmnt>
# 5 Endothelial    <tibble [2,200 x 20,950]> <cv.glmnt>
# 6 Neuroblast     <tibble [2,200 x 20,950]> <cv.glmnt>

custom_add_predictions <- function(data, model) {
    predictions <- predict(model, newx = as.matrix(data[,-(1:2)]), s = "lambda.min")
    p <- as_tibble(data.frame("Pred" = predictions[,1],
                              "Sample" = data[,1],
                              "Mouse" = data[,2]))
    return(p)
}

# Add predictions 
by_celltype <- by_celltype %>%
             mutate(Predictions = map2(data, model, custom_add_predictions))
output <- by_celltype %>% select(-c(data, model)) %>% unnest(Predictions)
range(output$Pred) # 0.008666655 0.316443725

#==================================================================================================
saveRDS(output, "data/parabiosis_bootstrapCell_predictions_bio.rds")
