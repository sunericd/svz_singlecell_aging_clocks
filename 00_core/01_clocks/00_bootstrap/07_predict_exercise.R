library(tidyverse)
library(glmnet)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")

models <- readRDS("data/models_all_bootstrap.rds")
colnames(models)[1] <- "Celltype"
models$lognormalized <- NULL

#==================================================================================================
# Load exercise bootstrapCells

ex_data <- readRDS("data/bootstrap_pseudocell_15_exercise.rds")

meta <- as.data.frame(ex_data[, c(1:4)])
umi <- ex_data[, -c(1:4)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)

by_celltype <- as_tibble(cbind(meta, logged)) %>%
            group_by(Celltype) %>%
            nest()

by_celltype <-  dplyr::inner_join(by_celltype, models, by = "Celltype")

by_celltype[1,2][[1]][[1]][1:5,1:4]

custom_add_predictions <- function(data, model) {
    predictions <- predict(model, newx = as.matrix(data[,-(1:3)]), s = "lambda.min")
    p <- as_tibble(data.frame("Pred" = predictions[,1],
                              "Sample" = data[,1],
                              "Mouse" = data[,2],
                              "Year" = data[,3]))
    return(p)
}

# Add predictions 
by_celltype <- by_celltype %>%
             mutate(Predictions = map2(data, model, custom_add_predictions))


output <- by_celltype %>% select(-c(data, model)) %>% unnest(Predictions)

range(output$Pred)

saveRDS(output, "data/exercise_predictions.rds")
