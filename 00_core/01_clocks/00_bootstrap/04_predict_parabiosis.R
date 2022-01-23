

library(tidyverse)
library(glmnet)
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/00_bootstrap")


models <- readRDS("data/models_all_bootstrap.rds")
colnames(models)[1] <- "Celltype"
models$lognormalized <- NULL

#==================================================================================================
# Load parabiosis bootstrapped cells

parabiosis_data <- readRDS("data/bootstrap_pseudocell_15_parabiosis.rds")

meta <- as.data.frame(parabiosis_data[, c(1:3)])
umi <- parabiosis_data[, -c(1:3)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)

by_celltype <- as_tibble(cbind(meta, logged)) %>%
            group_by(Celltype) %>%
            nest()

by_celltype <-  dplyr::inner_join(by_celltype, models, by = "Celltype")

by_celltype[1,2][[1]][[1]][1:5,1:4]

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

range(output$Pred) # -1.435905 32.640380
head(output)
#   Celltype        Pred Sample    Mouse                                  
#   <chr>          <dbl> <chr>     <chr>                                  
# 1 Astrocyte_qNSC  8.10 Young-Het BC2-Young-Het-5.10-1B-B1-14300-TCAATGGC
# 2 Astrocyte_qNSC 10.7  Young-Het BC2-Young-Het-5.10-1B-B1-14300-TCAATGGC
# 3 Astrocyte_qNSC 10.5  Young-Het BC2-Young-Het-5.10-1B-B1-14300-TCAATGGC
# 4 Astrocyte_qNSC 12.0  Young-Het BC2-Young-Het-5.10-1B-B1-14300-TCAATGGC
# 5 Astrocyte_qNSC 11.0  Young-Het BC2-Young-Het-5.10-1B-B1-14300-TCAATGGC
# 6 Astrocyte_qNSC  9.84 Young-Het BC2-Young-Het-5.10-1B-B1-14300-TCAATGGC

saveRDS(output, "data/parabiosis_predictions.rds")
