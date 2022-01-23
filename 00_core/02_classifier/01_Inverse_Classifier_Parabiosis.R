# Train on 75% of Parabiosis Data
# Test on 25% of Parabiosis Data
# Test on Multiseq Timecourse Data

library(tidyverse)
library(glmnet)
library(ggthemes)
library(viridis)
library(ggridges)


#==================================================================================================
# Load Data
#==================================================================================================
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/02_classifier")

parabiosis_data <- readRDS("../01_clocks/00_bootstrap/data/bootstrap_pseudocell_15_parabiosis.rds")

meta <- as.data.frame(parabiosis_data[, c(1:3)])
umi <- parabiosis_data[, -c(1:3)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)

pb  <- as_tibble(cbind(meta, logged)) %>% filter(Sample %in% c("Old-Het", "Old-Iso"))
pb$target <- pb$Sample
pb$target[pb$target == "Old-Het"] <- 0
pb$target[pb$target == "Old-Iso"] <- 1
pb <- pb %>% select(target, everything())
pb$target <- as.numeric(pb$target)


sample <- sample.int(n = nrow(pb), size = floor(.75*nrow(pb)), replace = F)
pb_train <- pb[sample, ]
pb_test  <- pb[-sample, ]

pb_train_bycelltype <- pb_train %>%
                       group_by(Celltype) %>%
                       nest()
pb_test_bycelltype <- pb_test %>%
                       group_by(Celltype) %>%
                       nest()

#==================================================================================================
# Train Logistic regression model with cross validation
#==================================================================================================
celltype_model <- function(input) {

    cv.glmnet(x = as.matrix(input[, -c(1:3)]),
              y = as.matrix(input[, 1]),
              type.measure = "mse",
              family = "binomial",
              standardize = F,
              relax = F, 
              nfolds = 5)
}

models <- pb_train_bycelltype %>% mutate(model = map(data, celltype_model))
models$data <- NULL

#==================================================================================================
# Evaluate on held out parabiosis data (sanity check)
#==================================================================================================
custom_add_predictions <- function(data, model) {
    predictions <- predict(model, newx = as.matrix(data[,-(1:3)]), s = "lambda.min")
    p <- as_tibble(data.frame("Pred" = predictions[,1],
                              "Sample" = data[,2],
                              "Mouse" = data[,3]))
    return(p)
}

# Add predictions
by_celltype <-  dplyr::inner_join(pb_test_bycelltype, models, by = "Celltype")
by_celltype <- by_celltype %>%
             mutate(Predictions = map2(data, model, custom_add_predictions))

output <- by_celltype %>% select(-c(data, model)) %>% unnest(Predictions)

ggplot(filter(output, Celltype == "aNSC_NPC"), aes(x=Pred, fill=Sample)) +
    geom_density(alpha = 0.75) +
    theme_classic() +
    scale_fill_tableau()


#==================================================================================================
# Evaluate on Timecourse
#==================================================================================================

df <- readRDS("../01_clocks/00_bootstrap/data/bootstrap_pseudocell_15.rds")
df <- df %>% ungroup %>%
            select(-c(hash.ID, orig.ident, data))

lognorm <- function(input) {
    norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
    log1p(norm * 10000)
}

df2 <- df %>% mutate(lognorm = map(pseudocell_all, lognorm))
df2$pseudocell_all <- NULL
df2 <- unnest(df2, lognorm)
df <- NULL # Clear memory

by_celltype <- df2 %>% group_by(Celltype.LowRes) %>% nest()
colnames(by_celltype) <- c("Celltype", "data")

# Combine with classifiers
by_celltype <-  dplyr::inner_join(by_celltype, models, by = "Celltype")

# Predict
custom_add_predictions2 <- function(data, model) {
    predictions <- predict(model, newx = as.matrix(data[,-1]), s = "lambda.min")
    p <- as_tibble(data.frame("Pred" = predictions[,1],
                              "Age" = data[,1]))
    return(p)
}

by_celltype <- by_celltype %>%
             mutate(Predictions = map2(data, model, custom_add_predictions2))

output <- by_celltype %>% select(-c(data, model)) %>% unnest(Predictions)
saveRDS(output, "data/parabiosis_classifier_on_multi.rds")

