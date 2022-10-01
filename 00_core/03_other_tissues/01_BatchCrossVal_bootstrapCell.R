# Predict Bootstrap Pseudocells using leave-one-batch-out cross validation
library(tidyverse)
library(glmnet)
library(viridis)
library(ggpubr)


df <- readRDS("data/bootstrap_pseudocell_15_seed42.rds")
df <- df %>% ungroup
df$age <- strtoi(strsplit(df$age, "m"))


# Lognormalize counts, including a pseudocount
lognorm <- function(input) {
    norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
    log1p(norm * 10000)
}

df2 <- df %>% mutate(lognorm = map(pseudocell_all, lognorm)) 
df2$pseudocell_all <- NULL
df2 <- unnest(df2, lognorm)
dim(df2) 
df <- NULL # Clear memory

#==================================================================================================
# Leave One Batch Out 
#==================================================================================================

# Initialize main df
full_df <- c()
Celltypes <- c("endothelial cell", "podocyte", "brush cell",
               "mature NK T cell")

for (CT in Celltypes) {
    print(CT)
    df_CT <- filter(df2, cell_ontology_class == CT)
    batches <- unique(df_CT$mouse.id)
    for (batch in batches) {
        print(batch)
        # Split object
        df_fmo <- filter(df_CT, mouse.id != batch)
        df_ss <- filter(df_CT, mouse.id == batch)
        # Train
        model <- cv.glmnet(x = as.matrix(df_fmo[,-c(1:4)]),
                           y = as.matrix(df_fmo[,2]),
                           type.measure="mae",
                           standardize=F,
                           relax=F,
                           nfolds = 5)
        # Predict
        testPredictions <- predict(model, newx=as.matrix(df_ss[,-c(1:4)]), s="lambda.min")
        df_ss$Predictions <- testPredictions[,1]
        test_df <- df_ss %>% select("mouse.id", "age", "cell_ontology_class",
                                    "tissue", "Predictions")
        test_df <- as.data.frame(test_df)
        colnames(test_df) <- c("Mouse", "Age", "Celltype", "Batch", "Prediction")
        # Update
        full_df <- rbind(full_df, test_df)
    }
}


dim(full_df)
saveRDS(full_df, "data/leaveOneBatchPredictions_bootstrapPseudocell.rds")
full_df <- readRDS("data/leaveOneBatchPredictions_bootstrapPseudocell.rds")

#==================================================================================================
# Viz and stats
#==================================================================================================

d <- as_tibble(full_df)

sumStats <- c()
for (cell in unique(d$Celltype)) {
    print(cell)
    cell_df <- filter(d, Celltype == cell)
    mae <- round(median(abs(cell_df$Age - cell_df$Prediction)), 3)
    meanae <- round(mean(abs(cell_df$Age - cell_df$Prediction)), 3)
    rho <- round(cor(cell_df$Age, cell_df$Prediction, method = "spearman"), 3)
    r2 <- round(cor(cell_df$Age, cell_df$Prediction, method = "pearson") ** 2, 3)
    print(paste0("MedianAbsErr : ", mae))
    print(paste0("MeanAbsErr : ", meanae))
    print(paste0("r2 : ", r2))
    print(paste0("spearman's rho : ", rho))
    print("------------------------- ")
    sumStats <- rbind(sumStats, c(cell, mae, meanae, rho, r2))
}

sumStats <- as.data.frame(sumStats)
colnames(sumStats) <- c("Celltype", "MedianAbsErr", "MeanAbsErr", "rho", "r2")
sumStats


#==================================================================================================
# Stats on the median values

d2 <- d %>% group_by(Celltype, Age) %>% mutate(med = median(Prediction)) %>% select(-Prediction) %>% distinct()
sumStats <- c()
for (cell in unique(d2$Celltype)) {
    print(cell)
    cell_df <- filter(d2, Celltype == cell)
    mae <- round(median(abs(cell_df$Age - cell_df$med)), 3)
    meanae <- round(mean(abs(cell_df$Age - cell_df$med)), 3)
    rho <- round(cor(cell_df$Age, cell_df$med, method = "spearman"), 3)
    r2 <- round(cor(cell_df$Age, cell_df$med, method = "pearson") ** 2, 3)
    print(paste0("MedianAbsErr : ", mae))
    print(paste0("MeanAbsErr : ", meanae))
    print(paste0("r2 : ", r2))
    print(paste0("spearman's rho : ", rho))
    print("------------------------- ")
    sumStats <- rbind(sumStats, c(cell, mae, meanae, rho, r2))
}
sumStats <- as.data.frame(sumStats)
colnames(sumStats) <- c("Celltype", "MedianAbsErr", "MeanAbsErr", "rho", "r2")
sumStats


celltype_order <- c("endothelial cell", "mature NK T cell", "podocyte")#, "brush cell")
d <- subset(d, Celltype %in% celltype_order)
d2 <- subset(d2, Celltype %in% celltype_order)
d$Celltype <- factor(d$Celltype, levels = celltype_order)
d2$Celltype <- factor(d2$Celltype, levels = celltype_order)

ggplot(d, aes(x = Age, y = Prediction)) +
    geom_bin2d(binwidth = 1.5) +
    scale_fill_distiller(palette="PuBuGn", direction = 1, limits = c(0, 70), oob = scales::squish) +
    facet_wrap(Celltype~.) +
    theme_classic() +
    stat_cor(label.y = 28, size = 2.5, color = "darkblue") +
    geom_point(data = d2, aes(x=Age, y=med), size = 1.5, color = "black") +
    geom_smooth(data = d2, aes(x=Age, y=med), method = "lm", color = "black", size = 1) +
    stat_cor(data = d2, aes(x=Age, y=med), label.y = 31, size = 2.7) # pearson r not rho
ggsave("plots/tabula_muris_senis_CV_predictions.pdf", width = 6.36, height = 3.29, useDingbats = F)










