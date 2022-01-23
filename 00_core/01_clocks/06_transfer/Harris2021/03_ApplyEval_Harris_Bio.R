library(Seurat)
library(tidyverse)
library(glmnet)
library(ggthemes)
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/06_transfer/Harris2021")

# Load all-mice-trained cell type specific bootstrap models and data
celltype_models <- readRDS("../../01_bootstrap_bioage/data/models_all_bootstrapCell_bio.rds")
genes <- colnames(celltype_models[1,2][[1]][[1]])[-1]

# Load Harris Pseudocell dataframe
ham <- readRDS("data/bootstrap_pseudocell_15_harris.rds")
ham.genes <- colnames(ham[1,4][[1]][[1]])

lognorm <- function(input) {
    norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
    log1p(norm * 10000)
}
df2 <- ham %>% mutate(lognormed = map(pseudocell_all, lognorm))

# Reorder and resize genes
regene <- function(input) {
    missing <- setdiff(genes, ham.genes)
    missing_df <- matrix(0, 100, length(missing))
    colnames(missing_df) <- missing
    rownames(missing_df) <- rownames(input)
    output <- cbind(input, missing_df)
    output[, genes] # Reorder and cut to size to match clock training data format
}

ham <- df2 %>% mutate(lognorm2 = map(lognormed, regene))
ham$pseudocell_all <- NULL
ham$lognormed <- NULL


ham <- ham %>% unnest()
ham <- ham %>% select(Celltype, orig.ident, Age, everything())
data <- ham

prediction_df <- c()
for (celltype in unique(celltype_models$Celltype.LowRes)) {
    print(celltype)
    model_celltype <- filter(celltype_models, Celltype.LowRes == celltype) # Subset model dataframe
    model_celltype <- model_celltype[1,3][[1]][[1]] # Extract model
    print(model_celltype)
    prediction <- predict(model_celltype, newx=as.matrix(data[,-c(1:4)]), s="lambda.min")
    output_df <- data.frame("Model" = rep(celltype, length(prediction)),
                            "Predicted_Celltype" = data$Celltype,
                            "ID" = data$orig.ident,
                            "minus_Predicted_Fraction" = 100*(0.35-prediction[,1]),
                            "Age" = data$Age)
    prediction_df <- rbind(prediction_df, output_df)
}

d <- prediction_df
celltypes <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Niche")
d$Model <- factor(d$Model, levels = celltypes)
d$Age <- factor(d$Age)

ggplot(d, aes(x= minus_Predicted_Fraction, color = Age, fill = Age)) +
    geom_density(alpha = 0.7) +
    facet_grid(Model ~ Predicted_Celltype, scales = "free_y") +
    scale_fill_tableau() + scale_color_tableau() +
    #xlim(c(-3,33)) +
    theme_classic() + theme(panel.background = element_rect(colour = "black", size=.6)) +
    xlab("Biological Age Score: 35 - (Proliferative %)") +
    ggtitle("Harris et al (2021) Dentate Gyrus Biological Age Prediction")

ggsave("plots/density_grid_harris_bio.pdf", useDingbats=F)


# Filter and plot matching only

d_filter <- filter(d, as.character(Model) == as.character(Predicted_Celltype))
ggplot(d_filter, aes(x= minus_Predicted_Fraction, color = Age, fill = Age)) +
    geom_density(alpha = 0.7) +
    facet_wrap(.~Model) +
    scale_fill_tableau() + scale_color_tableau() +
    xlim(c(0,30)) +
    theme_classic() + theme(panel.background = element_rect(colour = "black", size=.4)) +
    #theme_bw() + theme(panel.grid.minor = element_blank()) +
    xlab("Biological Age Score: 35 - (Proliferative %)") +
    ylab("") +
    ggtitle("Harris et al (2021) Dentate Gyrus Biological Age Prediction")
ggsave("plots/density_grid_matching_harris_bio.pdf", width = 5.59, height = 3,  useDingbats=F)


