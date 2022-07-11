library(tidyverse)
library(glmnet)
library(viridis)
library(RColorBrewer)
library(ggthemes)
library(cowplot)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/06_transfer")

# Load all mice trained cell type specific bootstrap models and data
celltype_models <- readRDS("../../00_bootstrap/data/models_all_bootstrap.rds")
genes <- colnames(celltype_models[1,2][[1]][[1]])[-1]

# Load DBN Pseudocell dataframe
dbn <- readRDS("data/bootstrap_pseudocell_15.rds")
dbn.genes <- colnames(dbn[1,4][[1]][[1]])

lognorm <- function(input) {
    norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
    log1p(norm * 10000)
}
df2 <- dbn %>% mutate(lognormed = map(pseudocell_all, lognorm))

# Reorder and resize genes
regene <- function(input) {
    missing <- setdiff(genes, dbn.genes)
    missing_df <- matrix(0, 100, length(missing))
    colnames(missing_df) <- missing
    rownames(missing_df) <- rownames(input)
    output <- cbind(input, missing_df)
    output[, genes] # Reorder and cut to size to match clock training data format
}

dbn <- df2 %>% mutate(lognorm2 = map(lognormed, regene))
dbn$pseudocell_all <- NULL
dbn$lognormed <- NULL


# Fix things
unique(dbn$Celltype)
dbn$Celltype <- as.character(dbn$Celltype)
dbn$Celltype <- plyr::mapvalues(dbn$Celltype,
                    from = c("Astrocytes_qNSCs", "aNSCs_NPCs", "Neuroblasts","Oligodendrocytes"),
                    to = c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast","Oligodendro"))
dbn$Age <- as.character(dbn$Age)
dbn$Age <- plyr::mapvalues(dbn$Age,
                            from = c("y", "o"),
                            to = c("3.5", "28"))
dbn$Age <- as.numeric(dbn$Age)                  

dbn <- dbn %>% unnest()
dbn <- dbn %>% select(Celltype, orig.ident, Age, everything())
data <- dbn

prediction_df <- c()
for (celltype in unique(data$Celltype)) {
    print(celltype)
    model_celltype <- filter(celltype_models, Celltype.LowRes == celltype) # Subset model dataframe
    model_celltype <- model_celltype[1,3][[1]][[1]] # Extract model
    print(model_celltype)
    prediction <- predict(model_celltype, newx=as.matrix(data[,-c(1:3)]), s="lambda.min")
    output_df <- data.frame("Model" = rep(celltype, length(prediction)),
                            "Predicted_Celltype" = data$Celltype,
                            "ID" = data$orig.ident,
                            "Predicted_Age" = prediction[,1],
                            "Age" = data$Age)
    prediction_df <- rbind(prediction_df, output_df)
}

d <- prediction_df
celltypes <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Niche")
d$Predicted_Celltype <- factor(d$Predicted_Celltype, levels = celltypes)
d$Model <- factor(d$Model, levels = celltypes)
d$Age <- factor(d$Age)

# Plot full matrix
ggplot(d, aes(x= Predicted_Age, color = Age, fill = Age)) +
    geom_density(alpha = 0.7) +
    facet_grid(Model ~ Predicted_Celltype, scales = "free_y") +
    scale_fill_tableau() + scale_color_tableau() +
    xlim(c(-3,33)) +
    theme_classic() + theme(panel.background = element_rect(colour = "black", size=.4)) +
    #theme_bw() + theme(panel.grid.minor = element_blank()) +
    ylab("Tested On Dulken2019") + xlab("Trained On")

ggsave("plots/density_grid_all.pdf", useDingbats=F)


# Filter and plot matching only
d_filter <- filter(d, Model == Predicted_Celltype)
ggplot(d_filter, aes(x= Predicted_Age, color = Age, fill = Age)) +
    geom_density(alpha = 0.7) +
    facet_wrap(.~Model) +
    scale_fill_tableau() + scale_color_tableau() +
    xlim(c(-3,33)) +
    theme_classic() + theme(panel.background = element_rect(colour = "black", size=.4)) +
    #theme_bw() + theme(panel.grid.minor = element_blank()) +
    ylab("Tested On Dulken2019") + xlab("Trained On")

ggsave("plots/density_grid_matching.pdf", width = 4.59, height = 2.42,  useDingbats=F)






