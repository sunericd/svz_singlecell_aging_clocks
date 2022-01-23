library(tidyverse)
library(glmnet)
library(viridis)
library(RColorBrewer)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/05_visualize")

# Load all mice trained cell type specific bootstrap models and data
celltype_models <- readRDS("../00_bootstrap/data/models_70pTrain.rds")
genes <- colnames(celltype_models[1,2][[1]][[1]])[-1]

test_data <- readRDS("../00_bootstrap/data/30pTest.rds")
data <- test_data %>% select(Celltype.LowRes, lognormalized) %>% unnest()

prediction_df <- c()
for (celltype in unique(data$Celltype.LowRes)) {
    print(celltype)
    model_celltype <- filter(celltype_models, Celltype.LowRes == celltype) # Subset model dataframe
    model_celltype <- model_celltype[1,3][[1]][[1]] # Extract model
    print(model_celltype)
    prediction <- predict(model_celltype, newx=as.matrix(data[,-c(1:2)]), s="lambda.min")
    output_df <- data.frame("Model" = rep(celltype, length(prediction)),
                            "Predicted_Celltype" = data$Celltype.LowRes,
                            "Predicted_Age" = prediction[,1],
                            "Age" = data$Age)
    prediction_df <- rbind(prediction_df, output_df)
}


prediction_df2 <- prediction_df %>% group_by(Model, Predicted_Celltype, Age) %>%
    summarise(Median_Prediction = median(Predicted_Age))
                    
prediction_df3 <- prediction_df2 %>% summarise(r = cor(Age, Median_Prediction), MAE = median(abs(Age - Median_Prediction)))

d <- prediction_df3
celltypes <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")
d$Predicted_Celltype <- factor(d$Predicted_Celltype, levels = rev(celltypes))
d$Model <- factor(d$Model, levels = celltypes)


ggplot(data = d, mapping = aes(y = Predicted_Celltype, x = Model, size = r, color = MAE)) +
    geom_point() +
    scale_color_viridis(option = "viridis", begin = 0, end = 1, direction = -1) +
    geom_point(shape = 1, color = "black") +
    scale_size(breaks = c(.3,.5,.7,.9)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30, hjust=1)) +
    theme(legend.position = "right") +
    ylab("Tested On") + xlab("Trained On")
ggsave(paste0("plots/CrossCelltypePerformance_Dotplot_", Sys.Date(), ".pdf"), useDingbats=F, width=3.56, height=2.29)



