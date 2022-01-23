library(vip)
library(caret)
library(glmnet)
library(tidyverse)
library(ggthemes)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/01_clocks/05_visualize")
models <- readRDS("../00_bootstrap/data/models_all_bootstrap.rds")
models$lognormalized <- NULL

celltypes = c("Oligodendro", "Microglia", "Endothelial",
            "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")

# Coefficient Donut Plots
for (i in c(1:6)) {
    ct <- as.character(models[i,1][[1]][[1]])
    print(ct)
    lasso <- models[i,2][[1]][[1]]
    print(lasso)
    c1 = vip::vi_model(lasso$glmnet.fit, s = lasso$lambda.min, method = "shap")
    d <- c1 %>% filter(Importance != 0)

    d$ImportanceSigned <- d$Importance
    d$ImportanceSigned[d$Sign == 'NEG'] <- -d$ImportanceSigned[d$Sign == 'NEG']

    # Compute percentages
    d <- d[order(d$ImportanceSigned),]
    d$fraction <- abs(d$Importance) / sum(abs(d$Importance))

    # Compute the cumulative percentages (top of each rectangle)
    d$ymax <- cumsum(d$fraction)

    # Compute the bottom of each rectangle
    d$ymin <- c(0, head(d$ymax, n=-1))

    # Compute label position
    d$labelPosition <- (d$ymax + d$ymin) / 2
    write.csv(d, paste0("data/importance_", ct, ".csv"))
    
    # Do the do not
    p <- ggplot(d, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Sign)) +
      geom_rect(color = "black", size = 0.1) +
      geom_label(x=4.2, aes(y=labelPosition, label=Variable), size=2.5) +
      coord_polar(theta="y") +
      xlim(c(2, 4)) +
      theme_void() +
      theme(legend.position = "none") +
      scale_fill_tableau() +
      annotate(geom = 'text', size = 9, x = 2, y = 0, label = paste(ct, "\n", length(d$Sign)))
    #ggsave(paste0("plots/donut_", ct, ".pdf"), p, width=7.07, height=6.78)
}



