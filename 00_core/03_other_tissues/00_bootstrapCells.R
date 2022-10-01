library(Seurat)
library(tidyverse)


# Read in CSV files
svz <- read.csv("data/tabula-muris-senis-droplet-official-raw-obj_DATA.csv", header=F)
meta <- read.csv("data/tabula-muris-senis-droplet-official-raw-obj_META.csv")


Convert_to_Dataframe <- function(svz) {
    meta <- meta[, c("mouse.id", "age", "cell_ontology_class", "tissue")]
    raw_counts <- as.data.frame(svz)
    df <- as_tibble(cbind(meta, raw_counts))
    return(df)
}

df <- Convert_to_Dataframe(svz) %>%
        group_by(cell_ontology_class, age, tissue, mouse.id) %>%
        nest()

head(df)


#===================================================================================================

bootstrap.pseudocells <- function(df, size=15, n=100, replace="dynamic") {
    pseudocells <- c()
    # If dynamic then only sample with replacement if required due to shortage of cells.
    if (replace == "dynamic") {
        if (nrow(df) <= size) {replace <- TRUE} else {replace <- FALSE}
    }
    for (i in c(1:n)) {
        batch <- df[sample(1:nrow(df), size = size, replace = replace), ]
        pseudocells <- rbind(pseudocells, colSums(batch))
    }
    colnames(pseudocells) <- colnames(df)
    return(as_tibble(pseudocells))
}
    
# Apply boostrap.pseudocells using map()
set.seed(42)
df2 <- df %>% mutate(pseudocell_all = map(data, bootstrap.pseudocells)) # ~15 minutes on a laptop

#==================================================================================================
# Remove single cell data; keep just key metadata and pseudocells
df2$data <- NULL
saveRDS(df2, "data/bootstrap_pseudocell_15_seed42.rds")


