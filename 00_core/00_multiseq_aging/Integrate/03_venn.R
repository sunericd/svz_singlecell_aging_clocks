library(tidyverse)
library(Seurat)
library(ggthemes)
library(ggvenn)

setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Integrate")

m <- readRDS("data/age_de_df_April2021.rds")
ex <- readRDS("data/ex_de_df_April2021.rds")
pb <- readRDS("data/pb_both_de_df_Sep2021.rds")
pb$experiment <- NULL

m$exp <- "Multi"
ex$exp <- "ex"
pb$exp <- "pb"
d <- rbind(m, ex, pb)

CELLTYPES <- c("Astrocyte_qNSC", "Oligodendro", "Microglia",
               "Endothelial", "Neuroblast", "aNSC_NPC")


# Apply cut offs
fdr_cut <- 0.1
fold_cut <- log(1.1)

df_all <- c()
for (cell_type in CELLTYPES) {
    print(cell_type)

    d_filter <- d %>% filter(abs(avg_logFC) > fold_cut, fdr < fdr_cut, celltype == cell_type)

    print(table(d_filter$comparision))

    age_genes <- filter(d_filter, comparision == "Old/Young")$gene
    pb_genes <- filter(d_filter, comparision == "OldIso/OldHet")$gene
    ex_genes <- filter(d_filter, comparision == "OldSed/OldEx")$gene

    all_olap <- length(intersect(age_genes, intersect(pb_genes, ex_genes)))
    age_pb_olap <- length(intersect(age_genes, pb_genes)) - all_olap
    age_ex_olap <- length(intersect(age_genes, ex_genes)) - all_olap
    pb_ex_olap <- length(intersect(ex_genes, pb_genes))
    age_only <- length(age_genes) - age_pb_olap - age_ex_olap - all_olap
    pb_only <- length(pb_genes) - age_pb_olap - pb_ex_olap - all_olap
    ex_only <- length(ex_genes) - age_ex_olap - pb_ex_olap - all_olap

    df_celltype <- c(cell_type, age_only, pb_only, ex_only,
                     age_pb_olap, age_ex_olap, pb_ex_olap,
                     all_olap)

    df_all <- rbind(df_all, df_celltype)

    venn_data <- list(
        "Aging" = age_genes,
        "Parabiosis" = pb_genes,
        "Exercise" = ex_genes
        )
    p <- ggvenn(
      venn_data,
      fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
      stroke_size = 0.5, set_name_size = 4
      )

    ggsave(paste0("plots/venn/", cell_type, "_fdr_", fdr_cut, "_fold_", exp(1)**fold_cut, "_gg_pb1.pdf"), p, width=4, height=4)
}