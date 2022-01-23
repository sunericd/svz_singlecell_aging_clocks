# Visualizing how clock genes are impacted by interventions.
# Specifically, what is the overlap in "rejuvenated" genes like?
# (Formerly Script 14.5 in vizclocks)

#================================================================================================
# PART ONE: Set up and load data
library(tidyverse)
library(ggthemes)
library(ComplexHeatmap)
library(ggvenn)
setwd("~/Dropbox/svz_singlecell_aging_clocks/00_core/00_multiseq_aging/Integrate")


# Clock Genes
clock_genes <- readRDS("../../01_clocks/00_bootstrap/data/clock_genes_all.rds")
colnames(clock_genes)[3] <- "celltype"
head(clock_genes)
#          gene     coef  Celltype
# 24      Crlf2 2.855687 Microglia
# 58    Naalad2 2.595200 Microglia
# 33       Klk8 2.204446 Microglia

# Exercise DEG # Generated in 02_Run_DE_Pb_Ex_Multi.R
ex <- readRDS("data/ex_de_df_April2021.rds")
ex <- ex %>% select(gene, celltype, fdr, avg_logFC)
colnames(ex)[3:4] <- c("fdr_ex", "avglogFC_ex")


pb <- readRDS("data/pb_both_de_df_Sep2021.rds") # ERIC: using combined
pb <- pb %>% select(gene, celltype, fdr, avg_logFC)
colnames(pb)[3:4] <- c("fdr_pb", "avglogFC_pb")

#================================================================================================
# PART TWO:

# Subset to clock genes by merging with clock_genes columns
combo <- merge(x = clock_genes, y = ex, by = c("celltype", "gene"))
combo <- merge(x = combo, y = pb, by = c("celltype", "gene"))


#================================================================================================
# Positive clock coeff + positive exercise DEG FC + positive pb DEG FC = 
# Up with age, high in sedentary than exercised, high in old iso than in old heterochronic
dim(filter(combo, coef > 0))[1]
dim(filter(combo, coef > 0, avglogFC_ex > 0, avglogFC_pb > 0))[1]
dim(filter(combo, coef > 0, avglogFC_ex < 0, avglogFC_pb > 0))[1]
dim(filter(combo, coef > 0, avglogFC_ex > 0, avglogFC_pb < 0))[1]
dim(filter(combo, coef > 0, avglogFC_ex < 0, avglogFC_pb < 0))[1]

#================================================================================================
# Effects on genes that go up with age

combo_df <- c()
for (ct in unique(combo$celltype)) {
	print(ct)
	combo_ct <- filter(combo, celltype == ct)
	all_up <- dim(filter(combo_ct, coef > 0))[1] # 673
	both_rev <- dim(filter(combo_ct, coef > 0, avglogFC_ex > 0, avglogFC_pb > 0))[1]
	pb_rev_only <- dim(filter(combo_ct, coef > 0, avglogFC_ex < 0, avglogFC_pb > 0))[1]
	ex_rev_only <- dim(filter(combo_ct, coef > 0, avglogFC_ex > 0, avglogFC_pb < 0))[1]
	neither_rev <- dim(filter(combo_ct, coef > 0, avglogFC_ex < 0, avglogFC_pb < 0))[1]
	add_row <- c(ct, all_up, both_rev, pb_rev_only, ex_rev_only, neither_rev)
	combo_df <- rbind(combo_df, add_row)
}

colnames(combo_df) <- c("Celltype", "all_up", "both_rev", "pb_rev_only", "ex_rev_only", "neither_rev")
rownames(combo_df) <- NULL
combo_df <- as.data.frame(combo_df)
combo_df

#         Celltype all_up both_rev pb_rev_only ex_rev_only neither_rev
# 1 Astrocyte_qNSC    211       57          86          32          36
# 2    Endothelial    100       19          45          15          21
# 3      Microglia     46       13          16           7          10
# 4     Neuroblast    172       67          37          35          33
# 5    Oligodendro     74       18          18          23          15
# 6       aNSC_NPC     75       24          20          17          14  <<<

#================================================================================================
# Stacked barplots with polar coords = Pie plots
combo_df$all_up <- NULL
combo_df <- pivot_longer(combo_df, cols = c(2:5), names_to = "group")
combo_df$value <- as.numeric(as.character(combo_df$value))
CELLTYPES <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")
combo_df$Celltype <- factor(combo_df$Celltype, levels = CELLTYPES)
combo_df$group <- factor(combo_df$group, levels = rev(c("both_rev", "pb_rev_only", "ex_rev_only", "neither_rev")))


p <- ggplot(data = combo_df, aes(x = factor(1), y = value, fill = group)) +
	geom_bar(stat = "identity", position=position_fill()) +
	facet_wrap(.~Celltype) +
	#coord_flip() +
	coord_polar(theta = "y") +
	theme_void() +
	scale_fill_manual(values = rev(c("#023047", "#219ebc", "#8ecae6", "#ffd7ba")))
p
ggsave("plots/Pie_Effect_Clock_Genes_Age_Increase_pb.pdf", width = 4.95, height = 2.31)


#================================================================================================
# Effects on genes that go DOWN with age

combo_df <- c()
for (ct in unique(combo$celltype)) {
	print(ct)
	combo_ct <- filter(combo, celltype == ct)
	all_up <- dim(filter(combo_ct, coef < 0))[1]
	both_rev <- dim(filter(combo_ct, coef < 0, avglogFC_ex < 0, avglogFC_pb < 0))[1]
	pb_rev_only <- dim(filter(combo_ct, coef < 0, avglogFC_ex > 0, avglogFC_pb < 0))[1]
	ex_rev_only <- dim(filter(combo_ct, coef < 0, avglogFC_ex < 0, avglogFC_pb > 0))[1]
	neither_rev <- dim(filter(combo_ct, coef < 0, avglogFC_ex > 0, avglogFC_pb > 0))[1]
	add_row <- c(ct, all_up, both_rev, pb_rev_only, ex_rev_only, neither_rev)
	combo_df <- rbind(combo_df, add_row)
}

colnames(combo_df) <- c("Celltype", "all_down", "both_rev", "pb_rev_only", "ex_rev_only", "neither_rev")
rownames(combo_df) <- NULL
combo_df <- as.data.frame(combo_df)
combo_df

#         Celltype all_down both_rev pb_rev_only ex_rev_only neither_rev
# 1 Astrocyte_qNSC      148       25          16          69          38
# 2    Endothelial      113       31          25          44          13
# 3      Microglia       50        7           6          24          13
# 4     Neuroblast      169       54          32          42          41
# 5    Oligodendro       59       18           8          23          10
# 6       aNSC_NPC       88       34          15          18          21

#================================================================================================
# Stacked barplots
combo_df$all_down <- NULL
combo_df <- pivot_longer(combo_df, cols = c(2:5), names_to = "group")
combo_df$value <- as.numeric(as.character(combo_df$value))
CELLTYPES <- c("Oligodendro", "Microglia", "Endothelial", "Astrocyte_qNSC", "aNSC_NPC", "Neuroblast")
combo_df$Celltype <- factor(combo_df$Celltype, levels = CELLTYPES)
combo_df$group <- factor(combo_df$group, levels = rev(c("both_rev", "pb_rev_only", "ex_rev_only", "neither_rev")))



p <- ggplot(data = combo_df, aes(x = factor(1), y = value, fill = group)) +
	geom_bar(stat = "identity", position=position_fill()) +
	facet_wrap(.~Celltype) +
	#coord_flip() +
	coord_polar(theta = "y") +
	theme_void() +
	scale_fill_manual(values = rev(c("#023047", "#219ebc", "#8ecae6", "#ffd7ba")))
p
ggsave("plots/Pie_Effect_Clock_Genes/Age_Decrease_pb.pdf", width = 4.95, height = 2.31)




#=======================================================
# ERIC 8/25/2021: aNSC-only look at genes with intervention effect and clock usage
combo_aNSC <- filter(combo, celltype == "aNSC_NPC")

# get individual gene lists for enrichment analysis
both_rev_aNSC <- filter(combo_aNSC, coef*avglogFC_ex > 0, coef*avglogFC_pb > 0)
write.csv(both_rev_aNSC$gene,"./data/pie_chart_both_aNSC.csv", row.names=F, quote=F)
pb_rev_only_aNSC <- filter(combo_aNSC, coef*avglogFC_ex < 0, coef*avglogFC_pb > 0)
write.csv(pb_rev_only_aNSC$gene,"./data/pie_chart_pb_only_aNSC.csv", row.names=F, quote=F)
ex_rev_only_aNSC <- filter(combo_aNSC, coef*avglogFC_ex > 0, coef*avglogFC_pb < 0)
write.csv(ex_rev_only_aNSC$gene,"./data/pie_chart_ex_only_aNSC.csv", row.names=F, quote=F)
neither_rev_aNSC <- filter(combo_aNSC, coef*avglogFC_ex < 0, coef*avglogFC_pb < 0)
write.csv(neither_rev_aNSC$gene,"./data/pie_chart_neither_aNSC.csv", row.names=F, quote=F)


# get two sets of genes for upset plot
pb_lt = filter(combo_aNSC, coef*avglogFC_pb > 0)
ex_lt = filter(combo_aNSC, coef*avglogFC_ex > 0)

# make upset plot
pb_genes <- pb_lt$gene
ex_genes <- ex_lt$gene
lt <- list(parabiosis = pb_genes, exercise = ex_genes)
m = make_comb_mat(lt)

pdf("./plots/Complex_Upset.pdf", width = 3, height = 3)
ht = draw(UpSet(m))
od = column_order(ht)
cs = comb_size(m)

decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), gp = gpar(fontsize = 8), rot = 45)
})
dev.off()

#========================================
# ERIC 8/30/2021: aNSC-only look at genes with intervention effect and clock usage

# pb and ex have DGE for all genes 
pb_aNSC <- filter(pb, celltype=="aNSC_NPC")
ex_aNSC <- filter(ex, celltype=="aNSC_NPC")

# get gene lists
fdr_cut <- 0.1
fold_cut <- log(1.1)
up_clock_ex <- filter(combo_aNSC, coef<0, avglogFC_ex<0)$gene # avglogFC defined w.r.t. to treatment so <0 is higher in treatment and >0 is less in treatment
down_clock_ex <- filter(combo_aNSC, coef>0, avglogFC_ex>0)$gene
up_clock_pb <- filter(combo_aNSC, coef<0, avglogFC_pb<0)$gene
down_clock_pb <- filter(combo_aNSC, coef>0, avglogFC_pb>0)$gene
up_DEG_ex <- filter(ex_aNSC, fdr_ex<fdr_cut, avglogFC_ex< -fold_cut)$gene
down_DEG_ex <- filter(ex_aNSC, fdr_ex<fdr_cut, avglogFC_ex>fold_cut)$gene
up_DEG_pb <- filter(pb_aNSC, fdr_pb<fdr_cut, avglogFC_pb< -fold_cut)$gene
down_DEG_pb <- filter(pb_aNSC, fdr_pb<fdr_cut, avglogFC_pb>fold_cut)$gene

# make venn diagrams
x <- list(
  up_clock_ex = sapply(up_clock_ex, as.character), 
  up_DEG_ex = up_DEG_ex 
)
ggvenn(x,fill_color = c("#0073C2FF", "#EFC000FF"),
       stroke_size = 0.5, set_name_size = 4)

x <- list(
  down_clock_ex = sapply(down_clock_ex, as.character), 
  down_DEG_ex = down_DEG_ex 
)
ggvenn(x,fill_color = c("#0073C2FF", "#EFC000FF"),
       stroke_size = 0.5, set_name_size = 4)

x <- list(
  up_clock_pb = sapply(up_clock_pb, as.character), 
  up_DEG_pb = up_DEG_pb 
)
ggvenn(x,fill_color = c("#0073C2FF", "#EFC000FF"),
       stroke_size = 0.5, set_name_size = 4)


x <- list(
  down_clock_pb = sapply(down_clock_pb, as.character), 
  down_DEG_pb = down_DEG_pb 
)
ggvenn(x,fill_color = c("#0073C2FF", "#EFC000FF"),
       stroke_size = 0.5, set_name_size = 4)

