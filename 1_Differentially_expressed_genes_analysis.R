library(tidyverse)
library(org.Mm.eg.db)
library(pheatmap)
library(clusterProfiler)
library(PCAtools)

#Identifying the differentially expressed genes (DEGs) in middle cerebral artery occlusion (MCAO) compared with sham animals.

#1- Loading dataset
nano.counts <- read.csv("Nanostring_counts.csv")
row.names(nano.counts) <- nano.counts$Gene_symbol
nano.counts <- nano.counts[,-1]

#2- Heatmap plot
heat.map <- pheatmap(nano.counts,cluster_cols=T, cluster_rows=T, scale="row",show_rownames=F, 
                     show_colnames=T, angle_col="45", treeheight_row=20, border_color=NA, fontsize_row=11, legend=T)

#3- Identification of DEGs
t_test_results <- nano.counts %>%
  rowwise() %>% 
  mutate(p_value = t.test(c_across(starts_with("sham")),c_across(starts_with("mcao")))$p.value) %>%
  ungroup()
DEGs <- t_test_results %>% filter(p_value < 0.05) %>% rownames()

#4- Functional annotation of DEGs
GO.results <- enrichGO(gene = DEGs, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = c("BP","CC"))

#5- Principle Component Analysis (PCA)
pca.results <- pca(nano.counts)



