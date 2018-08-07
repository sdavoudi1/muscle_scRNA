# Young_noimmune_schwann_diff_genes
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/Other/Young_noimmune_schwann_diff_genes.r")

# This file contains the scripts I used to prepare the files I sent to Alaura for the Schwann cell markers.
# Results were saved here: C:\Users\sadeg\Google Drive\scRNA\Results\Schwann Cells Alaura

# -----------------------------------------------------------------------------------------------------------

if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
library(Seurat)

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")

# -----------------------------------------------------------------------------------------------------------

young_noimmune_diff_markers <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_diff_markers.rds")
young_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune.rds")

# -----------------------------------------------------------------------------------------------------------

# This section of the script uses the differentially expressed genes in the Schwann cell cluster,
# and determines the normalized and raw expressions of those genes in all clusters, as well as the 
# percentage of cells in each cluster, expressing those genes.

# Next we find the differential genes expressed in the Schwann cell cluster
Schwann_noimmune_diff_genes <- young_noimmune_diff_markers$gene[which(young_noimmune_diff_markers$cluster == 8)]

# Next we find the average expression of the genes in all clusters.
young_noimmune_schwanngene_avg_exp <- Average_gene_exp_per_cluster(dataset = young_noimmune,
										gene_list = Schwann_noimmune_diff_genes, reorder_cluster = T, 
										cluster_order = c(9, 2, 5, 3, 1, 4, 6, 8, 7), 
										normalize = F, rename_cluster = T, 
										cluster_names = c("Schwann", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "FAP_4", "Tenocytes"))
write.csv(young_noimmune_schwanngene_avg_exp, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Schwann Cells Alaura/Young_noimmune_Schwann_Average_expression.csv")

# Here we create the normalized expression matrix
young_noimmune_schwanngene_avg_exp_norm <- Average_gene_exp_per_cluster(dataset = young_noimmune,
										gene_list = Schwann_noimmune_diff_genes, reorder_cluster = T, 
										cluster_order = c(9, 2, 5, 3, 1, 4, 6, 8, 7), 
										normalize = T, rename_cluster = T, 
										cluster_names = c("Schwann", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "FAP_4", "Tenocytes"))
write.csv(young_noimmune_schwanngene_avg_exp_norm, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Schwann Cells Alaura/Young_noimmune_Schwann_Average_expression_normalized.csv")

# Next we find the percentagte expression of the genes in all clusters.
young_noimmune_schwanngene_pct_exp <- pct_exp_per_cluster(dataset = young_noimmune, gene_list = Schwann_noimmune_diff_genes,
										reorder_cluster = T, cluster_order = c(9, 2, 5, 3, 1, 4, 6, 8, 7),
										rename_cluster = T, cluster_names = c("Schwann", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "FAP_4", "Tenocytes")) 
write.csv(young_noimmune_schwanngene_pct_exp, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Schwann Cells Alaura/Young_noimmune_Schwann_percent_expression.csv")