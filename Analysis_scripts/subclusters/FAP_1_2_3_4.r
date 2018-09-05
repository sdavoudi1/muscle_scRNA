# FAP_1_2_3_4
#

# This script compares the 4 FAP subclusters in the young population.
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/subclusters/FAP_1_2_3_4.r")

# ------------------------------------------------------------------------------------------

# First we load the required functions and libraries.

if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
library(Seurat)

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization_functions.r")

# ------------------------------------------------------------------------------------------

# Load the datasets
young_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune.rds")

# ------------------------------------------------------------------------------------------

# We determine the differentially expressed genes in each subcluster.
FAP1_genes <- FindMarkers(object = young_noimmune, ident.1 = 0, ident.2 = c(3,5,7), min.pct = 0.25, only.pos = T)
FAP2_genes <- FindMarkers(object = young_noimmune, ident.1 = 3, ident.2 = c(0,5,7), min.pct = 0.25, only.pos = T)
FAP3_genes <- FindMarkers(object = young_noimmune, ident.1 = 5, ident.2 = c(0,3,7), min.pct = 0.25, only.pos = T)
FAP4_genes <- FindMarkers(object = young_noimmune, ident.1 = 7, ident.2 = c(0,3,5), min.pct = 0.25, only.pos = T)

# Next we visualize the difference in the top 15 genes for the two clusters

gene_list = c(rownames(FAP1_genes)[1:15], rownames(FAP2_genes)[1:15], rownames(FAP3_genes)[1:15], rownames(FAP4_genes)[1:15])

circle_heatmap_genelist_spec_cluster(dataset = young_noimmune, gene_list, n_genes = 15,
								clusters = c(1,4,6,8), rename_cluster = T, cluster_names = c("FAP_1", "FAP_2", "FAP_3", "FAP_4"),
								image_dpi = 600, image_width = 20, image_height = 5,chart_name = "",
								image_name = "FAP_comparisons (top 15 DEG).png")
								
# ------------------------------------------------------------------------------------------

# Next we look at the differentially expressed genes in each subcluster that has less than 25% expression
# in the other one.

pct_exp_FAP1 <- pct_exp_per_cluster(dataset = young_noimmune, gene_list = rownames(FAP1_genes))
pct_exp_FAP1 <- pct_exp_FAP1[,c(1,4,6,8)]
FAP1_genes_unique <- pct_exp_FAP1[intersect(intersect(which(pct_exp_FAP1[,2] < 0.25), which(pct_exp_FAP1[,3] < 0.25)), which(pct_exp_FAP1[,4] < 0.25)),]

pct_exp_FAP2 <- pct_exp_per_cluster(dataset = young_noimmune, gene_list = rownames(FAP2_genes))
pct_exp_FAP2 <- pct_exp_FAP2[,c(1,4,6,8)]
FAP2_genes_unique <- pct_exp_FAP2[intersect(intersect(which(pct_exp_FAP2[,1] < 0.25), which(pct_exp_FAP2[,3] < 0.25)), which(pct_exp_FAP2[,4] < 0.25)),]

pct_exp_FAP3 <- pct_exp_per_cluster(dataset = young_noimmune, gene_list = rownames(FAP3_genes))
pct_exp_FAP3 <- pct_exp_FAP3[,c(1,4,6,8)]
FAP3_genes_unique <- pct_exp_FAP3[intersect(intersect(which(pct_exp_FAP3[,1] < 0.33), which(pct_exp_FAP3[,2] < 0.33)), which(pct_exp_FAP3[,4] < 0.33)),]
# FAP3 didn't have 15 unique differentially expressed genes. Had to raise the bar to 33% expression in other clusters.
# FAP3_genes_unique <- pct_exp_FAP3[intersect(intersect(which(pct_exp_FAP3[,1] < 0.25), which(pct_exp_FAP3[,2] < 0.25)), which(pct_exp_FAP3[,4] < 0.25)),]

pct_exp_FAP4 <- pct_exp_per_cluster(dataset = young_noimmune, gene_list = rownames(FAP4_genes))
pct_exp_FAP4 <- pct_exp_FAP4[,c(1,4,6,8)]
FAP4_genes_unique <- pct_exp_FAP4[intersect(intersect(which(pct_exp_FAP4[,1] < 0.25), which(pct_exp_FAP4[,2] < 0.25)), which(pct_exp_FAP4[,3] < 0.25)),]

gene_list = c(rownames(FAP1_genes_unique)[1:10], rownames(FAP2_genes_unique)[1:10], rownames(FAP3_genes_unique)[1:10], rownames(FAP4_genes_unique)[1:10])
circle_heatmap_genelist_spec_cluster(dataset = young_noimmune, gene_list, n_genes = 10,
								clusters = c(1,4,6,8), rename_cluster = T, cluster_names = c("FAP_1", "FAP_2", "FAP_3", "FAP_4"),
								image_dpi = 600, image_width = 20, image_height = 5,chart_name = "",
								image_name = "FAP_comparisons (unique).png")