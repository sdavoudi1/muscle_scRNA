# EC1_vs_EC2
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/subclusters/EC1_vs_EC2.r")

# This script compares the EC1 and EC2 clusters to look at the differentially expressed genes.

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

EC1_genes <- FindMarkers(object = young_noimmune, ident.1 = 1, ident.2 = 4, min.pct = 0.25, only.pos = T)
EC2_genes <- FindMarkers(object = young_noimmune, ident.1 = 4, ident.2 = 1, min.pct = 0.25, only.pos = T)

# Next we visualize the difference in the top 15 genes for the two clusters

gene_list = c(rownames(EC1_genes)[1:15], rownames(EC2_genes)[1:15])

circle_heatmap_genelist_spec_cluster(dataset = young_noimmune, gene_list, n_genes = 15,
								clusters = c(2,5), rename_cluster = T, cluster_names = c("EC_1", "EC_2"),
								image_dpi = 600, image_width = 10, image_height = 5,chart_name = "",
								image_name = "EC1_vs_EC2 (top 15 DEG).png")
								
# ------------------------------------------------------------------------------------------

# Next we look at the differentially expressed genes in each subcluster that has less than 25% expression
# in the other one.

pct_exp_EC1 <- pct_exp_per_cluster(dataset = young_noimmune, gene_list = rownames(EC1_genes))
pct_exp_EC1 <- pct_exp_EC1[,c(2,5)]
EC1_genes_unique <- pct_exp_EC1[which(pct_exp_EC1[,2] < 0.25),]

pct_exp_EC2 <- pct_exp_per_cluster(dataset = young_noimmune, gene_list = rownames(EC2_genes))
pct_exp_EC2 <- pct_exp_EC2[,c(2,5)]
EC2_genes_unique <- pct_exp_EC2[which(pct_exp_EC2[,1] < 0.25),]

gene_list = c(rownames(EC1_genes_unique)[1:15], rownames(EC2_genes_unique)[1:15])
circle_heatmap_genelist_spec_cluster(dataset = young_noimmune, gene_list, n_genes = 15,
								clusters = c(2,5), rename_cluster = T, cluster_names = c("EC_1", "EC_2"),
								image_dpi = 600, image_width = 10, image_height = 5,chart_name = "",
								image_name = "EC1_vs_EC2 (unique).png")