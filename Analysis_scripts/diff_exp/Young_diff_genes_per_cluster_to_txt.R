# This script is used to perform a portion of the enrichment analysis on the young dataset

# Young analysis

# ----------------------------------------------------------------------------------

#REQUIREMENT: Load the young Seurat object first.

# First we load the library and source code for some of the functions we're using.
if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
library(Seurat)

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.R")

# Next, we load the differentially expressed genes in each cluster.
young_diff_markers <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_diff_markers.rds"))

# Next, we determine the differentially expressed genes for each cluster using the
# following function.

diff_genes_per_cluster_txt(dataset_name = "young", diff_genes = young_diff_markers ,output_dir = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/Young_Enrichment_Analysis/_Enrichment_Analysis/")

# ----------------------------------------------------------------------------------