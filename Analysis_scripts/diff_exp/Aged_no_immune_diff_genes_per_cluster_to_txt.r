# This script is used to perform a portion of the enrichment analysis on the aged_noimmune dataset
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/diff_exp/Aged_no_immune_diff_genes_per_cluster_to_txt.r")

# Aged_noimmune analysis

# ----------------------------------------------------------------------------------

# First we load the library and source code for some of the functions we're using.
if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
library(Seurat)

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.R")

# Next, we load the differentially expressed genes in each cluster.
aged_noimmune_diff_markers <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune_diff_markers.rds")


# Next, we determine the differentially expressed genes for each cluster using the
# following function.

diff_genes_per_cluster_txt(dataset_name = "aged_noimmune", dataset = aged_noimmune_diff_markers ,output_dir = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/aged_noimmune_enrichment_analysis/")

# ----------------------------------------------------------------------------------