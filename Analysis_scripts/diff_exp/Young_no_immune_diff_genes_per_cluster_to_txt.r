# This script is used to perform a portion of the enrichment analysis on the young dataset
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/diff_exp/Young_no_immune_diff_genes_per_cluster_to_txt.r")

# Young_noimmune analysis

# ----------------------------------------------------------------------------------

# First we load the library and source code for some of the functions we're using.
if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
library(Seurat)

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.R")

# Next, we load the differentially expressed genes in each cluster.
young_noimmune_diff_markers <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_diff_markers.rds")

# Next, we determine the differentially expressed genes for each cluster using the
# following function.

diff_genes_per_cluster_txt(dataset_name = "young_noimmune", dataset = young_noimmune_diff_markers ,output_dir = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/young_noimmune_enrichment_analysis/")

# ----------------------------------------------------------------------------------