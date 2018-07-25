# This script is used to perform a portion of the enrichment analysis on the young dataset

# Young analysis

# ----------------------------------------------------------------------------------

#REQUIREMENT: Load the young Seurat object first.

# First we load the library and source code for some of the functions we're using.
if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
library(Seurat)

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.R")

# First, we use the FindAllMarkers from the Seurat package to find the differentially expressed,
# genes in each cluster.
aged_diff_markers <- FindAllMarkers(object = aged, only.pos = T, min.pct = 0.25, thresh.use = 0.25)

# Next, we determine the differentially expressed genes for each cluster using the
# following function.

diff_genes_per_cluster_txt(dataset_name = "aged", diff_genes = aged_diff_markers ,output_dir = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/Aged_Enrichment_Analysis/")

# ----------------------------------------------------------------------------------