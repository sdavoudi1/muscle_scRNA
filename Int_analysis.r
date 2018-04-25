# Int_analysis.R
# Integrated Analysis of young and aged muscle scRNA

# This scripts takes the young and aged muscle scRNA data and combines them using the Seurat algorithm
# outlined in Butler et al, Nat Biotech, 2018.

# --------------------------------------------------------------------------------------------------------

# We load the required libraries.
library(Seurat)

# Next we load the functions we need
source("Int_analysis_functions.R")

# Use the Seurat_from_10Xfile function to Read in the young (healthy) and aged expression datasets and create
# Seurat objects based on the filters
young <- Seurat_from_10Xfile(directory = "C:/Users/sadeg/Desktop/Muscle sc-RNA Seq/Healthy/Gilbert_Sadegh__Healthy/filtered_gene_bc_matrices/mm10/")
aged <- Seurat_from_10Xfile(directory = "C:/Users/sadeg/Desktop/Muscle sc-RNA Seq/Aged/filtered_gene_bc_matrices/mm10/")
