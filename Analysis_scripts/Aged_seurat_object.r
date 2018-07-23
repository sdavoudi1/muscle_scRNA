# Aged_seurat_object.R
# Reads in and creates Seurat object from aged muscle scRNA

# This scripts takes the aged muscle scRNA data creates a Seurat object using the Seurat algorithm
# outlined in Butler et al, Nat Biotech, 2018.

# --------------------------------------------------------------------------------------------------------

# This script requires the Seurat library.
# library(Seurat)

# Next we load the functions we need
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.R")

# Read in the data from the aged expression dataset. We change the colnames since there are some cells that have
# the same cell.id in the aged and aged conditions. As a result, in the integrated analysis we change the 
# name of the columns here to avoid confusion later on. Additionally, the cell names will be the same in both
# the integrated analysis and here. "https://github.com/satijalab/seurat/issues/135"
aged.data <- Read10X(data.dir = "C:/Users/sadeg/Desktop/Muscle sc-RNA Seq/Aged/filtered_gene_bc_matrices/mm10/")
colnames(aged.data) = paste0("aged_", colnames(aged.data))

# Next we use the Seurat_from_10Xfile function to Read in the & expression datasets and create
# Seurat objects based on the filters, and normalize, scale and find the variable genes in each dataset.
aged <- Seurat_from_10Xfile(object = aged.data)
aged@meta.data$group <- "aged"