# Aged_analysis.R
# Analysis of aged muscle scRNA

# This scripts takes the aged muscle scRNA data and combines them using the Seurat algorithm
# outlined in Butler et al, Nat Biotech, 2018.

# --------------------------------------------------------------------------------------------------------

# We load the required libraries.
library(Seurat)

# # Next we load the functions we need
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.R")

# Next we create the aged and young Seurat objects using the scripts below:
cat("Creating aged Seurat Object", "\n")
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Aged_seurat_object.r")
cat("----", "\n")

# Perform linear dimensional reduction
cat("Performing linear dimensional reduction", "\n")
aged <- RunPCA(object = aged, pc.genes = aged@var.genes, do.print = FALSE)
aged <- ProjectPCA(object = aged, do.print = FALSE)
cat("----", "\n")

# Determine statistically significant principal components
cat("Determining statistically significant PCs (Jackstraw)", "\n")
aged <- JackStraw(object = aged, num.replicate = 100, display.progress = TRUE)
cat("----", "\n")

# # Plot the JackStraw plots to visualize the correct plots
# JackStrawPlot(object = aged, PCs = 1:20)

# Cluster cells
cat("Clustering cells", "\n")
aged <- FindClusters(object = aged, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
cat("----", "\n")

# # Clustering parameters
# PrintFindClustersParams(object = aged)

# Non-linear dimensional reduction (tSNE)
cat("Performing non-linear dimensional reduction", "\n")
aged <- RunTSNE(object = aged, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = aged)

# To save the results, we use the following code:
saveRDS(aged, file = "C:/Users/sadeg/Google Drive/scRNA/data/Aged_analysis/aged.rds")