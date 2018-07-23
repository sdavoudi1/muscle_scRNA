# Young_analysis.R
# Analysis of young muscle scRNA

# This scripts takes the young muscle scRNA data and combines them using the Seurat algorithm
# outlined in Butler et al, Nat Biotech, 2018.

# --------------------------------------------------------------------------------------------------------

# We load the required libraries.
if(!require(Seurat)) {
install.packages("Seurat"); require(Seurat)}

library(Seurat)

# # Next we load the functions we need
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.R")

# Next we create the young Seurat objects using the scripts below:
cat("Creating young Seurat Object", "\n")
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/Young_seurat_object.r")
cat("----", "\n")

# Perform linear dimensional reduction
cat("Performing linear dimensional reduction", "\n")
young <- RunPCA(object = young, pc.genes = young@var.genes, do.print = FALSE)
young <- ProjectPCA(object = young, do.print = FALSE)
cat("----", "\n")

# Determine statistically significant principal components
cat("Determining statistically significant PCs (Jackstraw)", "\n")
young <- JackStraw(object = young, num.replicate = 100, display.progress = TRUE)
cat("----", "\n")

# # Plot the JackStraw plots to visualize the correct plots
# JackStrawPlot(object = young, PCs = 1:20)

# Cluster cells
cat("Clustering cells", "\n")
young <- FindClusters(object = young, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
cat("----", "\n")

# # Clustering parameters
# PrintFindClustersParams(object = young)

# Non-linear dimensional reduction (tSNE)
cat("Performing non-linear dimensional reduction", "\n")
young <- RunTSNE(object = young, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = young)

# To save the results, we use the following code:
saveRDS(young, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young.rds")