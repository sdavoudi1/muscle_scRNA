# This script analyzes the EC_1 and EC_2 subclusters for further comparison.
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/subclusters/EC_subclustering.r")

# ----------------------------------------------------------------------------------------------

# First we load the required libraries.
if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
library(Seurat)
if(!require(plyr)) {install.packages("plyr"); require(plyr)}
library(plyr)

# Next we source the functions we need.
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")

# ----------------------------------------------------------------------------------------------

# Next we load the dataset
young_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune.rds")

# ----------------------------------------------------------------------------------------------

# IF WE DON'T WANT SUBCLUSTERING, JUST USE THESE 2 LINES. THE REST OF THE SCRIPT RECLUSTERS THE WHOLE
# DATASET USING THE EC_1 AND EC_2 POPULATIONS.

# Next we determine the differentially expressed genes between EC_1 and EC_2
# EC1_diff <- FindMarkers(object = young_noimmune, ident.1 = 1, ident.2 = 4, only.pos = T, min.pct = 0.25)
# EC2_diff <- FindMarkers(object = young_noimmune, ident.1 = 4, ident.2 = 1, only.pos = T, min.pct = 0.25)

# ----------------------------------------------------------------------------------------------

# Next we subset the EC clusters (cluster 1,4)
cat("Subsetting Endothelial cells", "\n")
Endothelial_cells <- SubsetData(object = young_noimmune, subset.name = "res.0.6", accept.value = c(1,4), do.clean = T)
cat("----", "\n")

# Now all we need to do, is the normalization and the Seurat_from_10Xfile.R script from the normalization step onward.
# The normalize_scale_fvg function will do that.
cat("Rescaling, normalizing the data, and finding variable genes", "\n")
Endothelial_cells <- normalize_scale_fvg(Endothelial_cells)
cat("----", "\n")

# We save the raw un-analyzed Seurat object
# saveRDS(Endothelial_cells, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/Endothelial_cells_preanalysis.rds")

# ----------------------------------------------------------------------------------------------

# Perform linear dimensional reduction
cat("Performing linear dimensional reduction", "\n")
Endothelial_cells <- RunPCA(object = Endothelial_cells, pc.genes = Endothelial_cells@var.genes, do.print = FALSE)
Endothelial_cells <- ProjectPCA(object = Endothelial_cells, do.print = FALSE)
cat("----", "\n")

# Determine statistically significant principal components
cat("Determining statistically significant PCs (Jackstraw)", "\n")
Endothelial_cells <- JackStraw(object = Endothelial_cells, num.replicate = 100, display.progress = TRUE)
cat("----", "\n")

# # Plot the JackStraw plots to visualize the significant plots. The last significant PC is PC11.
# JackStrawPlot(object = Endothelial_cells, PCs = 1:20)

# Cluster cells
cat("Clustering cells", "\n")
Endothelial_cells <- FindClusters(object = Endothelial_cells, reduction.type = "pca", dims.use = 1:11, resolution = 0.6, print.output = 0, save.SNN = TRUE)
cat("----", "\n")

# ----------------------------------------------------------------------------------------------

# Non-linear dimensional reduction (tSNE)
cat("Performing non-linear dimensional reduction", "\n")
Endothelial_cells <- RunTSNE(object = Endothelial_cells, dims.use = 1:11, do.fast = TRUE)
TSNEPlot(object = Endothelial_cells)

# To save the results, we use the following code:
saveRDS(Endothelial_cells, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/subclusters/Endothelial_cells_reclustering.rds")