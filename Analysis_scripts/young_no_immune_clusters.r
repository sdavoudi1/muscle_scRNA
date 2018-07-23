# young_no_immune_clusters.r
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/young_no_immune_clusters.r")

# This script takes the clustered cells prepared using the Young_cluster.r script,
# removes the immune cells, and prepares new Seurat object to be clustered. Then analyzes the data to
# create unbiased clusters.

# ------------------------------------------------------------------------------------------------------

if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
library(Seurat)
if(!require(plyr)) {install.packages("plyr"); require(plyr)}
library(plyr)

# Next we source the functions we need.
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")

# ------------------------------------------------------------------------------------------------------

# Begin with young data. Need to remove cluster #6 (B_cells).

# First we load the data.
cat("Loading in young Seurat object", "\n")
young <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young.rds")
cat("----", "\n")

# Next we subset all clusters but the immune cells (cluster 6)
cat("Removing immune cell clusters", "\n")
young_noimmune <- SubsetData(object = young, subset.name = "res.0.6", accept.value = c(0,1,2,3,4,5,7,8), do.clean = T)
cat("----", "\n")

# Now all we need to do, is the normalization and the Seurat_from_10Xfile.R script from the normalization step onward.
# The normalize_scale_fvg function will do that.
cat("Rescaling, normalizing the data, and finding variable genes", "\n")
young_noimmune <- normalize_scale_fvg(young_noimmune)
cat("----", "\n")

# We save the raw un-analyzed Seurat object
saveRDS(young_noimmune, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_preanalysis.rds")

# ------------------------------------------------------------------------------------------------------

# Perform linear dimensional reduction
cat("Performing linear dimensional reduction", "\n")
young_noimmune <- RunPCA(object = young_noimmune, pc.genes = young_noimmune@var.genes, do.print = FALSE)
young_noimmune <- ProjectPCA(object = young_noimmune, do.print = FALSE)
cat("----", "\n")

# Determine statistically significant principal components
cat("Determining statistically significant PCs (Jackstraw)", "\n")
young_noimmune <- JackStraw(object = young_noimmune, num.replicate = 100, display.progress = TRUE)
cat("----", "\n")

# # Plot the JackStraw plots to visualize the significant plots
# JackStrawPlot(object = young_noimmune, PCs = 1:20)

# Cluster cells
cat("Clustering cells", "\n")
young_noimmune <- FindClusters(object = young_noimmune, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
cat("----", "\n")

# # Clustering parameters
# PrintFindClustersParams(object = young_noimmune)

# Non-linear dimensional reduction (tSNE)
cat("Performing non-linear dimensional reduction", "\n")
young_noimmune <- RunTSNE(object = young_noimmune, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = young_noimmune)

# To save the results, we use the following code:
saveRDS(young_noimmune, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune.rds")

# # To rename the clusters to what we want:
# current.cluster.ids <- c(0,1,2,3,4,5,6,7,8)
# new.cluster.ids <- c("FAP_1", "EC_1", "MuSC", "FAP_2", "EC_2", "FAP_3", "Tenocyte", "FAP_4", "Schwann")
# young_noimmune@ident <- plyr::mapvalues(x = young_noimmune@ident, from = current.cluster.ids, to = new.cluster.ids)

# # To save the results, we use the following code:
# saveRDS(young_noimmune, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_labeled.rds")