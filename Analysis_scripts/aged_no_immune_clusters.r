# aged_no_immune_clusters.r
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/aged_no_immune_clusters.r")

# This script takes the clustered cells prepared using the Aged_clusters.r script
# removes the immune cells, and prepares new Seurat objects to be clustered. Then analyzes the data to
# create unbiased clusters.

# ------------------------------------------------------------------------------------------------------

if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
library(Seurat)
if(!require(plyr)) {install.packages("plyr"); require(plyr)}
library(plyr)

# Next we source the functions we need.
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")

# ------------------------------------------------------------------------------------------------------

# Begin with aged data. Need to remove cluster #6 (B_cells).

# First we load the data.
cat("Loading in aged Seurat object", "\n")
aged <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged.rds")
cat("----", "\n")

# Next we subset all clusters but the immune cells (cluster 6)
cat("Removing immune cell clusters", "\n")
aged_noimmune_wRBC <- SubsetData(object = aged, subset.name = "res.0.6", accept.value = c(0,1,2,3,4,5,7,8,9), do.clean = T)
cat("----", "\n")

# Now all we need to do, is the normalization and the Seurat_from_10Xfile.R script from the normalization step onward.
# The normalize_scale_fvg function will do that.
cat("Rescaling, normalizing the data, and finding variable genes", "\n")
aged_noimmune_wRBC <- normalize_scale_fvg(aged_noimmune_wRBC)
cat("----", "\n")

# To save the results preanalysis, we use the following code:
saveRDS(aged_noimmune_wRBC, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune_wRBC_preanalysis.rds")

# ------------------------------------------------------------------------------------------------------

# Perform linear dimensional reduction
cat("Performing linear dimensional reduction", "\n")
aged_noimmune_wRBC <- RunPCA(object = aged_noimmune_wRBC, pc.genes = aged_noimmune_wRBC@var.genes, do.print = FALSE)
aged_noimmune_wRBC <- ProjectPCA(object = aged_noimmune_wRBC, do.print = FALSE)
cat("----", "\n")

# Determine statistically significant principal components
cat("Determining statistically significant PCs (Jackstraw)", "\n")
aged_noimmune_wRBC <- JackStraw(object = aged_noimmune_wRBC, num.replicate = 100, display.progress = TRUE)
cat("----", "\n")

# # Plot the JackStraw plots to visualize the significant plots
# JackStrawPlot(object = aged_noimmune_wRBC, PCs = 1:20)

# Cluster cells
cat("Clustering cells", "\n")
aged_noimmune_wRBC <- FindClusters(object = aged_noimmune_wRBC, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
cat("----", "\n")

# # Clustering parameters
# PrintFindClustersParams(object = aged_noimmune_wRBC)

# Non-linear dimensional reduction (tSNE)
cat("Performing non-linear dimensional reduction", "\n")
aged_noimmune_wRBC <- RunTSNE(object = aged_noimmune_wRBC, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = aged_noimmune_wRBC)

# To save the results, we use the following code:
saveRDS(aged_noimmune_wRBC, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune_wRBC.rds")

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Upon completing the analysis, a new cluster emerged (cluster 8) which had high levels of expression of 
# Red blood cell genes (e.g. Gypa, Hbb-bt, Hba-a2, ...). Therefore, we go through the entire process again,
# this time excluding the RBCs.

# We subset all clusters but the immune cells (cluster 8)
cat("Removing RBC clusters", "\n")
aged_noimmune <- SubsetData(object = aged_noimmune_wRBC, subset.name = "res.0.6", accept.value = c(0,1,2,3,4,5,6,7), do.clean = T)
cat("----", "\n")

# Now all we need to do, is the normalization and the Seurat_from_10Xfile.R script from the normalization step onward.
# The normalize_scale_fvg function will do that.
cat("Rescaling, normalizing the data, and finding variable genes", "\n")
aged_noimmune <- normalize_scale_fvg(aged_noimmune)
cat("----", "\n")

saveRDS(aged_noimmune, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune_preanalysis.rds")

# ------------------------------------------------------------------------------------------------------

# Perform linear dimensional reduction
cat("Performing linear dimensional reduction", "\n")
aged_noimmune <- RunPCA(object = aged_noimmune, pc.genes = aged_noimmune@var.genes, do.print = FALSE)
aged_noimmune <- ProjectPCA(object = aged_noimmune, do.print = FALSE)
cat("----", "\n")

# Determine statistically significant principal components
cat("Determining statistically significant PCs (Jackstraw)", "\n")
aged_noimmune <- JackStraw(object = aged_noimmune, num.replicate = 100, display.progress = TRUE)
cat("----", "\n")

# # Plot the JackStraw plots to visualize the significant plots
# JackStrawPlot(object = aged_noimmune, PCs = 1:20)

# Cluster cells
cat("Clustering cells", "\n")
aged_noimmune <- FindClusters(object = aged_noimmune, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
cat("----", "\n")

# # Clustering parameters
# PrintFindClustersParams(object = aged_noimmune)

# Non-linear dimensional reduction (tSNE)
cat("Performing non-linear dimensional reduction", "\n")
aged_noimmune <- RunTSNE(object = aged_noimmune, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = aged_noimmune)

# To save the results, we use the following code:
saveRDS(aged_noimmune, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune.rds")

# # To rename the clusters to what we want:
# current.cluster.ids <- c(0,1,2,3,4,5,6)
# new.cluster.ids <- c("Myo_1", "FAP_1", "Myo_2", "EC", "FAP_2", "Tenocyte", "MuSC")
# aged_noimmune@ident <- plyr::mapvalues(x = aged_noimmune@ident, from = current.cluster.ids, to = new.cluster.ids)

# # To save the results, we use the following code:
# saveRDS(aged_noimmune, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune_labeled.rds")