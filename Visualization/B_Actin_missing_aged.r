# B_Actin_missing.r
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization/B_Actin_missing_aged.r")

# This script creates the dot plot for the missing B_Actin in the myogenic cells
# in the aged sample.

# ---------------------------------------------------------------------------------------------------

aged_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune.rds")

aged_b_actin <- FeaturePlot(object = aged_noimmune, features.plot = c("Actb"), cols.use = c("blue", "yellow"), 
							reduction.use = "tsne")
							
young_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune.rds")

young_b_actin <- FeaturePlot(object = young_noimmune, features.plot = c("Actb"), cols.use = c("blue", "yellow"), 
							reduction.use = "tsne")