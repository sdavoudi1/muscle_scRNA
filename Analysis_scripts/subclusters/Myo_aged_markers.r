# Myo_aged_markers.R

# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/subclusters/Myo_aged_markers.r")

# ---------------------------------------------------------------------------------------

# We first load the aged sample.
aged <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune.rds")

# ---------------------------------------------------------------------------------------

# For the initial analysis, we will combine the 2 myo populations.

current.cluster.ids <- c(0,1,2,3,4,5,6)
new.cluster.ids <- c(0,1,0,3,4,5,6)
aged@ident <- plyr::mapvalues(x = aged@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = aged, do.label = T, pt.size = 0.5)

# Next we find the differentially expressed genes in the Myo population
Myo.markers <- FindMarkers(object = aged, ident.1 = 0, min.pct = 0.25)