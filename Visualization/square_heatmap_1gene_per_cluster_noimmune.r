# heatmap_1gene_per_cluster_noimmune


# This script creates the square heatmaps with 1 gene for each cluster, normalized between 0 and 1.


# -----------------------------------------------------------------------------------------------------

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization_functions.r")

# -----------------------------------------------------------------------------------------------------

# YOUNG_NOIMMUNE

young_noimmune_labeled <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_labeled.rds")
genes <- c("Cdh5", "Pecam1", "Pax7", "Cd34", "Col1a1", "Ly6a", "Pdgfra", "Plp1", "Tnmd")
cluster_order <- c("EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "FAP_4", "Schwann", "Tenocyte")
square_hmap_young_noimmune <- square_heatmap(young_noimmune_labeled, genes = genes, cluster_order = cluster_order)
saveRDS(square_hmap_young_noimmune, "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/graphs/young_noimmune_1gene_heatmap.RDS")

# -----------------------------------------------------------------------------------------------------

# AGED_NOIMMUNE

aged_noimmune_labeled <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune_labeled.rds")
genes <- c("Pecam1", "Pax7", "Myod1", "Des", "Col1a1", "Pdgfra", "Tnmd")
cluster_order <- c("EC", "MuSC", "Myo_1", "Myo_2", "FAP_1", "FAP_2", "Tenocyte")
square_hmap_aged_noimmune <- square_heatmap(aged_noimmune_labeled, genes = genes, cluster_order = cluster_order)
saveRDS(square_hmap_aged_noimmune, "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/graphs/aged_noimmune_1gene_heatmap.RDS")

# -----------------------------------------------------------------------------------------------------