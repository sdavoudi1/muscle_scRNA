# square_heatmap_growth_factors_young_noimmune.r
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization/square_heatmap_growth_factors_young_noimmune.r")


# This script creates the square heatmaps with list of genes for each cluster, normalized between 0 and 1.


# -----------------------------------------------------------------------------------------------------

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization_functions.r")

# -----------------------------------------------------------------------------------------------------

# YOUNG_NOIMMUNE

young_noimmune_labeled <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_labeled.rds")
genes <- c("Fst", "Igf1", "Fgf2", "Ngf", "Angpt1", "Vegfa", "Thbs1")
cluster_order <- c("EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "FAP_4", "Schwann", "Tenocyte")
square_hmap_growth_factors_young_noimmune <- square_heatmap(young_noimmune_labeled, genes = genes, cluster_order = cluster_order)
saveRDS(square_hmap_growth_factors_young_noimmune, "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/graphs/young_noimmune_growth_factors_heatmap.RDS")

# -----------------------------------------------------------------------------------------------------