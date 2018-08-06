# circle_heatmap_5gene_noimmune.R

# -----------------------------------------------------------------------------------------------------

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization_functions.r")
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")


# -----------------------------------------------------------------------------------------------------

# Young_noimmune

# First we load the required datasets
young_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune.rds")
young_noimmune_diff_markers <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_diff_markers.rds")

# Then we assign cluster order and names
cluster_order <- c(2, 5, 3, 1, 4, 6, 8, 9, 7)
cluster_names <- c("EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "FAP_4", "Schwann", "Tenocyte")

# Then we use the circle_heatmap function to draw the graph.
young_noimmune_circle_heatmap_5genes <- circle_heatmap(dataset = young_noimmune, diff_markers = young_noimmune_diff_markers,
							image_name = "young_noimmune_circle_heatmap.png", image_dpi = 600, 
							image_width = 10, image_height = 5, n_genes = 5, cluster_names = cluster_names, 
							cluster_order = cluster_order)
								
# -----------------------------------------------------------------------------------------------------

# Aged_noimmune

# First we load the required datasets
aged_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune.rds")
aged_noimmune_diff_markers <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune_diff_markers.rds")

# Then we assign cluster order and names
cluster_order <- c(4, 7, 1, 3, 2, 5, 6)
cluster_names <- c("EC", "MuSC", "Myo_1", "Myo_2", "FAP_1", "FAP_2", "Tenocyte")

# Then we use the circle_heatmap function to draw the graph.
aged_noimmune_circle_heatmap_5genes <- circle_heatmap(dataset = aged_noimmune, diff_markers = aged_noimmune_diff_markers,
							image_name = "aged_noimmune_circle_heatmap.png", image_dpi = 600, 
							image_width = 10, image_height = 5, n_genes = 5, cluster_names = cluster_names, 
							cluster_order = cluster_order)