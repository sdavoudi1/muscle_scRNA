# circle_heatmap_young_noimmune_Schwann_top10genes.R
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization/circle_heatmap_young_noimmune_Schwann_top10genes.r")

# -----------------------------------------------------------------------------------------------------

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization_functions.r")
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")


# -----------------------------------------------------------------------------------------------------

# Young_noimmune

# First we load the required datasets
young_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune.rds")
young_noimmune_diff_markers <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_diff_markers.rds")

# Then we assign cluster order and names
cluster_order <- c(9, 2, 5, 3, 1, 4, 6, 8, 7)
cluster_names <- c("Schwann", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "FAP_4", "Tenocytes")

# Then we determine the list of genes we want to draw. We use the top 10 differentially expressed genes
# in the Schwann cell cluster.
Schwann_noimmune_diff_genes <- young_noimmune_diff_markers$gene[which(young_noimmune_diff_markers$cluster == 8)]
gene_list <- Schwann_noimmune_diff_genes[1:10]

# Then we use the circle_heatmap function to draw the graph.
circle_heatmap_young_noimmune_Schwann_top10 <- circle_heatmap_genelist(dataset = young_noimmune, gene_list = gene_list,
											image_name = "image.png", image_dpi = 600, image_width = 5,
											image_height = 5, reorder_cluster = T, cluster_order = cluster_order, 
											rename_cluster = T, cluster_names = cluster_names)