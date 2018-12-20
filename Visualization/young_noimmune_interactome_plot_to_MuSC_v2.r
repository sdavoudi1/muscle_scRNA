# young_noimmune_interactome_plot_v2.R
# source("C:/users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization/young_noimmune_interactome_plot_v2.R")

# This script creates the young_noimmune interactome network.

# ----------------------------------------------------------------------------------------------------

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization_functions.r")


# First we load the directed network matrix.
load("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/interactome/20% expression/young_interactome_dir_v2.rdata")

# Next we rearrange the interactome network in the order that we want.
order <- c(1,4,2,0,3,6,5,7) + 1
interactome_dir2 <- interactome_dir_v2[order,]
interactome_dir2 <- interactome_dir2[, order]

# Then we draw the figures.

# For the interactome of the MuSCs we use the following
igraph_circle_net_lig_to_cluster(interactome_dir_v2, order = order, cluster = 2,
									label_v = T, save_pdf = F)
									
igraph_circle_net_lig_to_cluster(interactome_dir_v2, order = order, cluster = 2, label_v = F,
									save_pdf = T, plot_name = "C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/young_noimmune_lig_to_MuSC.pdf")

# ----------