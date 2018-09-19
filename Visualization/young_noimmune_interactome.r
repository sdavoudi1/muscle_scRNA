# young_noimmune_interactome.R
# source("C:/users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization/young_noimmune_interactome.R

# This script creates the young_noimmune interactome network.

# ----------------------------------------------------------------------------------------------------

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization_functions.r")


# First we load the directed network matrix.
load("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/interactome/25% expression/young_interactome_dir.rdata")

# Then we draw the figures.
igraph_circle_network_full(interactome_dir,save_pdf = T, plot_name = "C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/young_noimmune_full_interactome.pdf")