# young_noimmune_interactome.R
# source("C:/users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization/young_noimmune_interactome.R")

# This script creates the young_noimmune interactome network.

# ----------------------------------------------------------------------------------------------------

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization_functions.r")


# First we load the directed network matrix.
load("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/interactome/25% expression/young_interactome_dir.rdata")

# Next we rearrange the interactome network in the order that we want.
order <- c(1,4,2,0,3,5,7,8,6) + 1
interactome_dir2 <- interactome_dir[order,]
interactome_dir2 <- interactome_dir2[, order]

# Then we draw the figures.

# For the full interactome we use the following
igraph_circle_network_full(interactome_dir2, label_v = F,
							save_pdf = T, plot_name = "C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/young_noimmune_full_interactome.pdf")
							
# For the single cell interactome we have 2 versions, the signals received by each 
# cluster (columns) vs the signals given by each cluster (rows).

# Signals received
for (i in 1:9) {
	interactome_dir_received <- data.frame(matrix(0, nrow = 9, ncol = 9))
	interactome_dir_received[,i] <- interactome_dir2[,i]
	rownames(interactome_dir_received) <- rownames(interactome_dir2)
	colnames(interactome_dir_received) <- colnames(interactome_dir2)
	plotname <- paste("C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/young_recipient_interactome_", as.character(i), ".pdf", sep = "")
	igraph_circle_network_full(interactome_dir_received, label_v = T, 
							   start_curve = 0, save_pdf = T,
							   plot_name = plotname)
}