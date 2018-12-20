# young_interactome_v2.r
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/Interactome/young_interactome_v2.r")


# This script imports and analyzes the young interactome after merging FAP_1 and FAP_3

# ---------------------------------------------------------------------------------

# First we load the data
load("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/interactome_v4/_inx.RData")

# First, we take the list of the node edges ("0~0", "0~1", ...)
age <- "young"
connection_list <- names(inx[[age]])

# ----------------------------------------------------------------------------------

# We want to create 2 dataframes, one in which the direction of the interaction
# is important, i.e. "0 -> 1", is not the same as "1 -> 0", and the other on in which
# the direction is not important. At the same time, we want to save the results
# so we can analyze them later.

# in the case of the young sample, there are 8 cell clusters. So we make the 
# final dataframes
## for no directionality:
interactome_no_dir_v2 <- data.frame(matrix(0, nrow = 8, ncol = 8))
colnames(interactome_no_dir_v2) <- c(0:7)
rownames(interactome_no_dir_v2) <- c(0:7)

## directional interactome
interactome_dir_v2 <- interactome_no_dir_v2

for (i in 1:length(connection_list)) {

	comp <- connection_list[i]
	
	# First, we find node names of genes detected in > 20% of cells in the cluster
	DRnodes <- inxNode[[age]][[comp]]$node[inxNode[[age]][[comp]]$DR > 0.2]
	# Ligand-receptor interactions where both genes were detected in > 20% of cells
	# in the cluster.
	interact_subset <- inx[[age]][[comp]][inx[[age]][[comp]]$nodeA %in% DRnodes & 
						inx[[age]][[comp]]$nodeB %in% DRnodes & 
						inx[[age]][[comp]]$direction %in% c("LtoR","RtoL"),]
	
	cluster_numbers <- strsplit(comp, "")
	nodeA <- as.integer(cluster_numbers[[1]][1])
	nodeB <- as.integer(cluster_numbers[[1]][3])
	
	# for the no direction interactome numbers:
	# for non-autocrine
	if (nodeA != nodeB) {
		interactome_no_dir_v2[nodeA + 1, nodeB + 1] <- nrow(interact_subset)
		interactome_no_dir_v2[nodeB + 1, nodeA + 1] <- nrow(interact_subset)
	}
	# for autocrine (the reason it is divided by 2 is that it will count the 
	# interactions twice, once as LtoR and once as RtoL.
	if (nodeA == nodeB) {
		interactome_no_dir_v2[nodeA + 1, nodeB + 1] <- nrow(interact_subset)/2
	}
	
	# for the directional interactome numbers:
	# for non-autocrine
	if (nodeA != nodeB) {
		interactome_dir_v2[nodeA + 1, nodeB + 1] <- nrow(interact_subset[interact_subset$direction == "LtoR",])
		interactome_dir_v2[nodeB + 1, nodeA + 1] <- nrow(interact_subset[interact_subset$direction == "RtoL",])
	}
	# for autocrine (the reason it is divided by 2 is that it will count the 
	# interactions twice, once as LtoR and once as RtoL.
	if (nodeA == nodeB) {
		interactome_dir_v2[nodeA + 1, nodeB + 1] <- nrow(interact_subset)/2
	}
	
}

save(interactome_dir_v2, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/interactome/young_interactome_dir_v2.rdata")
save(interactome_no_dir_v2, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/interactome/young_interactome_no_dir_v2.rdata")

# ----------------------------------------------------------------------------------