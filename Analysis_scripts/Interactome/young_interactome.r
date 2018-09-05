# young_interactome.r
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/Interactome/young_interactome.r")


# This script imports and analyzes the young interactome.

# ---------------------------------------------------------------------------------

# First we load the data
load("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/interactome_v3/_inx.RData")

# First, we take the list of the node edges ("0~0", "0~1", ...)
age <- "young"
connection_list <- names(inx[[age]])

# We want to create 2 dataframes, one in which the direction of the interaction
# is important, i.e. "0 -> 1", is not the same as "1 -> 0", and the other on in which
# the directin is not important. At the same time, we want to save the results
# so we can analyze them later.

# in the case of the young sample, there are 9 cell clusters. So we make the 
# final dataframes
## for no directionality:
interactome_no_dir <- data.frame(matrix(0, nrow = 9, ncol = 9))
colnames(interactome_no_dir) <- c(0:8)
rownames(interactome_no_dir) <- c(0:8)

## directional interactome
interactome_dir <- interactome_no_dir

for (i in 1:length(connection_list)) {

	comp <- connection_list[i]
	
	# First, we find node names of genes detected in > 20% of cells in the cluster
	DRnodes <- inxNode[[age]][[comp]]$node[inxNode[[age]][[comp]]$detectRate > 0.2]
	# Ligand-receptor interactions where both genes were detected in > 20% of cells
	# in the cluster.
	interact_subset <- inx[[age]][[comp]][inx[[age]][[comp]]$nodeA %in% DRnodes & 
						inx[[age]][[comp]]$nodeB %in% DRnodes & 
						inx[[age]][[comp]]$direction %in% c("LtoR","RtoL"),]
	
	cluster_numbers <- strsplit(comp, "")
	nodeA <- as.integer(cluster_numbers[[1]][1])
	nodeB <- as.integer(cluster_numbers[[1]][3])
	
	# for the no direction interactome numbers:
	interactome_no_dir[nodeA + 1, nodeB + 1] <- nrow(interact_subset)
	interactome_no_dir[nodeB + 1, nodeA + 1] <- nrow(interact_subset)
	
	# for the directional interactome numbers:
	# for non-autocrine
	if (nodeA != nodeB) {
		interactome_dir[nodeA + 1, nodeB + 1] <- nrow(interact_subset[interact_subset$direction == "LtoR",])
		interactome_dir[nodeB + 1, nodeA + 1] <- nrow(interact_subset[interact_subset$direction == "RtoL",])
	}
	# for autocrine
	if (nodeA == nodeB) {
		interactome_dir[nodeA + 1, nodeB + 1] <- nrow(interact_subset)
	}
	
}

save(interactome_dir, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/interactome/young_interactome_dir.rdata")
save(interactome_no_dir, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/interactome/young_interactome_no_dir.rdata")