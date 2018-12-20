# MuSC_ligand_source_heatmap.R
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/Interactome/MuSC_ligand_source_heatmap.R")

# --------------------------------------------------------------------------------------------------------

# First we load the the functions we need.
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization_functions.r")

# --------------------------------------------------------------------------------------------------------

# Next we load the dataset
young <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_labeled.rds")

# Then we load the interactome data
load("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/interactome_v4/_inx.RData")

# Next we extract the ligands presented to the MuSCs from interactome
young_MuSC_interact <- cluster_interactome_rec_v2(inx, inxNode, cluster_id = "2", age = "young")

MuSC_heatmap <- ligand_source_heatmap(interac = young_MuSC_interact, labeled_dataset = young,
									cluster = "MuSCs")

#ds