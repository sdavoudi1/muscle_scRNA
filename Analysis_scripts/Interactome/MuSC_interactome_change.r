# MuSC_interactome_change.R
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/Interactome/MuSC_interactome_change.R")

# This script analyzes the changes in the MuSC interactome.

# ---------------------------------------------------------------------------------------------------

if (!require("ks")) {install.packages("ks"); require(ks)}
library(ks)

if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
	library(Seurat)

# First we load the data
load("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/interactome_v3/_inx.RData")

# Then, we take name the connection we want in the case of the young and aged (i.e. for the MuSCs)
young_connection <- "2~2"
aged_connection <- "6~6"

# Next we take the list of the nodes that are present in 
DRnodes_young_MuSC <- inxNode[["young"]][[young_connection]]$node[inxNode[["young"]][[young_connection]]$detectRate > 0.25]
DRnodes_aged_MuSC <- inxNode[["aged"]][[aged_connection]]$node[inxNode[["aged"]][[aged_connection]]$detectRate > 0.25]

# Next we find the ligand receptor interactions where both genes were detected in > 25% of the cells
# cluster
young_MuSC_subset <- inx[["young"]][[young_connection]][inx[["young"]][[young_connection]]$nodeA %in% DRnodes_young_MuSC & 
													   inx[["young"]][[young_connection]]$nodeB %in% DRnodes_young_MuSC & 
													   inx[["young"]][[young_connection]]$direction %in% c("LtoR","RtoL"),]

aged_MuSC_subset <- inx[["aged"]][[aged_connection]][inx[["aged"]][[aged_connection]]$nodeA %in% DRnodes_aged_MuSC &
													inx[["aged"]][[aged_connection]]$nodeB %in% DRnodes_aged_MuSC &
													inx[["aged"]][[aged_connection]]$direction %in% c("LtoR", "RtoL"),]

# Next we need to compare the 2 lists to see which interactions are missing in each condition

# To do that, we create new dataframes with ligand -> receptor directionality.

# First we create a dataframe with the same number of rows as the young autocrine
# signaling dataset have. It will have 3 columns: ligand gene, receptor gene, 
# and the combined name. The data from young_MuSC_subset will be reorganized
# so they are all LtoR.
young_MuSC_autocrine <- data.frame(matrix(data = NA, nrow = nrow(young_MuSC_subset), ncol = 3))
colnames(young_MuSC_autocrine) <- c("Ligand","Receptor","LtoR")

for (i in 1:nrow(young_MuSC_subset)) {

	if (young_MuSC_subset$direction[i] == "LtoR") {
	
		young_MuSC_autocrine$Ligand[i] <- young_MuSC_subset$geneA[i]
		young_MuSC_autocrine$Receptor[i] <- young_MuSC_subset$geneB[i]
		young_MuSC_autocrine$LtoR[i] <- paste(young_MuSC_subset$geneA[i], young_MuSC_subset$geneB[i], sep = "_")
	
	}
	else {
	
		young_MuSC_autocrine$Ligand[i] <- young_MuSC_subset$geneB[i]
		young_MuSC_autocrine$Receptor[i] <- young_MuSC_subset$geneA[i]
		young_MuSC_autocrine$LtoR[i] <- paste(young_MuSC_subset$geneB[i], young_MuSC_subset$geneA[i], sep = "_")
	
	}

}

young_MuSC_autocrine <- unique(young_MuSC_autocrine)
rownames(young_MuSC_autocrine) <- young_MuSC_autocrine$LtoR

# We do the same for the aged MuSC autocrine signaling dataset as well.

aged_MuSC_autocrine <- data.frame(matrix(data = NA, nrow = nrow(aged_MuSC_subset), ncol = 3))
colnames(aged_MuSC_autocrine) <- c("Ligand", "Receptor", "LtoR")

for (i in 1:nrow(aged_MuSC_subset)) {

	if (aged_MuSC_subset$direction[i] == "LtoR") {
	
		aged_MuSC_autocrine$Ligand[i] <- aged_MuSC_subset$geneA[i]
		aged_MuSC_autocrine$Receptor[i] <- aged_MuSC_subset$geneB[i]
		aged_MuSC_autocrine$LtoR[i] <- paste(aged_MuSC_subset$geneA[i], aged_MuSC_subset$geneB[i], sep = "_")
	
	}
	else {
	
		aged_MuSC_autocrine$Ligand[i] <- aged_MuSC_subset$geneB[i]
		aged_MuSC_autocrine$Receptor[i] <- aged_MuSC_subset$geneA[i]
		aged_MuSC_autocrine$LtoR[i] <- paste(aged_MuSC_subset$geneB[i], aged_MuSC_subset$geneA[i], sep = "_")
	
	}

}

aged_MuSC_autocrine <- unique(aged_MuSC_autocrine)
rownames(aged_MuSC_autocrine) <- aged_MuSC_autocrine$LtoR

# -------------------------------------------------------------------------------

# Now, we go ahead and find the unique elements in each of the young and aged
# autocrine interactions.

unique_young_MuSC_autocrine <- young_MuSC_autocrine[!(young_MuSC_autocrine$LtoR %in% aged_MuSC_autocrine$LtoR),]
unique_aged_MuSC_autocrine <- aged_MuSC_autocrine[!(aged_MuSC_autocrine$LtoR %in% young_MuSC_autocrine$LtoR),]

# Next, we determine the differentially expressed genes in the MuSCs of aged and young mice.
# We load in the labeled muscle_combined_noimmune data.
muscle.combined <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/muscle_noimmune_combined_labeled.rds")
muscle.combined2 <- muscle.combined
# ----------------------------------------------------------------------------------

# To analyze 
muscle.combined2@meta.data$celltype.age <- paste0(muscle.combined2@ident, "_", muscle.combined2@meta.data$group)
muscle.combined2 <- StashIdent(muscle.combined2, save.name = "celltype")
muscle.combined2 <- SetAllIdent(muscle.combined2, id = "celltype.age")
aged.increase <- FindMarkers(muscle.combined2, ident.1 = "MuSC_aged", ident.2 = "MuSC_young", only.pos = T)
young.increase <- FindMarkers(muscle.combined2, ident.1 = "MuSC_young", ident.2 = "MuSC_aged", only.pos = T)
aged.increase$gene <- rownames(aged.increase)
young.increase$gene <- rownames(young.increase)

# ----------------------------------------------------------------------------------

# Now, for each of the unique_aged_MuSC_autocrine and unique_young_MuSC_autocrine dataframes,
# we want to check to see whether any of the genes are actually differntially expressed.

unique_DE_young_MuSC_autocrine <- unique_young_MuSC_autocrine[unique_young_MuSC_autocrine$Ligand %in% young.increase$gene & 
																unique_young_MuSC_autocrine$Receptor %in% young.increase$gene,]

# unique_DE_aged_MuSC_autocrine <- unique_aged_MuSC_autocrine[unique_aged_MuSC_autocrine$Ligand %in% aged.increase$gene & 
																# unique_aged_MuSC_autocrine$Receptor %in% aged.increase$gene,]
unique_RDE_aged_MuSC_autocrine <- unique_aged_MuSC_autocrine[unique_aged_MuSC_autocrine$Receptor %in% aged.increase$gene,]