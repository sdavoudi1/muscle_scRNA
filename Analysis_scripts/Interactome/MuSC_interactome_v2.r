# MuSC_interactome_v2.R
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/Interactome/MuSC_interactome_v2.R")

# -----------------------------------------------------------------------------------------------------

# First we load the functions we need
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")
	
# Then we load the interactome data
load("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/interactome_v4/_inx.RData")

# then we use the cluster_interactome_rec functions to determine all the signals the MuSCs are receiving
young_MuSC_interact <- cluster_interactome_rec_v2(inx, inxNode, cluster_id = "2", age = "young")
#aged_MuSC_interact <- cluster_interactome_rec(inx, inxNode, cluster_id = "6", age = "aged")

# To save these results, we replace the cluster number with their ids
# young
young_MuSC_interact2 <- young_MuSC_interact
if (nrow(young_MuSC_interact2) > 0) {
	young_MuSC_interact2[young_MuSC_interact2$cellTypeB == 2, "cellTypeB"] <- "MuSC"

	young_MuSC_interact2[young_MuSC_interact2$cellTypeA == 0, "cellTypeA"] <- "FAP_1"
	young_MuSC_interact2[young_MuSC_interact2$cellTypeA == 1, "cellTypeA"] <- "EC_1"
	young_MuSC_interact2[young_MuSC_interact2$cellTypeA == 2, "cellTypeA"] <- "MuSC"
	young_MuSC_interact2[young_MuSC_interact2$cellTypeA == 3, "cellTypeA"] <- "FAP_2"
	young_MuSC_interact2[young_MuSC_interact2$cellTypeA == 4, "cellTypeA"] <- "EC_2"
	young_MuSC_interact2[young_MuSC_interact2$cellTypeA == 5, "cellTypeA"] <- "Tenocyte"
	young_MuSC_interact2[young_MuSC_interact2$cellTypeA == 6, "cellTypeA"] <- "FAP_3"
	young_MuSC_interact2[young_MuSC_interact2$cellTypeA == 7, "cellTypeA"] <- "Schwann"
}
write.csv(young_MuSC_interact2, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/young_MuSC_interact_v2.csv")
rm(young_MuSC_interact2)

# aged
# aged_MuSC_interact2 <- aged_MuSC_interact
# if (nrow(aged_MuSC_interact2 > 0)) {
	# aged_MuSC_interact2[aged_MuSC_interact2$cellTypeB == 6, "cellTypeB"] <- "MuSC"

	# aged_MuSC_interact2[aged_MuSC_interact2$cellTypeA == 0, "cellTypeA"] <- "Myo_1"
	# aged_MuSC_interact2[aged_MuSC_interact2$cellTypeA == 1, "cellTypeA"] <- "FAP_1"
	# aged_MuSC_interact2[aged_MuSC_interact2$cellTypeA == 2, "cellTypeA"] <- "Myo_2"
	# aged_MuSC_interact2[aged_MuSC_interact2$cellTypeA == 3, "cellTypeA"] <- "EC"
	# aged_MuSC_interact2[aged_MuSC_interact2$cellTypeA == 4, "cellTypeA"] <- "FAP_2"
	# aged_MuSC_interact2[aged_MuSC_interact2$cellTypeA == 5, "cellTypeA"] <- "Tenocyte"
	# aged_MuSC_interact2[aged_MuSC_interact2$cellTypeA == 6, "cellTypeA"] <- "MuSC"
# }
# write.csv(aged_MuSC_interact2, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/aged_MuSC_interact.csv")
# rm(aged_MuSC_interact2)

# ------------------------------------------------------------------------------------------------------

# Next, we find the unique interactions in each of the young and aged datasets
# unq_young_MuSC_interact <- young_MuSC_interact[!(young_MuSC_interact$interaction %in% aged_MuSC_interact$interaction),]
# unq_young_MuSC_interact2 <- unq_young_MuSC_interact
# unq_aged_MuSC_interact <- aged_MuSC_interact[!(aged_MuSC_interact$interaction %in% young_MuSC_interact$interaction),]
# unq_aged_MuSC_interact2 <- unq_aged_MuSC_interact


# # To save these results, we replace the cluster number with their ids
# # young
# if (nrow(unq_young_MuSC_interact2 > 0)) {
	# unq_young_MuSC_interact2[unq_young_MuSC_interact2$cellTypeB == 2, "cellTypeB"] <- "MuSC"

	# unq_young_MuSC_interact2[unq_young_MuSC_interact2$cellTypeA == 0, "cellTypeA"] <- "FAP_1"
	# unq_young_MuSC_interact2[unq_young_MuSC_interact2$cellTypeA == 1, "cellTypeA"] <- "EC_1"
	# unq_young_MuSC_interact2[unq_young_MuSC_interact2$cellTypeA == 2, "cellTypeA"] <- "MuSC"
	# unq_young_MuSC_interact2[unq_young_MuSC_interact2$cellTypeA == 3, "cellTypeA"] <- "FAP_2"
	# unq_young_MuSC_interact2[unq_young_MuSC_interact2$cellTypeA == 4, "cellTypeA"] <- "EC_2"
	# unq_young_MuSC_interact2[unq_young_MuSC_interact2$cellTypeA == 5, "cellTypeA"] <- "FAP_3"
	# unq_young_MuSC_interact2[unq_young_MuSC_interact2$cellTypeA == 6, "cellTypeA"] <- "Tenocyte"
	# unq_young_MuSC_interact2[unq_young_MuSC_interact2$cellTypeA == 7, "cellTypeA"] <- "FAP_4"
	# unq_young_MuSC_interact2[unq_young_MuSC_interact2$cellTypeA == 8, "cellTypeA"] <- "Schwann"
# }
# write.csv(unq_young_MuSC_interact2, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/unique_young_MuSC_interact.csv")
# rm(unq_young_MuSC_interact2)

# # aged
# if (nrow(unq_aged_MuSC_interact2 > 0)) {
	# unq_aged_MuSC_interact2[unq_aged_MuSC_interact2$cellTypeB == 6, "cellTypeB"] <- "MuSC"

	# unq_aged_MuSC_interact2[unq_aged_MuSC_interact2$cellTypeA == 0, "cellTypeA"] <- "Myo_1"
	# unq_aged_MuSC_interact2[unq_aged_MuSC_interact2$cellTypeA == 1, "cellTypeA"] <- "FAP_1"
	# unq_aged_MuSC_interact2[unq_aged_MuSC_interact2$cellTypeA == 2, "cellTypeA"] <- "Myo_2"
	# unq_aged_MuSC_interact2[unq_aged_MuSC_interact2$cellTypeA == 3, "cellTypeA"] <- "EC"
	# unq_aged_MuSC_interact2[unq_aged_MuSC_interact2$cellTypeA == 4, "cellTypeA"] <- "FAP_2"
	# unq_aged_MuSC_interact2[unq_aged_MuSC_interact2$cellTypeA == 5, "cellTypeA"] <- "Tenocyte"
	# unq_aged_MuSC_interact2[unq_aged_MuSC_interact2$cellTypeA == 6, "cellTypeA"] <- "MuSC"
# }
# write.csv(unq_aged_MuSC_interact2, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/unique_aged_MuSC_interact.csv")
# rm(unq_aged_MuSC_interact2)

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------

# Next, we want to look at the interactions in which the receptor in the MuSCs is DE in the aged or young.

# We determine the differentially expressed genes in the MuSCs of aged and young mice.
# We load in the labeled muscle_combined_noimmune data.
# muscle.combined <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/muscle_noimmune_combined_labeled.rds")
# muscle.combined2 <- muscle.combined

# # To analyze 
# muscle.combined2@meta.data$celltype.age <- paste0(muscle.combined2@ident, "_", muscle.combined2@meta.data$group)
# muscle.combined2 <- StashIdent(muscle.combined2, save.name = "celltype")
# muscle.combined2 <- SetAllIdent(muscle.combined2, id = "celltype.age")
# aged.increase <- FindMarkers(muscle.combined2, ident.1 = "MuSC_aged", ident.2 = "MuSC_young", only.pos = T)
# young.increase <- FindMarkers(muscle.combined2, ident.1 = "MuSC_young", ident.2 = "MuSC_aged", only.pos = T)
# aged.increase$gene <- rownames(aged.increase)
# young.increase$gene <- rownames(young.increase)

# -----------------------------------------------------------------------------------------------------------

# Now, for each of the unq_young_MuSC_interact and unq_aged_MuSC_interact dataframes,
# we want to check to see whether any of them have receptor that genes are actually differntially expressed.

# unq_R_DE_young_MuSC_interact <- unq_young_MuSC_interact[unq_young_MuSC_interact$geneB %in% young.increase$gene,]
# unq_R_DE_aged_MuSC_interact <- unq_aged_MuSC_interact[unq_aged_MuSC_interact$geneB %in% aged.increase$gene,]

# # To save these results, we replace the cluster number with their ids
# # young
# unq_R_DE_young_MuSC_interact2 <- unq_R_DE_young_MuSC_interact
# if (nrow(unq_R_DE_young_MuSC_interact2) > 0) {
	# unq_R_DE_young_MuSC_interact2 <- unq_R_DE_young_MuSC_interact
	# unq_R_DE_young_MuSC_interact2[unq_R_DE_young_MuSC_interact2$cellTypeB == 2, "cellTypeB"] <- "MuSC"

	# unq_R_DE_young_MuSC_interact2[unq_R_DE_young_MuSC_interact2$cellTypeA == 0, "cellTypeA"] <- "FAP_1"
	# unq_R_DE_young_MuSC_interact2[unq_R_DE_young_MuSC_interact2$cellTypeA == 1, "cellTypeA"] <- "EC_1"
	# unq_R_DE_young_MuSC_interact2[unq_R_DE_young_MuSC_interact2$cellTypeA == 2, "cellTypeA"] <- "MuSC"
	# unq_R_DE_young_MuSC_interact2[unq_R_DE_young_MuSC_interact2$cellTypeA == 3, "cellTypeA"] <- "FAP_2"
	# unq_R_DE_young_MuSC_interact2[unq_R_DE_young_MuSC_interact2$cellTypeA == 4, "cellTypeA"] <- "EC_2"
	# unq_R_DE_young_MuSC_interact2[unq_R_DE_young_MuSC_interact2$cellTypeA == 5, "cellTypeA"] <- "FAP_3"
	# unq_R_DE_young_MuSC_interact2[unq_R_DE_young_MuSC_interact2$cellTypeA == 6, "cellTypeA"] <- "Tenocyte"
	# unq_R_DE_young_MuSC_interact2[unq_R_DE_young_MuSC_interact2$cellTypeA == 7, "cellTypeA"] <- "FAP_4"
	# unq_R_DE_young_MuSC_interact2[unq_R_DE_young_MuSC_interact2$cellTypeA == 8, "cellTypeA"] <- "Schwann"
# }
# write.csv(unq_R_DE_young_MuSC_interact2, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/unique_R_DE_young_MuSC_interact.csv")
# rm(unq_R_DE_young_MuSC_interact2)

# # aged
# unq_R_DE_aged_MuSC_interact2 <- unq_R_DE_aged_MuSC_interact
# if (nrow(unq_R_DE_aged_MuSC_interact2) > 0) {
	# unq_R_DE_aged_MuSC_interact2[unq_R_DE_aged_MuSC_interact2$cellTypeB == 6, "cellTypeB"] <- "MuSC"

	# unq_R_DE_aged_MuSC_interact2[unq_R_DE_aged_MuSC_interact2$cellTypeA == 0, "cellTypeA"] <- "Myo_1"
	# unq_R_DE_aged_MuSC_interact2[unq_R_DE_aged_MuSC_interact2$cellTypeA == 1, "cellTypeA"] <- "FAP_1"
	# unq_R_DE_aged_MuSC_interact2[unq_R_DE_aged_MuSC_interact2$cellTypeA == 2, "cellTypeA"] <- "Myo_2"
	# unq_R_DE_aged_MuSC_interact2[unq_R_DE_aged_MuSC_interact2$cellTypeA == 3, "cellTypeA"] <- "EC"
	# unq_R_DE_aged_MuSC_interact2[unq_R_DE_aged_MuSC_interact2$cellTypeA == 4, "cellTypeA"] <- "FAP_2"
	# unq_R_DE_aged_MuSC_interact2[unq_R_DE_aged_MuSC_interact2$cellTypeA == 5, "cellTypeA"] <- "Tenocyte"
	# unq_R_DE_aged_MuSC_interact2[unq_R_DE_aged_MuSC_interact2$cellTypeA == 6, "cellTypeA"] <- "MuSC"
# }
# write.csv(unq_R_DE_aged_MuSC_interact2, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/unique_R_DE_aged_MuSC_interact.csv")
# rm(unq_R_DE_aged_MuSC_interact2)

# # -----------------------------------------------------------------------------------------------------------

# # We also want to check for simply DE receptors in the total interactome of the MuSCs (i.e. they don't have 
# # to be unique.

# R_DE_young_MuSC_interact <- young_MuSC_interact[young_MuSC_interact$geneB %in% young.increase$gene,]
# R_DE_aged_MuSC_interact <- aged_MuSC_interact[aged_MuSC_interact$geneB %in% aged.increase$gene,]

# R_DE_young_MuSC_interact2 <- R_DE_young_MuSC_interact
# if (nrow(R_DE_young_MuSC_interact2) > 0) {
	# R_DE_young_MuSC_interact2 <- R_DE_young_MuSC_interact
	# R_DE_young_MuSC_interact2[R_DE_young_MuSC_interact2$cellTypeB == 2, "cellTypeB"] <- "MuSC"

	# R_DE_young_MuSC_interact2[R_DE_young_MuSC_interact2$cellTypeA == 0, "cellTypeA"] <- "FAP_1"
	# R_DE_young_MuSC_interact2[R_DE_young_MuSC_interact2$cellTypeA == 1, "cellTypeA"] <- "EC_1"
	# R_DE_young_MuSC_interact2[R_DE_young_MuSC_interact2$cellTypeA == 2, "cellTypeA"] <- "MuSC"
	# R_DE_young_MuSC_interact2[R_DE_young_MuSC_interact2$cellTypeA == 3, "cellTypeA"] <- "FAP_2"
	# R_DE_young_MuSC_interact2[R_DE_young_MuSC_interact2$cellTypeA == 4, "cellTypeA"] <- "EC_2"
	# R_DE_young_MuSC_interact2[R_DE_young_MuSC_interact2$cellTypeA == 5, "cellTypeA"] <- "FAP_3"
	# R_DE_young_MuSC_interact2[R_DE_young_MuSC_interact2$cellTypeA == 6, "cellTypeA"] <- "Tenocyte"
	# R_DE_young_MuSC_interact2[R_DE_young_MuSC_interact2$cellTypeA == 7, "cellTypeA"] <- "FAP_4"
	# R_DE_young_MuSC_interact2[R_DE_young_MuSC_interact2$cellTypeA == 8, "cellTypeA"] <- "Schwann"
# }
# write.csv(R_DE_young_MuSC_interact2, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/R_DE_young_MuSC_interact.csv")
# rm(R_DE_young_MuSC_interact2)

# # aged
# R_DE_aged_MuSC_interact2 <- R_DE_aged_MuSC_interact
# if (nrow(R_DE_aged_MuSC_interact2) > 0) {
	# R_DE_aged_MuSC_interact2[R_DE_aged_MuSC_interact2$cellTypeB == 6, "cellTypeB"] <- "MuSC"

	# R_DE_aged_MuSC_interact2[R_DE_aged_MuSC_interact2$cellTypeA == 0, "cellTypeA"] <- "Myo_1"
	# R_DE_aged_MuSC_interact2[R_DE_aged_MuSC_interact2$cellTypeA == 1, "cellTypeA"] <- "FAP_1"
	# R_DE_aged_MuSC_interact2[R_DE_aged_MuSC_interact2$cellTypeA == 2, "cellTypeA"] <- "Myo_2"
	# R_DE_aged_MuSC_interact2[R_DE_aged_MuSC_interact2$cellTypeA == 3, "cellTypeA"] <- "EC"
	# R_DE_aged_MuSC_interact2[R_DE_aged_MuSC_interact2$cellTypeA == 4, "cellTypeA"] <- "FAP_2"
	# R_DE_aged_MuSC_interact2[R_DE_aged_MuSC_interact2$cellTypeA == 5, "cellTypeA"] <- "Tenocyte"
	# R_DE_aged_MuSC_interact2[R_DE_aged_MuSC_interact2$cellTypeA == 6, "cellTypeA"] <- "MuSC"
# }
# write.csv(R_DE_aged_MuSC_interact2, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Interactome/R_DE_aged_MuSC_interact.csv")
# rm(R_DE_aged_MuSC_interact2)