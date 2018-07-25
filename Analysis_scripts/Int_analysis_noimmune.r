# Int_analysis_noimmune.R
# Integrated Analysis of young and aged muscle scRNA
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/Int_analysis_noimmune.r")

# This scripts takes the young and aged muscle scRNA data and combines them using the Seurat algorithm
# outlined in Butler et al, Nat Biotech, 2018.

# --------------------------------------------------------------------------------------------------------

# We load the required libraries.
if(!require(Seurat)) {
install.packages("Seurat"); require(Seurat)}

# Next we load the functions we need
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.R")

# Next we load the aged and young Seurat objects with no immune or Red Blood cells using the scripts below:
young_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_preanalysis.rds")
aged_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune_preanalysis.rds")

# Next we determine genes used for CCA. These genes will be genes that are highly variable in both datasets
g.1 <- head(rownames(young_noimmune@hvg.info), 1000)
g.2 <- head(rownames(aged_noimmune@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(young_noimmune@scale.data))
genes.use <- intersect(genes.use, rownames(aged_noimmune@scale.data))

# Next we run the canonical correlation analysis (CCA) and combine the 2 objects into 1
muscle_noimmune.combined <- RunCCA(young_noimmune, aged_noimmune, genes.use = genes.use, num.cc = 30)

# Next to choose the CCs we use downstream for analysis, we use the "MetageneBicorPlot" which examines 
# a measure of correlation strength for each CC and find that this statistic generally staurates after a
# reasonable number of CCs.
# For us it seems to be ~20.
p1 <- MetageneBicorPlot(muscle_noimmune.combined, grouping.var = "group", dims.eval = 1:30, display.progress = FALSE)

# Next we align the CCA subspaces so we can use them for clustering.
muscle_noimmune.combined <- AlignSubspace(muscle_noimmune.combined, reduction.type = "cca", grouping.var = "group", dims.align = 1:20)

# Next we run a single integrated analysis on the cells to cluster them and perform tSNE
muscle_noimmune.combined <- RunTSNE(muscle_noimmune.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
muscle_noimmune.combined <- FindClusters(muscle_noimmune.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)

# # Next we plot the tSNE based on both clusters and age group
p1 <- TSNEPlot(muscle_noimmune.combined, do.return = T, pt.size = 0.5, group.by = "group")
p2 <- TSNEPlot(muscle_noimmune.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)

# Next we save the results
saveRDS(muscle_noimmune.combined, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/muscle_noimmune_combined.rds")

# Next we look at the expression of certain markers for each cluster to determine their identity.
# FeaturePlot(object = muscle.combined, features.plot = c("Pax7", "Myod1", "Des", "Cd34", "Ly6a", "Pdgfra", "Pecam1", "Plp1", "Tnmd"),min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)

# # To plot the tSNE plots of the young and aged side by side use the following snippet of code
# p1 <- TSNEPlot(muscle.combined, do.return = T, pt.size = 0.5, cells.use = muscle.combined@cell.names[1:length(which(muscle.combined@meta.data$orig.ident == "young"))]
# p2 <- TSNEPlot(muscle.combined, do.return = T, pt.size = 0.5, cells.use = muscle.combined@cell.names[length(which(muscle.combined@meta.data$orig.ident == "young"))+1:length(muscle.combined@cell.names)])
# plot_grid(p1,p2)

# To rename the clusters to what they are:
# new.ident <- c("Myo", "FAP_1", "EC_1", "MuSC", "FAP_2", "EC_2", "FAP_3", "Tenocyte", "FAP_4", "Schwann")
# for (i in 0:9) {
	# muscle_noimmune.combined <- RenameIdent(object = muscle_noimmune.combined, old.ident.name = i, new.ident.name = new.ident[i+1])
# }
# saveRDS(muscle_noimmune.combined, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/muscle_noimmune_combined_labeled.rds")