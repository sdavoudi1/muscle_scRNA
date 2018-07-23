# Int_analysis.R
# Integrated Analysis of young and aged muscle scRNA

# This scripts takes the young and aged muscle scRNA data and combines them using the Seurat algorithm
# outlined in Butler et al, Nat Biotech, 2018.

# --------------------------------------------------------------------------------------------------------

# We load the required libraries.
if(!require(Seurat)) {
install.packages("Seurat"); require(Seurat)}

# Next we load the functions we need
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.R")

# Next we create the aged and young Seurat objects using the scripts below:
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/Young_seurat_object.r")
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/Aged_seurat_object.r")

# Next we determine genes used for CCA. These genes will be genes that are highly variable in both datasets
g.1 <- head(rownames(young@hvg.info), 1000)
g.2 <- head(rownames(aged@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(young@scale.data))
genes.use <- intersect(genes.use, rownames(aged@scale.data))

# Next we run the canonical correlation analysis (CCA) and combine the 2 objects into 1
muscle.combined <- RunCCA(young, aged, genes.use = genes.use, num.cc = 30)

# Next to choose the CCs we use downstream for analysis, we use the "MetageneBicorPlot" which examines 
# a measure of correlation strength for each CC and find that this statistic generally staurates after a
# reasonable number of CCs.
# For us it seems to be ~20.
p1 <- MetageneBicorPlot(muscle.combined, grouping.var = "group", dims.eval = 1:30, display.progress = FALSE)

# Next we align the CCA subspaces so we can use them for clustering.
muscle.combined <- AlignSubspace(muscle.combined, reduction.type = "cca", grouping.var = "group", dims.align = 1:20)

# Next we run a single integrated analysis on the cells to cluster them and perform tSNE
muscle.combined <- RunTSNE(muscle.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
muscle.combined <- FindClusters(muscle.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)

# # Next we plot the tSNE based on both clusters and age group
# p1 <- TSNEPlot(muscle.combined, do.return = T, pt.size = 0.5, group.by = "group")
# p2 <- TSNEPlot(muscle.combined, do.label = T, do.return = T, pt.size = 0.5)
# plot_grid(p1, p2)

# Next we look at the expression of certain markers for each cluster to determine their identity.
# FeaturePlot(object = muscle.combined, features.plot = c("Pax7", "Myod1", "Cd34", "Ly6a", "Cd79a", "Cd79b", "Kdr", "Pecam1", "Cd68"),min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)

# # To plot the tSNE plots of the young and aged side by side use the following snippet of code
# p1 <- TSNEPlot(muscle.combined, do.return = T, pt.size = 0.5, cells.use = muscle.combined@cell.names[1:length(which(muscle.combined@meta.data$orig.ident == "young"))]
# p2 <- TSNEPlot(muscle.combined, do.return = T, pt.size = 0.5, cells.use = muscle.combined@cell.names[length(which(muscle.combined@meta.data$orig.ident == "young"))+1:length(muscle.combined@cell.names)])
# plot_grid(p1,p2)