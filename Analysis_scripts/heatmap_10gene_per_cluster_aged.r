# Generates heatmap with 1 key gene per cluster (Fig. 1c of Skelly et al)

# THIS REQUIRES THE Seurat object

# ----------------------------------------------------------------------------

# First we load the required source functions and libraries
if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
library(Seurat)
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")
source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization_functions.r")

# ---------------------------------------------------------------------------

# Next we load the differentially expressed genes for each cluster.
aged_diff_markers <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_diff_markers.rds")
aged <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged.rds")

# Next we extract the list of the top 10 differentially expressed genes in each cluster.
top_10_genes <- diff_genes_per_cluster(aged_diff_markers, n_genes = 10)

# Next we extract the average gene expression of each cluster from the dataset and
# set it up in a new dataframe with the clusters as columns and genes as the rows.

aged_avg_exp <- data.frame()

for (i in c(7,3,10,1,4,6,2,5,8,9)) {
	
	temp <- AverageExpression(aged, genes.use = top_10_genes[,i], show.progress = F)
	aged_avg_exp <- rbind(aged_avg_exp, temp)
	
}
aged_avg_exp <- aged_avg_exp[c(7,3,10,1,4,6,2,5,8,9)]
saveRDS(young_avg_exp, "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/heatmap_10gene_avg_exp_aged.rds")

# Next we normalize the expression in each row so the max is set as 1.
aged_avg_exp_norm <- data.frame(t(apply(aged_avg_exp, 1, function(x)(x)/(max(x)))))
names(aged_avg_exp_norm) <- c("Immune_cells", "EC", "MuSC", "Myogenic_1", "Myogenic_2", "Myogenic_3", "FAP_1", "FAP_2", "FAP_3", "Tenocytes")

saveRDS(aged_avg_exp_norm, "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/heatmap_10gene_avg_exp_norm_aged.rds")

# Next we extract the percentage of cells expressing each of the genes in each cluster
aged_pct_exp_all <- AverageDetectionRate(aged)

aged_pct_exp <- data.frame()

for (i in c(7,3,10,1,4,6,2,5,8,9)) {
	
	temp <- aged_pct_exp_all[top_10_genes[,i],]
	aged_pct_exp <- rbind(aged_pct_exp, temp)
	
}
aged_pct_exp <- aged_pct_exp[c(7,3,10,1,4,6,2,5,8,9)]
names(aged_pct_exp) <- c("Immune_cells", "EC", "MuSC", "Myogenic_1", "Myogenic_2", "Myogenic_3", "FAP_1", "FAP_2", "FAP_3", "Tenocytes")

saveRDS(aged_pct_exp, "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/heatmap_10gene_pct_exp_aged.rds")

# ---------------------------------------------------------------------------

# Function to create heatmap.