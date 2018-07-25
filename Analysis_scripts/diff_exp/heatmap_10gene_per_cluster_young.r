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
young_diff_markers <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_diff_markers.rds")
young <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young.rds")

# Next we extract the list of the top 10 differentially expressed genes in each cluster.
top_10_genes <- diff_genes_per_cluster(young_diff_markers, n_genes = 10)

# Next we extract the average gene expression of each cluster from the dataset and
# set it up in a new dataframe with the clusters as columns and genes as the rows.

young_avg_exp <- data.frame()

for (i in c(7,2,6,3,1,4,5,9,8)) {
	
	temp <- AverageExpression(young, genes.use = top_10_genes[,i], show.progress = F)
	young_avg_exp <- rbind(young_avg_exp, temp)
	
}
young_avg_exp <- young_avg_exp[c(7, 2, 6, 3, 1, 4, 5, 9, 8)]
saveRDS(young_avg_exp, "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/heatmap_10gene_avg_exp_young.rds")

# Next we normalize the expression in each row so the max is set as 1.
young_avg_exp_norm <- data.frame(t(apply(young_avg_exp, 1, function(x)(x)/(max(x)))))
names(young_avg_exp_norm) <- c("B_cells", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "Schwann", "Tenocytes")

saveRDS(young_avg_exp_norm, "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/heatmap_10gene_avg_exp_norm_young.rds")

# Next we extract the percentage of cells expressing each of the genes in each cluster
young_pct_exp_all <- AverageDetectionRate(young)

young_pct_exp <- data.frame()

for (i in c(7,2,6,3,1,4,5,9,8)) {
	
	temp <- young_pct_exp_all[top_10_genes[,i],]
	young_pct_exp <- rbind(young_pct_exp, temp)
	
}
young_pct_exp <- young_pct_exp[c(7, 2, 6, 3, 1, 4, 5, 9, 8)]
names(young_pct_exp) <- c("B_cells", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "Schwann", "Tenocytes")

saveRDS(young_pct_exp, "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/heatmap_10gene_pct_exp_young.rds")

# ---------------------------------------------------------------------------

# Function to create heatmap.