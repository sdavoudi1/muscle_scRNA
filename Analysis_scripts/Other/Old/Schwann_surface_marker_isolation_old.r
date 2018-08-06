# This script was used to extract and put together information for Schwann cell isolation

# Used this site to determine the location of the proteins created by the genes:
# https://www.uniprot.org/uniprot/P16388

# --------------------------------------------------------------------------------------

# First we load the young Seurat object and the differentially expressed genes in each cluster.
young_diff_markers <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_diff_markers.rds")
young <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young.rds")

# -------------------------------------------------------------------------------------

gene_list <- c("Gpr37l1", "Cadm1")

# -------------------------------------------------------------------------------------

# Here we create the matrix with the average gene expression
young_schwann_avg_exp <- data.frame()
young_schwann_avg_exp <- rbind(young_schwann_avg_exp, AverageExpression(young, genes.use = gene_list, show.progress = F))
young_schwann_avg_exp <- young_schwann_avg_exp[c(9, 7, 2, 6, 3, 1, 4, 5, 8)]
names(young_schwann_avg_exp) <- c("Schwann", "B_cells", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "Tenocytes")

# Here we create the normalized expression matrix
young_schwann_avg_exp_norm <- data.frame(t(apply(young_schwann_avg_exp, 1, function(x)(x)/(max(x)))))
names(young_schwann_avg_exp_norm) <- c("Schwann", "B_cells", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "Tenocytes")

# Here we determine the percent expression in the clusters
young_pct_exp_all <- AverageDetectionRate(young)
young_schwann_pct_exp <- data.frame()
young_schwann_pct_exp <- rbind(young_schwann_pct_exp, young_pct_exp_all[gene_list,])
young_schwann_pct_exp <- young_schwann_pct_exp[c(9, 7, 2, 6, 3, 1, 4, 5, 8)]
names(young_schwann_pct_exp) <- c("Schwann", "B_cells", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "Tenocytes")