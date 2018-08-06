# Schwann_cell_Alaura
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/Other/Schwann_cell_Alaura.r")

# This file contains the scripts I used to prepare the files I sent to Alaura for the Schwann cell markers.
# Results were saved here: C:\Users\sadeg\Google Drive\scRNA\Results\Schwann Cells Alaura

# -----------------------------------------------------------------------------------------------------------

library(Seurat)
young_diff_markers <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_diff_markers.rds")
young <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young.rds")

# -----------------------------------------------------------------------------------------------------------

# Next we find the differential genes expressed in the Schwann cell cluster
Schwann_diff_genes <- young_diff_markers$gene[which(young_diff_markers$cluster == 8)]

# Next we find the average expression of the genes in all clusters.

young_schwann_avg_exp <- data.frame()
young_schwann_avg_exp <- rbind(young_schwann_avg_exp, AverageExpression(young, genes.use = Schwann_diff_genes, show.progress = F))
young_schwann_avg_exp <- young_schwann_avg_exp[c(9, 7, 2, 6, 3, 1, 4, 5, 8)]
names(young_schwann_avg_exp) <- c("Schwann", "B_cells", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "Tenocytes")
write.csv(young_schwann_avg_exp, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Schwann Cells Alaura/Young_Schwann_Average_expression.csv")

# Here we create the normalized expression matrix
young_schwann_avg_exp_norm <- data.frame(t(apply(young_schwann_avg_exp, 1, function(x)(x)/(max(x)))))
names(young_schwann_avg_exp_norm) <- c("Schwann", "B_cells", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "Tenocytes")
write.csv(young_schwann_avg_exp_norm, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Schwann Cells Alaura/Young_Schwann_Average_expression_normalized.csv")

# Next we find the percentagte expression of the genes in all clusters.
young_pct_exp_all <- AverageDetectionRate(young)
young_schwann_pct_exp <- data.frame()
young_schwann_pct_exp <- rbind(young_schwann_pct_exp, young_pct_exp_all[Schwann_diff_genes,])
young_schwann_pct_exp <- young_schwann_pct_exp[c(9, 7, 2, 6, 3, 1, 4, 5, 8)]
names(young_schwann_pct_exp) <- c("Schwann", "B_cells", "EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "Tenocytes")
write.csv(young_schwann_pct_exp, file = "C:/Users/sadeg/Google Drive/scRNA/Results/Schwann Cells Alaura/Young_Schwann_percent_expression.csv")