# Generates heatmap with 1 key gene per cluster (Fig. 1c of Skelly et al)

# THIS script REQUIRES THE SEURAT object

# ----------------------------------------------------------------------------

# first we determine the genes we are interested in:
genes <- c("Cd79a", "Cdh5", "Pecam1", "Pax7", "Col1a1", "Ly6a", "Pdgfra", "Plp1", "Tnmd")

# Next we get the average expression per cluster for the genes of interest:
avg_exp <- AverageExpression(young, genes.use = genes, show.progress = F)

# Next we reorder the columns in the avg_exp dataframe to be in the order that we want
avg_exp <- avg_exp[c(7, 2, 6, 3, 1, 4, 5, 9, 8)]

# Next we rename the cluster numbers to the cell types
names(avg_exp) <- c("B_Cells", "EC_1", "EC_2", "MuSCs", "FAP_1", "FAP_2", "FAP_3", "Schwann", "Tenocytes")

# Next we normalize the values in each row (gene) between 0 and 1.
avg_exp_norm <- t(apply(avg_exp, 1, function(x)(x)/(max(x))))

# Next we draw the heatmap
if (!require("plotly")) {install.packages("plotly"); require(plotly)}
p <- plot_ly(x=colnames(avg_exp_norm), y=row.names(avg_exp_norm), z = avg_exp_norm, type = "heatmap", showscale = T)