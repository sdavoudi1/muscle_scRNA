# Analysis_functions

# This file contains the functions used in the various scripts.

#-----------------------------------------------------------------------------------------------------

# This function reads imported RNA Seq data and returns a Seurat object created from the data

# REQUIRES: library(Seurat)

Seurat_from_10Xfile <- function(object, min_cells = 1, min_genes = 200, max_genes = 5000, max_nUMI = 20000, max_percent_mito = 0.25) {
	
	# Create Seurat object and filter out the data
	datasample <- CreateSeuratObject(raw.data = object, min.cells = min_cells, min.genes = min_genes)
	mito.genes <- grep(pattern = "mt-", x = rownames(x=datasample@data), value = TRUE)
	percent.mito <- Matrix::colSums(datasample@raw.data[mito.genes, ])/Matrix::colSums(datasample@raw.data)
	datasample <- AddMetaData(object = datasample, metadata = percent.mito, col.name = "percent.mito")
	datasample <- FilterCells(object = datasample, subset.names = c("nGene", "nUMI","percent.mito"), low.thresholds = c(min_genes, -Inf, -Inf), high.thresholds = c(max_genes, max_nUMI, max_percent_mito))
	datasample <- NormalizeData(object = datasample, normalization.method = "LogNormalize", scale.factor = 10000)
	datasample <- ScaleData(object = datasample, vars.to.regress = c("nUMI"))
	
	# Next we find the variable genes in the dataset
	datasample <- FindVariableGenes(object = datasample, do.plot = F, display.progress = F)
		
	# Next we return the Seurat object
	return(datasample)
}

# -----------------------------------------------------------------------------------------------------

# This function calculates the % of cells that are above a certain threshold in a list.

pct.above <- function(x, threshold) {
	return(length(x = x[x > threshold]) / length(x = x))
}

# -----------------------------------------------------------------------------------------------------

# This function takes the Seurat object and a gene list, and determines the average expression and the
# percentage of the cells expressing the genes on the list for each cluster of cells.

gene_per_cluster_exp_pct <- function(object, gene.list) {

	library(dplyr)
	library(tidyr)

	# first we extract the information regarding the genes of interest from the Seurat object and store
	# it in a dataframe. Next we create 2 new categories: the cell identifier, and cluster.
	data.to.analyze <- data.frame(FetchData(object = object, vars.all = gene.list))
	data.to.analyze$cell <- rownames(x = data.to.analyze)
	data.to.analyze$id <- object@ident
	
	# Next, we use the gather function (from dplyr and tidyr) to reorganize the dataframe based on genes
	# and the value we care about is the expression of that gene in the cell. We dont want to categorize
	# based on cell and cluster (hence they are put in as "-".
	data.to.analyze %>% gather(key = gene.list, value = expression, -c(cell, id)) -> data.to.analyze
	
	# Next we use group_by to group the dataframe based on cluster (id) and gene.list so we get the average
	# expression of each gene in our list in each of the clusters and the percentage of the cells in that
	# cluster expressing the gene.
	data.to.analyze %>% group_by(id, gene.list) %>% summarize(avg.exp = ExpMean(x = expression),
		pct.exp = pct.above(x = expression, threshold = 0)) -> data.to.analyze
	
	return(data.to.analyze)
}