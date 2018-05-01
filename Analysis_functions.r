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