# Int_analysis_functions

# This file contains the functions used in the Int_analysis script.

#-----------------------------------------------------------------------------------------------------

# This function reads in 10X data from CellRanger and returns a Seurat object created from the data

Seurat_from_10Xfile <- function(directory, min_cells = 3, min_genes = 200, max_genes = 5000, max_nUMI = 20000, max_percent_mito = 0.25) {
	
	# Read in the data from expression dataset
	datasample.data <- Read10X(data.dir = directory)
	
	# Create Seurat object and filter out the data
	datasample <- CreateSeuratObject(raw.data = datasample.data, min.cells = min_cells, min.genes = min_genes)
	mito.genes <- grep(pattern = "mt-", x = rownames(x=datasample@data), value = TRUE)
	percent.mito <- Matrix::colSums(datasample@raw.data[mito.genes, ])/Matrix::colSums(datasample@raw.data)
	datasample <- AddMetaData(object = datasample, metadata = percent.mito, col.name = "percent.mito")
	datasample <- FilterCells(object = datasample, subset.names = c("nGene", "nUMI","percent.mito"), low.thresholds = c(min_genes, -Inf, -Inf), high.thresholds = c(max_genes, max_nUMI, max_percent_mito))
	
	# Next we return the Seurat object
	return(datasample)
}

# -----------------------------------------------------------------------------------------------------