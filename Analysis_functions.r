# Analysis_functions
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")

# This file contains the functions used in the various scripts.

#-----------------------------------------------------------------------------------------------------

# This function reads imported RNA Seq data and returns a Seurat object created from the data

Seurat_from_10Xfile <- function(object, min_cells = 1, min_genes = 200, max_genes = 5000, max_nUMI = 20000, max_percent_mito = 0.25) {
	
	# Make sure packages exist
	if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
	library(Seurat)
	
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

# This function is similar to the Seurat_from_10Xfile function with the exception that it only performs the
# normalization, scaling, and finding variable genes steps. The input is a Seurat object that we want
# to work on.

normalize_scale_fvg <- function(datasample) {

	# Make sure packages exist
	if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
	library(Seurat)

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

	if(!require(dplyr)) {
		install.packages("dplyr"); require(dplyr)}
	library(dplyr)

	if(!require(tidyr)) {
		install.packages("tidyr"); require(tidyr)}
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

# -----------------------------------------------------------------------------------------------------

# This function determines and returns the differentially expressed genes in each cluster and 
# exports the list of those expressed genes for each cluster of a Seurat object into separate
# txt files which can be input into for example gprofiler.

# THIS FUNCTION USED TO BE diff_genes_per_cluster

diff_genes_per_cluster_txt <- function(dataset_name, dataset, output_dir = "") {
	
	for (i in 1:length(unique(dataset$cluster))) {
	
		# Next we run gprofiler analysis on the select gene list and save the results as a text file.
		write.table(dataset$gene[which(dataset$cluster == i-1)], file = paste(paste(output_dir, paste(dataset_name,"cluster", sep="_"), sep = ""), i-1, "diff_gene_list.txt", sep = "_"), sep = "\t",row.names = F, col.names = F, quote = F)
	}
	
}

# -----------------------------------------------------------------------------------------------------

# This function determines and returns the top N number of differentially expressed genes in each cluster and 
# return them as a dataframe

diff_genes_per_cluster <- function(dataset, n_genes = 10) {
	
	df <- data.frame(matrix(0, ncol = length(unique(dataset$cluster)), nrow = n_genes))
	names(df) <- as.character(seq(0,length(unique(dataset$cluster))-1))
	
	for (i in 1:length(unique(dataset$cluster))) {
	
		# Next we run gprofiler analysis on the select gene list and save the results as a text file.
		gene_list <- dataset$gene[which(dataset$cluster == i-1)]
		df[,i] <- gene_list[1:n_genes]
	}
	
	return(df)
	
}

# -----------------------------------------------------------------------------------------------------

# This function generates gprofiler enrichment data. 
# THIS FUNCTION DOESNT WORK PROPERLY YET. ONCE INPUT INTO CYTOSCAPE, IT GIVES ERRORS.

if(!require(gProfileR)) {
	install.packages("gProfileR"); require(gProfileR)}
library(gProfileR)

runGprofiler <- function(genes, current_organism = "mmusculus", set_size_max = 500, set_size_min = 3, filter_min_overlap_size = 5 , exclude_iea = T, query_ordered = F){
  
    # Run gprofiler
    gprofiler_results <- gprofiler(query = genes ,organism = current_organism, ordered_query = query_ordered, significant = T,
								exclude_iea = exclude_iea, max_p_value = 1, min_set_size = set_size_min, max_set_size = set_size_max,
								correction_method = "fdr", src_filter = c("GO:BP","REAC"))
  
    #filter results so we exclude GO terms smaller than 3 genes and GO terms which our genes and the genes in it 
    # have less than filter_min_overlap_size in common.
    gprofiler_results <- gprofiler_results[which(gprofiler_results[,'term.size'] >= 3
                                        & gprofiler_results[,'overlap.size'] >= filter_min_overlap_size ),]
  
    # gProfileR returns corrected p-values only.  Set p-value to corrected p-value
    if(dim(gprofiler_results)[1] > 0){
    em_results <- cbind(gprofiler_results[,c("term.id","term.name","p.value","p.value")], 1,
                                gprofiler_results[,"intersection"])
    colnames(em_results) <- c("Name","Description", "pvalue","qvalue","phenotype","genes")
  
    return(em_results)
    } else {
      return("no gprofiler results for supplied query")
    }
}

# -----------------------------------------------------------------------------------------------------

# This function automatically creates an enrichment map in cytoscape. Make sure cytoscape is open.

# THIS FUNCTION DOES NOT WORK. NEEDS FIXING.

if(!require(RCy3)) {
	source("https://bioconductor.org/biocLite.R")
	biocLite("RCy3"); require(RCy3)}
library(RCy3)

create_em <- function(pvalue_threshold = 0.01, qvalue_threshold = 0.01, similarity_metric = "COMBINED",
                      similarity_threshold = 0.375, results_filename, model_name, gmt_file){
  
    current_network_name <- paste(model_name,pvalue_threshold,qvalue_threshold,sep="_")
  
    em_command = paste('enrichmentmap build analysisType="generic" ', "gmtFile=",gmt_file,
                     'pvalue=',pvalue_threshold, 'qvalue=',qvalue_threshold,
                     'similaritycutoff=',similarity_threshold,
                     'coefficients=',similarity_metric,
                     'enrichmentsDataset1=',results_filename,
                     sep=" ")
  
    #enrichment map command will return the suid of newly created network.
    response <- commandsGET(em_command)
  
    current_network_suid <- 0
  
    #enrichment map command will return the suid of newly created network unless it Failed.  If it failed it will contain the word failed
    if(grepl(pattern="Failed", response)){
      paste(response)
    } else {
      current_network_suid <- response
    }
  
    response <- renameNetwork(current_network_name, network = as.numeric(current_network_suid))
  
    return(current_network_suid);
  
}

