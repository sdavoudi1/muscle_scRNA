# This file contains all the functions that were used to create the images for the 
# scRNA Seq project.
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Visualization_functions.r")

# This file contains the functions used in the various scripts.
# List of functions:
# 	-square_heatmap
# 	-circle_heatmap
#   -circle_heatmap_genelist
#	-circle_heatmap_genelist_spec_cluster
#	-igraph_circle_network_full
# 	-igraph_circle_net_lig_to_cluster
#	-ligand_source_heatmap


# ----------------------------------------------------------------------------------------------

if (!require("plotly")) {install.packages("plotly"); require(plotly)}
library(plotly)

if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
library(Seurat)

if(!require(dplyr)) {install.packages("dplyr"); require(dplyr)}
library(dplyr)

if (!require("igraph")) install.packages("igraph")
library(igraph)
 
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(RColorBrewer)
 
if (!require("colorspace")) install.packages("colorspace")
library(colorspace)

if (!require("xlsx")) install.packages("xlsx")
library(xlsx)

source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_functions.r")
source("C:/users/sadeg/Google Drive/scRNA/muscle_scRNA/General_visualization_functions.r")

# ----------------------------------------------------------------------------------------------

# This function creates a heatmap, with the average expression of the gene list for each cluster in a 
# normalized fashion. It takes the dataset, the gene list of interest, and the order of the clusters
# and return a heatmap for all of the genes in the list which are normalized between 0 and 1.

# The cluster order can be either column number or name. If the dataset is unlabeled, the function takes
# the cluster names as well. However, the cluster names are added after reordering, so add cluster names in
# the final order.

square_heatmap <- function(dataset, genes = c(), cluster_order, cluster_names,
						xaxis_font_size = 15, yaxis_font_size = 15) {

	# We first get the average expression for each of the genes in gene list.
	avg_exp <- AverageExpression(dataset, genes.use = genes, show.progress = F)
	
	# Next we reorder the columns in the avg_exp dataframe to be in the order that we want and
	# name them if they don't have names yet.
	if (is.numeric(cluster_order) == T) {
		avg_exp <- avg_exp[cluster_order]
		if ((is.character(cluster_names) == T) & (length(cluster_names) > 0)) {
			names(avg_exp) <- cluster_names
		}
	}
	if (is.character(cluster_order) == T) {
		avg_exp <- avg_exp[, cluster_order]
	}
	
	# Next we normalize the values in each row (gene) between 0 and 1.
	avg_exp_norm <- t(apply(avg_exp, 1, function(x)(x)/(max(x))))
	
	# Next we draw the heatmap
	p <- plot_ly(
		x=colnames(avg_exp_norm),
		y=row.names(avg_exp_norm),
		z = avg_exp_norm,
		type = "heatmap", showscale = T) %>%
		layout(xaxis = list(tickfont = list(size = xaxis_font_size)), 
        yaxis = list(tickfont = list(size = yaxis_font_size)))
	
	return(p)
}

# ----------------------------------------------------------------------------------------------

# This function creates a heatmap with circles, in which the circle size depicts ratio of cells in cluster
# expressing that gene, and the color is a normalized expression of that gene.

# USE UNLABELED DATASET WITH THIS FUNCTION

# Inputs: dataset, differential markers for each cluster in dataset.

circle_heatmap <- function(dataset, diff_markers, image_name = "image.png", image_dpi = 600, 
								image_width = 30, image_height = 5, n_genes = 10, cluster_names, 
								cluster_order, chart_name = "") {
		
	# First, we figure out the top n_genes expressed in each cluster.
	top_genes <- diff_genes_per_cluster(diff_markers, n_genes = n_genes)
	
	# Next we extract the average gene expression of each cluster from the dataset and
	# set it up in a new dataframe with the clusters as columns and genes as the rows.
	avg_exp <- data.frame()
	for (i in cluster_order) {
		temp <- AverageExpression(dataset, genes.use = top_genes[,i], show.progress = F)
		avg_exp <- rbind(avg_exp, temp)
	}
	avg_exp <- avg_exp[cluster_order]
	
	# Next we normalize the expression in each row so the max is set as 1.
	avg_exp_norm <- data.frame(t(apply(avg_exp, 1, function(x)(x)/(max(x)))))
	
	# Next we rename the columns to the names we want.
	if (length(cluster_names) > 0) {
		names(avg_exp_norm) <- cluster_names
	}
	
	# Next we extract the percentage of cells expressing each of the genes in each cluster
	pct_exp_all <- AverageDetectionRate(dataset)
	pct_exp <- data.frame()
	for (i in cluster_order) {
		temp <- pct_exp_all[top_genes[,i],]
		pct_exp <- rbind(pct_exp, temp)
	}
	pct_exp <- pct_exp[cluster_order]
	names(pct_exp) <- cluster_names

	#---------------------- Next we begin drawing the heatmap with circles ---------------------
	
	# First, we set up the vectors (x & y labels)
	cell_type <- c(colnames(pct_exp))
	genes <- c(rownames(pct_exp))
	
	# Create the data frame (add pct_exp & avg_exp_norm to the dataframe)
	df <- expand.grid(genes, cell_type)
	per_expression <-  vector(mode="double", length=0)
	for (item in pct_exp){
		per_expression <- c(per_expression, item)
	}
	df$pct_exp <- per_expression

	avg_expression <-  vector(mode="double", length=0)
	for (item in avg_exp_norm){
		avg_expression <- c(avg_expression, item)
	}
	df$avg_exp_norm <- avg_expression
	
	#Setting up background (set color interval for every n_genes)
	rect_left <- vector(mode="integer", length=0)
	x_pos = n_genes
	i <- 1
	while (i <= (length(cell_type)/2)){
		rect_left <- c(rect_left, x_pos)
		x_pos <- x_pos + 2*n_genes
		i <- i + 1
	}
	
	rectangles <- data.frame(
		xmin = rect_left+0.5,
		xmax = rect_left+n_genes+0.5,
		ymin = -Inf,
		ymax = +Inf
	)
	
	g <- ggplot(df, aes(Var1, Var2)) + labs(title=chart_name) + geom_point(colour=NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 70, hjust = 1)) + xlab("Cell Type") + ylab("Genes")
	g + scale_size_continuous(range=c(0,5)) + scale_color_gradient(low = 'blue', high = 'red') + coord_fixed(ratio = 1.8) + geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill='grey', alpha=0.3,inherit.aes = FALSE) + geom_point(data = df, aes(x=Var1, y=Var2,color = avg_exp_norm, size = pct_exp), alpha = 1.0)
	
	# Save data with given height and width
	ggsave(file=image_name, width=image_width, height=image_height, dpi=image_dpi)
	
}

# ----------------------------------------------------------------------------------------------

# This function creates a heatmap with circles, in which the circle size depicts ratio of cells in cluster
# expressing that gene, and the color is a normalized expression of that gene.

# USE UNLABELED DATASET WITH THIS FUNCTION

# Inputs: dataset, genelist for which we want to draw the heatmap for.

circle_heatmap_genelist <- function(dataset, gene_list = "",image_name = "image.png", image_dpi = 600, 
								image_width = 30, image_height = 5, reorder_cluster = F, cluster_order = "", 
								rename_cluster = F, cluster_names = "", chart_name = "") {
		
	# We extract the average gene expression of each cluster from the dataset and
	# set it up in a new dataframe with the clusters as columns and genes as the rows.
	avg_exp_norm <- Average_gene_exp_per_cluster(dataset = dataset, gene_list = gene_list, 
						reorder_cluster = reorder_cluster, cluster_order = cluster_order, normalize = T,
						rename_cluster = rename_cluster, cluster_names = cluster_names)
	
	# Next we extract the percentage of cells expressing each of the genes in each cluster
	pct_exp <- pct_exp_per_cluster(dataset = dataset, gene_list = gene_list,
					reorder_cluster = reorder_cluster, cluster_order = cluster_order, 
					rename_cluster = rename_cluster, cluster_names = cluster_names)

	#---------------------- Next we begin drawing the heatmap with circles ---------------------
	
	# First, we set up the vectors (x & y labels)
	cell_type <- c(colnames(pct_exp))
	genes <- c(rownames(pct_exp))
	
	# Create the data frame (add pct_exp & avg_exp_norm to the dataframe)
	df <- expand.grid(genes, cell_type)
	per_expression <-  vector(mode="double", length=0)
	for (item in pct_exp){
		per_expression <- c(per_expression, item)
	}
	df$pct_exp <- per_expression

	avg_expression <-  vector(mode="double", length=0)
	for (item in avg_exp_norm){
		avg_expression <- c(avg_expression, item)
	}
	df$avg_exp_norm <- avg_expression
	
	#Setting up background (set color interval for every 1 gene)
	rect_left <- vector(mode="integer", length=0)
	x_pos = 1
	i <- 1
	while (i <= (length(cell_type)/2)){
		rect_left <- c(rect_left, x_pos)
		x_pos <- x_pos + 2
		i <- i + 1
	}
	
	rectangles <- data.frame(
		xmin = rect_left+0.5,
		xmax = rect_left+0.5,
		ymin = -Inf,
		ymax = +Inf
	)
	
	g <- ggplot(df, aes(Var1, Var2)) + labs(title=chart_name) + geom_point(colour=NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 70, hjust = 1)) + xlab("Cell Type") + ylab("Genes")
	g + scale_size_continuous(range=c(0,5)) + scale_color_gradient(low = 'blue', high = 'red') + coord_fixed(ratio = 1.8) + geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill='grey', alpha=0.3,inherit.aes = FALSE) + geom_point(data = df, aes(x=Var1, y=Var2,color = avg_exp_norm, size = pct_exp), alpha = 1.0)
	
	# Save data with given height and width
	ggsave(file=image_name, width=image_width, height=image_height, dpi=image_dpi)
	
}

# ----------------------------------------------------------------------------------------------

# This function takes a gene list, and cluster numbers and returns the circle heatmap for that gene list
# and those numbers.

# The clusters have to be adjusted to start from 1, not from 0.

circle_heatmap_genelist_spec_cluster <- function(dataset, gene_list = "", n_genes = 10,
								clusters = c(), rename_cluster = F, cluster_names = "",
								image_dpi = 600, image_width = 10, image_height = 5,chart_name = "",
								image_name = "image.png"){
	if (length(clusters) > 0) {
		
		# We extract the average gene expression of each cluster from the dataset and
		# set it up in a new dataframe with the clusters as columns and genes as the rows.
		avg_exp_norm <- Average_gene_exp_per_cluster(dataset = dataset, gene_list = gene_list, 
							normalize = T)
		avg_exp_norm <- avg_exp_norm[,clusters]	
		if (rename_cluster == T) {
			if (length(cluster_names) > 0) {
				names(avg_exp_norm) <- cluster_names
			}
		}
		
		# Next we extract the percentage of cells expressing each of the genes in each cluster
		pct_exp <- pct_exp_per_cluster(dataset = dataset, gene_list = gene_list)
		pct_exp <- pct_exp[,clusters]
		if (rename_cluster == T) {
			if (length(cluster_names) > 0) {
				names(pct_exp) <- cluster_names
			}
		}
		

		#---------------------- Next we begin drawing the heatmap with circles ---------------------
		
		# First, we set up the vectors (x & y labels)
		cell_type <- c(colnames(pct_exp))
		genes <- c(rownames(pct_exp))
		
		# Create the data frame (add pct_exp & avg_exp_norm to the dataframe)
		df <- expand.grid(genes, cell_type)
		per_expression <-  vector(mode="double", length=0)
		for (item in pct_exp){
			per_expression <- c(per_expression, item)
		}
		df$pct_exp <- per_expression

		avg_expression <-  vector(mode="double", length=0)
		for (item in avg_exp_norm){
			avg_expression <- c(avg_expression, item)
		}
		df$avg_exp_norm <- avg_expression
		
		
		#Setting up background (set color interval for every n_genes)
		rect_left <- vector(mode="integer", length=0)
		x_pos = n_genes
		i <- 1
		while (i <= (length(cell_type)/2)){
			rect_left <- c(rect_left, x_pos)
			x_pos <- x_pos + 2*n_genes
			i <- i + 1
		}
		
		rectangles <- data.frame(
			xmin = rect_left+0.5,
			xmax = rect_left+n_genes+0.5,
			ymin = -Inf,
			ymax = +Inf
		)	
				
		g <- ggplot(df, aes(Var1, Var2)) + labs(title=chart_name) + geom_point(colour=NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 70, hjust = 1)) + xlab("Cell Type") + ylab("Genes")
		g + scale_size_continuous(range=c(0,5)) + scale_color_gradient(low = 'blue', high = 'red') + coord_fixed(ratio = 1.8) + geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill='grey', alpha=0.3,inherit.aes = FALSE) + geom_point(data = df, aes(x=Var1, y=Var2,color = avg_exp_norm, size = pct_exp), alpha = 1.0)
		
		# Save data with given height and width
		ggsave(file=image_name, width=image_width, height=image_height, dpi=image_dpi)
	}
	else {
		cat("No clusters selected")
	}	
}

# ----------------------------------------------------------------------------------------------

# This function takes a dataframe version of the network data, converts it into a adjacency matrix
# and then creates a circular network graph.

igraph_circle_network_full <- function(x, dir_mode = "directed", show_plot = T,
									   width_multiplier = 7.5, start_curve = 0.1,
									   label_v = T,
									   save_pdf = F, plot_name = "test.pdf") {

	# Import desire data file and make it a matrix
	xm <- as.matrix(x)

	# Create adjacency matrix (directed, weighted), and edge for the network
	gm <- graph.adjacency(xm, mode = dir_mode, weighted = TRUE, diag = TRUE)

	# Set circular layout
	layout <- layout.circle(gm)

	# We next set the edge weights to max 5.
	edge_width <- width_multiplier*E(gm)$weight/max(E(gm)$weight)

	# Next we setup our color palette
	cols <- brewer.pal(12, name = "Set3")
	cols <- readhex(file = textConnection(paste(cols, collapse = "\n")), class = "RGB")
	cols <- as(cols, "HLS")
	cols@coords[, "L"] <- cols@coords[, "L"] * 0.85
	cols <- as(cols, "RGB")
	cols <- hex(cols)
	cols_order <- c(1,3,4,5,6,7,8,9,10,2,11,12)
	cols <- cols[cols_order]

	# Next we set the vertices and edge colors
	#colr <- c("red", "orange", "cyan","green", "blue", "purple","gray40", "brown", "violet")
	V(gm)$color <- cols[V(gm)]
	V(gm)$size <- 23
	V(gm)$label.color <- "black"
	edge.start <- ends(gm, es=E(gm), names=F)[,1]
	edge.col <- V(gm)$color[edge.start]
	
	# We also choose to show or omit the vertex labels
	if (label_v == T) {
		vertex_labels <- V(gm)$name
	}
	if (label_v == F) {
		vertex_labels <- NA
	}
	
	# Next, we derive the angles for the auto-vertices
	n_v <- length(V(gm))
	loop_angle <- numeric(n_v*n_v)
	a <- 1
	for (i in 1:n_v) {
		loop_angle[a] <- (2*pi)*(1 - ((i-1)/n_v))
		a <- a + (n_v + 1)
	}

	# Next, we plot the results
	curves <- autocurve.edges2(gm, start = start_curve)
	
	if (save_pdf == T) {
		pdf(plot_name, 10, 10)
		plot(gm, layout = layout, 
			edge.width = edge_width, edge.arrow.size = 0,
			edge.color = edge.col, edge.curved = curves, edge.loop.angle = loop_angle,
			label = label_v,
			vertex.label = vertex_labels, margin = 0.5)
		dev.off()
	}
	if (show_plot == T) {
		plot(gm, layout = layout, 
			edge.width = edge_width, edge.arrow.size = 0,
			edge.color = edge.col, edge.curved = curves, edge.loop.angle = loop_angle, 
			vertex.label = vertex_labels, margin = 0.5)
	}
}

# ----------------------------------------------------------------------------------------------

# This function takes a dataframe version of the network data, converts it into a adjacency matrix
# and then creates a circular network graph that only shows the ligands presented to each of the clusters.

igraph_circle_net_lig_to_cluster <- function(x, cluster = 2, order = c(1,4,2,0,3,6,5,7) + 1, 
									   dir_mode = "directed", show_plot = T,
									   width_multiplier = 7.5, start_curve = 0.1,
									   label_v = T,
									   save_pdf = F, plot_name = "test.pdf") {
									   
	x <- x[order,]
	x <- x[,order]
	cluster_real <- cluster+1
	x[,order[!(order %in% cluster_real)]] <- 0
	
	xm <- as.matrix(x)
	gm <- graph.adjacency(xm, mode = dir_mode, weighted = TRUE, diag = TRUE)
	
	layout <- layout.circle(gm)
	edge_width <- width_multiplier*E(gm)$weight/max(E(gm)$weight)
	
	cols <- brewer.pal(12, name = "Set3")
	cols <- readhex(file = textConnection(paste(cols, collapse = "\n")), class = "RGB")
	cols <- as(cols, "HLS")
	cols@coords[, "L"] <- cols@coords[, "L"] * 0.85
	cols <- as(cols, "RGB")
	cols <- hex(cols)
	cols_order <- c(1,3,4,5,6,7,8,9,10,2,11,12)
	cols <- cols[cols_order]
	
	V(gm)$color <- cols[V(gm)]
	V(gm)$size <- 23
	V(gm)$label.color <- "black"
	edge.start <- ends(gm, es=E(gm), names=F)[,1]
	edge.col <- V(gm)$color[edge.start]
	
	if (label_v == T) {
		vertex_labels <- V(gm)$name
	}
	if (label_v == F) {
		vertex_labels <- NA
	}	
	
	n_v <- length(V(gm))
	loop_angle <- numeric(n_v)
	loop_angle[cluster_real] <- (2*pi)*(1 - ((cluster_real-1)/n_v))
	
	curves <- autocurve.edges2(gm, start = start_curve)
	
	if (save_pdf == T) {
		pdf(plot_name, 10, 10)
		plot(gm, layout = layout, 
			edge.width = edge_width, edge.arrow.size = 0,
			edge.color = edge.col, edge.curved = curves, edge.loop.angle = loop_angle,
			label = label_v,
			vertex.label = vertex_labels, margin = 0.5)
		dev.off()
	}
	if (show_plot == T) {
		plot(gm, layout = layout, 
			edge.width = edge_width, edge.arrow.size = 0,
			edge.color = edge.col, edge.curved = curves, edge.loop.angle = loop_angle, 
			vertex.label = vertex_labels, margin = 0.5)
	}
}

# ----------------------------------------------------------------------------------------------

# This function creates a heatmap with the source of the ligands presented to a specific cluster. 
# It takes the young seurat object, the receiving interactome dataset for the cluster in mind
# (using the cluster_interactome_rec_v2 function), the order of the dataset, and cluster of interest 
# as input.

ligand_source_heatmap <- function(interac, labeled_dataset = "young", 
								cluster = "MuSC", cluster_order = c("EC_1", "EC_2", "MuSC", "FAP_1", "FAP_2", "FAP_3", "Tenocyte", "Schwann")
								) {

	# First, we change the names in the interactome dataset columns so they are the cell types instead of numbers
	if (nrow(interac) > 0) {
		interac[interac$cellTypeB == 0, "cellTypeB"] <- "FAP_1"
		interac[interac$cellTypeB == 1, "cellTypeB"] <- "EC_1"
		interac[interac$cellTypeB == 2, "cellTypeB"] <- "MuSC"
		interac[interac$cellTypeB == 3, "cellTypeB"] <- "FAP_2"
		interac[interac$cellTypeB == 4, "cellTypeB"] <- "EC_2"
		interac[interac$cellTypeB == 5, "cellTypeB"] <- "Tenocyte"
		interac[interac$cellTypeB == 6, "cellTypeB"] <- "FAP_3"
		interac[interac$cellTypeB == 7, "cellTypeB"] <- "Schwann"

		interac[interac$cellTypeA == 0, "cellTypeA"] <- "FAP_1"
		interac[interac$cellTypeA == 1, "cellTypeA"] <- "EC_1"
		interac[interac$cellTypeA == 2, "cellTypeA"] <- "MuSC"
		interac[interac$cellTypeA == 3, "cellTypeA"] <- "FAP_2"
		interac[interac$cellTypeA == 4, "cellTypeA"] <- "EC_2"
		interac[interac$cellTypeA == 5, "cellTypeA"] <- "Tenocyte"
		interac[interac$cellTypeA == 6, "cellTypeA"] <- "FAP_3"
		interac[interac$cellTypeA == 7, "cellTypeA"] <- "Schwann"
	}
	
	# Next we extract the unique ligands presented to the cluster of interest, and determine the average 
	# expression of each of the clusters for those genes, and reorganize the results in order of our cluster order.
	ligands <- unique(interac[,1])
	avg_exp <- AverageExpression(object = labeled_dataset, genes.use = ligands, show.progress = F)
	
	if (is.character(cluster_order) == T) {
		avg_exp <- avg_exp[, cluster_order]
	}
	
	# Next we normalize the values so the cluster with the highest average is set to 1, and identify
	# which cluster has the highest expression for each of the genes and store it in max_col.
	avg_exp_norm <- t(apply(avg_exp, 1, function(x)(x)/(max(x))))
	max_col <- apply(avg_exp_norm, 1, which.max)
	avg_exp_norm <- data.frame(avg_exp_norm)
	
	# Next we create empty temporary datasets to store information in.
	EC_1 <- data.frame(matrix(ncol = 8, nrow = 0))
	names(EC_1) <- cluster_order
	EC_2 <- EC_1
	MuSC <- EC_1
	FAP_1 <- EC_1
	FAP_2 <- EC_1
	FAP_3 <- EC_1
	Tenocyte <- EC_1
	Schwann <- EC_1
	
	# Next we assign each ligand normalized expression values, to the dataset_cluster which has
	# the highest normalized expression value. At the end, we merge them all in the order we want.
	for (i in (1:nrow(avg_exp_norm))) {
		if (max_col[i] == 1) {
			EC_1 <- rbind(EC_1, avg_exp_norm[i,])
		}
		if (max_col[i] == 2) {
			EC_2 <- rbind(EC_2, avg_exp_norm[i,])
		}
		if (max_col[i] == 3) {
			MuSC <- rbind(MuSC, avg_exp_norm[i,])
		}
		if (max_col[i] == 4) {
			FAP_1 <- rbind(FAP_1, avg_exp_norm[i,])
		}
		if (max_col[i] == 5) {
			FAP_2 <- rbind(FAP_2, avg_exp_norm[i,])
		}
		if (max_col[i] == 6) {
			FAP_3 <- rbind(FAP_3, avg_exp_norm[i,])
		}
		if (max_col[i] == 7) {
			Tenocyte <- rbind(Tenocyte, avg_exp_norm[i,])
		}
		if (max_col[i] == 8) {
			Schwann <- rbind(Schwann, avg_exp_norm[i,])
		}
	}
	avg_exp_norm_reorganized <- rbind(EC_1, EC_2, MuSC, FAP_1, FAP_2, FAP_3, Tenocyte, Schwann)
	avg_exp_norm_reorg_matrix <- t(as.matrix(avg_exp_norm_reorganized))
	
	# Finally, we plot the heatmap.
	xaxis_font_size = 15
	yaxis_font_size = 15
	p <- plot_ly(
			y=rownames(avg_exp_norm_reorg_matrix),
			x=colnames(avg_exp_norm_reorg_matrix),
			z = avg_exp_norm_reorg_matrix,
			type = "heatmap", showscale = T) %>%
			layout(xaxis = list(tickfont = list(size = xaxis_font_size)), 
			yaxis = list(tickfont = list(size = yaxis_font_size)))
	
	return(p)
}	