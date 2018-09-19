network_graph_2 <- function(x, dir_mode = "undirected") {

  # Install and use igraph package
  if (!require("igraph")) install.packages("igraph")
  library(igraph)

  # Import desire data file and make it a matrix
  m=as.matrix(x)
  
  # Create a network (undirected, weighted), and edge for the network
  net=graph.adjacency(m,mode=dir_mode,weighted=TRUE,diag=TRUE)
  
  colrs <- c("red", "orange", "yellow","green", "blue", "purple","gray40", "brown", "violet")
  V(net)$color <- colrs[V(net)]
  V(net)$size <- 23
  V(net)$label.color <- "black"
  edge.start <- ends(net, es=E(net), names=F)[,1]
  edge.col <- V(net)$color[edge.start]
  
  #Plot graph (scale down the weight by 200 so the edges are shown appropricately)
  plot.igraph(net,vertex.label=V(net)$name, edge.color=edge.col,edge.width=E(net)$weight/150)
  
  
  # Pick out the data row by row
  for(i in 1:nrow(x)) {
    interact <-  vector(mode="double", length=0)
    edge_value <-  vector(mode="double", length=0)
    row <- x[i,]
    
    # Store the name of the row  
    for (i in names(row)){
      interact <- c(interact, rownames(row))
      interact <- c(interact, i)
    }
    
    # Store the value of the row
    for (i in row){
      edge_value <- c(edge_value, i)
    }
    
    # Align the name of the data with the edge and graph 
    g2 <- graph(as.numeric(interact)+1, directed=F) 
    V(g2)$label.color <- "black"
    V(g2)$label <- colnames(row)
    V(g2)$size <- 25
    colrs <- c("red", "orange", "yellow","green", "blue", "purple","gray40", "brown", "violet")
    V(g2)$color <- colrs[V(net)]
    E(g2)$weight <- edge_value
    plot(g2,edge.width=E(g2)$weight/80)
  }
}