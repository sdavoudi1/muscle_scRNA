network_graph <- function(x) {

  # Install and use igraph package
  if (!require("igraph")) install.packages("igraph")
  library(igraph)
  
  # Import desire data file and make it a matrix
  m=as.matrix(x)
  
  # Create a network (undirected, weighted), and edge for the network
  net=graph.adjacency(m,mode="undirected",weighted=TRUE,diag=FALSE)
  
  #Plot graph (scale down the weight by 200 so the edges are shown appropricately)
  plot.igraph(net,vertex.label=V(net)$name, edge.color="black",edge.width=E(net)$weight/120)
  
  
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
    E(g2)$weight <- edge_value
    V(g2)$label <- colnames(row)
    plot(g2,edge.width=E(g2)$weight/80)
  }
}