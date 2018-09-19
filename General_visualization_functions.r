# This script contains general visualization functions used in different areas.
# source("C:/users/sadeg/Google Drive/scRNA/muscle_scRNA/General_visualization_functions.r")

# -----------------------------------------------------------------------------------------------

if (!require("igraph")) install.packages("igraph")
library(igraph)

# -----------------------------------------------------------------------------------------------

# This function is to be used with igraph. It modifies autocurve.edges() so that it works also if 
# the edges are in different directions. 
# source: https://stackoverflow.com/questions/16875547/using-igraph-how-to-force-curvature-when-arrows-point-in-opposite-directions

autocurve.edges2 <- function(graph, start = 0.1)
{
    cm <- count.multiple(graph)
    mut <-is.mutual(graph)  #are connections mutual?
    el <- apply(get.edgelist(graph, names = FALSE), 1, paste,
        collapse = ":")
    ord <- order(el)
    res <- numeric(length(ord))
    p <- 1
    while (p <= length(res)) {
        m <- cm[ord[p]]
        mut.obs <-mut[ord[p]] #are the connections mutual for this point?
        idx <- p:(p + m - 1)
        if (m == 1 & mut.obs==FALSE) { #no mutual conn = no curve
            r <- 0
        }
        else {
            r <- seq(-start, start, length = m)
        }
        res[ord[idx]] <- r
        p <- p + m
    }
    res
}

# ------------------------------------------------------------------------------------------------