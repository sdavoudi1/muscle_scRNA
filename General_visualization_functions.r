# This script contains general visualization functions used in different areas.
# source("C:/users/sadeg/Google Drive/scRNA/muscle_scRNA/General_visualization_functions.r")

# List of functions:
#	- autocurve.edges2
#	- LabelPoint
#	- LabelUR
#	- LabelUL
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

# These few functions are created in Seurat page to help with drawing plots comparing
# aged and young clusters.

LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
    adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
    for (i in genes) {
        x1 <- exp.mat[i, 1]
        y1 <- exp.mat[i, 2]
        plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
            label = i, size = text.size)
        plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 + 
            adj.y.s, yend = y1, size = segment.size)
    }
    return(plot)
}

LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05, 
    adj.r.s = 0.05, ...) {
    return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t, 
        adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
    adj.l.s = 0.05, ...) {
    return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t, 
        adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}

LabelBR <- function(plot, genes, exp.mat, adj.b.t = 0.1, adj.r.t = 0.15, adj.b.s = 0.05, 
    adj.r.s = 0.05, ...) {
    return(LabelPoint(plot, genes, exp.mat, adj.y.t = -adj.b.t, adj.x.t = adj.r.t, 
        adj.y.s = -adj.b.s, adj.x.s = adj.r.s, ...))
}

LabelBL <- function(plot, genes, exp.mat, adj.b.t = 0.1, adj.l.t = 0.15, adj.b.s = 0.05, 
    adj.l.s = 0.05, ...) {
    return(LabelPoint(plot, genes, exp.mat, adj.y.t = -adj.b.t, adj.x.t = -adj.l.t, 
        adj.y.s = -adj.b.s, adj.x.s = -adj.l.s, ...))
}