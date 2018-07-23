# Simplified_young_cluster.R

# This script merges the sub-clusters in the young muscle analysis to make major cell types

#--------------------------------------------------------------------------------------------

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8)

new.cluster.ids <- c(0,1,2,0,0,1,3,4,5)

young_reduced@ident <- plyr::mapvalues(x = young_reduced@ident, from = current.cluster.ids, to = new.cluster.ids)