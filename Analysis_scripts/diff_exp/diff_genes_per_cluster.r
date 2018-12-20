# Differentially expressed genes in each cluster from the young and aged dataset
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/diff_exp/diff_genes_per_cluster.r")

# ----------------------------------------------------------------------------------------------

# young <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young.rds")
# young_diff_markers <- FindAllMarkers(object = young, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
# saveRDS(young_diff_markers, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_diff_markers.rds")

# aged <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged.rds")
# aged_diff_markers <- FindAllMarkers(object = aged, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
# saveRDS(aged_diff_markers, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_diff_markers.rds")

# young_reduced <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_reduced.rds")
# young_reduced_diff_markers <- FindAllMarkers(object = young_reduced, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
# saveRDS(young_reduced_diff_markers, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_reduced_diff_markers.rds")

# ----------------------------------------------------------------------------------------------

# young_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune.rds")
# young_noimmune_diff_markers <- FindAllMarkers(object = young_noimmune, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
# saveRDS(young_noimmune_diff_markers, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_diff_markers.rds")

aged_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune.rds")
aged_noimmune_diff_markers <- FindAllMarkers(object = aged_noimmune, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
saveRDS(aged_noimmune_diff_markers, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_noimmune_diff_markers.rds")

# young_noimmune_labeled <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_labeled.rds")
# young_noimmune_labeled_diff_markers <- FindAllMarkers(object = young_noimmune_labeled, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
# saveRDS(young_noimmune_labeled_diff_markers, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_labeled_diff_markers.rds")

# ----------------------------------------------------------------------------------------------

# After merging FAP_1 and FAP_3

young_noimmune_v2 <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_v2.rds")
young_noimmune_v2_diff_markers <- FindAllMarkers(object = young_noimmune_v2, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
saveRDS(young_noimmune__v2_diff_markers, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_v2_diff_markers.rds")

young_noimmune_v2_labeled <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_v2_labeled.rds")
young_noimmune_v2_labeled_diff_markers <- FindAllMarkers(object = young_noimmune_v2_labeled, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
saveRDS(young_noimmune_v2_labeled_diff_markers, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune_v2_labeled_diff_markers.rds")