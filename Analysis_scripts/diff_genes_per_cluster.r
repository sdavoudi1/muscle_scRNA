# Differentially expressed genes in each cluster from the young and aged dataset


young <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young.rds")
young_diff_markers <- FindAllMarkers(object = young, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
saveRDS(young_diff_markers, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_diff_markers.rds")

aged <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged.rds")
aged_diff_markers <- FindAllMarkers(object = aged, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
saveRDS(aged_diff_markers, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/aged_diff_markers.rds")

young_reduced <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_reduced.rds")
young_reduced_diff_markers <- FindAllMarkers(object = young_reduced, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
saveRDS(young_reduced_diff_markers, file = "C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_reduced_diff_markers.rds")