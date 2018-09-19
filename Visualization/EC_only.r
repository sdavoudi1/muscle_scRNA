# EC_only.r

# Visualizing EC cells

#--------------------------------------------------------------------------

young_noimmune <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/young_noimmune.rds")

young_EC <- SubsetData(object = young_noimmune, subset.name = "res.0.6", accept.value = c(1,4), do.clean = F)
TSNEPlot(object = young_EC) # 350 x 275

young_FAP <- SubsetData(object = young_noimmune, subset.name = "res.0.6", accept.value = c(0,3,5,7), do.clean = F)
TSNEPlot(object = young_FAP) # 475 x 375