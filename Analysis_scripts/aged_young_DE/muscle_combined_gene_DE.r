# muscle_combined_gene_DE.R
# source("C:/Users/sadeg/Google Drive/scRNA/muscle_scRNA/Analysis_scripts/aged_young_DE/muscle_combined_gene_DE.R")


# This script was used to look at changes in gene expression between the aged
# and young clusters.
# ----------------------------------------------------------------------------------

if (!require("ks")) {install.packages("ks"); require(ks)}
library(ks)

if(!require(Seurat)) {install.packages("Seurat"); require(Seurat)}
	library(Seurat)

# We load in the labeled muscle_combined_noimmune data.
muscle.combined <- readRDS("C:/Users/sadeg/Google Drive/scRNA/data/young_aged/objects/muscle_noimmune_combined_labeled.rds")
muscle.combined2 <- muscle.combined
# ----------------------------------------------------------------------------------

# To analyze 
muscle.combined2@meta.data$celltype.age <- paste0(muscle.combined2@ident, "_", muscle.combined2@meta.data$group)
muscle.combined2 <- StashIdent(muscle.combined2, save.name = "celltype")
muscle.combined2 <- SetAllIdent(muscle.combined2, id = "celltype.age")
aged.increase <- FindMarkers(muscle.combined2, ident.1 = "MuSC_aged", ident.2 = "MuSC_young", only.pos = T)
young.increase <- FindMarkers(muscle.combined2, ident.1 = "MuSC_young", ident.2 = "MuSC_aged", only.pos = T)

# ----------------------------------------------------------------------------------

MuSC <- SubsetData(muscle.combined, ident.use = "MuSC", subset.raw = T)
MuSC <- SetAllIdent(MuSC, id = "group")
avg.MuSC <- log1p(AverageExpression(MuSC, show.progress = FALSE))
avg.MuSC$gene <- rownames(avg.MuSC)
avg.MuSC$DE <- "n.s."
avg.MuSC$DE[avg.MuSC$gene %in% rownames(young.increase)] <- "DE.young"
avg.MuSC$DE[avg.MuSC$gene %in% rownames(aged.increase)] <- "DE.aged"
p1 <- ggplot(avg.MuSC, aes(young, aged)) + geom_point() + ggtitle("MuSC gene expression")
p2 <- ggplot(avg.MuSC, aes(aged, young, colour = DE)) + geom_point() + ggtitle("MuSC gene expression")

# Next, to plot the labels on genes that are visually different on the graph, we calculate the 
# absolute value difference of the log1p values of that gene for young and aged and take the ones that
# have a difference of more than 1.
avg.MuSC$Difference <- abs(avg.MuSC$aged[,] - avg.MuSC$young[,])
genes.to.label.aged <- rownames(avg.MuSC[(avg.MuSC$DE == "DE.aged") & (avg.MuSC$Difference > 1),])
genes.to.label.young <- rownames(avg.MuSC[(avg.MuSC$DE == "DE.young") & (avg.MuSC$Difference > 1),])

# Then, to avoid overlaps, we separate the aged values to be labeled Upper Right, Upper Left, 
# Bottom Right, or Bottom Left.
UR <- c("Fmod", "UR", "S100a6", "Atp2a1", "Cdkn1c", "Dcn", "Acta1", "Tpm1", "Tcap", "Igfbp6", "Apod")
UL <- c("Tnni2", "Mylpf", "Gm10076", "Clec3b", "Serpinf1", "Gm10020", "C3")
BR <- c("Myl1", "Prg4", "Cst3", "Tnnc2", "Ckm", "Eno3", "Pvalb")
BL <- c("Tpm2", "Tnnt3", "Comp", "Thbs4", "Kera")
p2 <- LabelUR(p2, genes = UR, avg.MuSC, adj.u.t = 0.3, adj.u.s = 0.23)
p2 <- LabelUL(p2, genes = UL, avg.MuSC, adj.u.t = 0.3, adj.u.s = 0.23)
p2 <- LabelBL(p2, genes = BL, avg.MuSC, adj.b.t = 0.3, adj.b.s = 0.23)
p2 <- LabelBR(p2, genes = BR, avg.MuSC, adj.b.t = 0.3, adj.b.s = 0.23)
p2 <- LabelUR(p2, genes = genes.to.label.young, avg.MuSC, adj.u.t = 0.3, adj.u.s = 0.23)


# ----------------------------------------------------------------------------------

# To plot the average expressions of genes in the young vs aged MuSCs and determine outliers, we
# use kernel density estimate to determine the outlier. However, this method doesn't work because of the
# concentration of genes in the low expression corner. As a result, almost half the cells will fall
# in the low density area and be flagged. A better approach would be for the deviation of the center
# or simply plotting these same data but using the DE between aged and young MuSCs.
# Inspired by: https://stats.stackexchange.com/questions/114214/finding-outliers-on-a-scatter-plot

# p <- cbind(avg.MuSC$young, avg.MuSC$aged)
# dens <- kde(p)
# n.levels <- 15
# colors <- gray(seq(1, 0, length.out=n.levels))
# plot(dens, display="filled.contour2", cont=seq(0, 100, length.out=n.levels),
     # col=colors, xlab="Young", ylab="Aged")
	 
# dens <- kde(p, eval.points = p)
# avg.MuSC$Density <- dens$estimate

# points(avg.MuSC$young, avg.MuSC$aged, pch=19, cex=sqrt(avg.MuSC$Density/8))

# m <- mean(avg.MuSC$Density)
# e <- subset(avg.MuSC, subset=(Density < m/20))
# points(e$young, e$aged, col="#00000080")