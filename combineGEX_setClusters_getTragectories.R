# Load libraries
library(Seurat)
library(SeuratDisk)
library(sctransform)
library(ggplot2)
library(glmGamPoi)
library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(grDevices)
library(RColorBrewer)

# Load data matrix files
CON1 <- Read10X_h5("D://single_cell/Williams/GEX_ATAC/TG3_CON1/outs/filtered_feature_bc_matrix.h5")
CON2 <- Read10X_h5("D://single_cell/Williams/GEX_ATAC/TG4_CON2/outs/filtered_feature_bc_matrix.h5")
DES1 <- Read10X_h5("D://single_cell/Williams/GEX_ATAC/TG1_DES3/outs/filtered_feature_bc_matrix.h5")
DES2 <- Read10X_h5("D://single_cell/Williams/GEX_ATAC/TG2_DES4/outs/filtered_feature_bc_matrix.h5")

# Create and combine Seurat objects and determine mito RNA content
ctl_org1 <- CreateSeuratObject(counts = CON1$"Gene Expression", assay = "RNA", project = "ctl1")
ctl_org2 <- CreateSeuratObject(counts = CON2$"Gene Expression", assay = "RNA", project = "ctl2")
des_org1 <- CreateSeuratObject(counts = DES1$"Gene Expression", assay = "RNA", project = "des1")
des_org2 <- CreateSeuratObject(counts = DES2$"Gene Expression", assay = "RNA", project = "des2")
allorg <- merge(ctl_org1, y=c(ctl_org2,des_org1,des_org2), add.cell.ids = c("ctl1","ctl2","des1","des2"), 
                 project = "organoid")
allorg <- PercentageFeatureSet(allorg, pattern = "^mt-", col.name = "percent.mt")

# SCT normalize data and remove cells with high mitochondrial RNA
allorg <- SCTransform(allorg, method = "glmGamPoi", vars.to.regress = "percent.mt")

# Run dimensional reduction and find clusters
allorg <- RunPCA(allorg, verbose = FALSE)
allorg <- RunUMAP(allorg, dims = 1:30, verbose = FALSE)
allorg <- FindNeighbors(allorg, dims = 1:30, verbose = FALSE)
allorg <- FindClusters(allorg, resolution = 0.1, verbose = FALSE)

# Create UMAP split by sample and save Seurat object
DimPlot(allorg, label = TRUE, split.by = "orig.ident")
SaveH5Seurat(allorg, filename = "allorg", overwrite = TRUE)

# Find variable features and make slingshot matrices
allorg <- FindVariableFeatures(allorg, nfeatures = 2000)
dimred <- allorg@reductions$umap@cell.embeddings
clustering <- allorg$SCT_snn_res.0.1
counts <- as.matrix(allorg@assays$SCT@counts[allorg@assays$SCT@var.features, ])

# Define cell lineages
lineages <- getLineages(data = dimred, clusterLabels = clustering)

par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = clustering, cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = clustering, cex = 0.5, pch = 16)
lines(lineages, lwd = 3, col = "black")







# Convert seurat object to single cell experiment and add slingshot data
sce <- as.SingleCellExperiment(allorg, assay = "SCT")
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')

# Use Kolmogorov-Smirnov test to find differencnes in cell distributions across samples
ks.test(slingPseudotime(sce)[colData(sce)$orig.ident == "ctl1", 1],
        slingPseudotime(sce)[colData(sce)$orig.ident == "des2", 1])

# Make slingshot plot of inferred cell trajectories
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

# Plot initial trajectories estimates
plot(reducedDims(sce)$UMAP, col = brewer.pal(9,'Set1')[sce$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')



icMat <- evaluateK(counts = counts, sds = sce, k = 3:7, nGenes = 100, plot = TRUE)

















# Fit negative binomial GAM
sce <- fitGAM(sce)

# test for dynamic expression
ATres <- associationTest(sce)
topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$GMM[pst.ord]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

# Get cell lineages
lineages <- getLineages(sce.nest = dimred, clusterLabels = clustering)

lineages





# Create heatmap across clusters of select genes 
wjgenes <- read.delim("Wendy_Gene_List.txt", header = TRUE)
neworg <- ScaleData(allorg, features = wjgenes$Gene)
DoHeatmap(neworg, features = wjgenes$Gene)

# Create Violin and feature plots of select genes
VlnPlot(allorg,features= c("Krt18"))
FeaturePlot(allorg,features= c("Krt14","Trp63"), split.by = "orig.ident")

# Find differentially expressed genes between cell clusters and write to file
mymarkers <- FindAllMarkers(allorg, min.pct = 0.25, min.diff.pct = 0.25)
write.table(mymarkers, "Cluster_markers.txt", row.names = TRUE, col.names = NA, sep = "\t")

# Retrieve metadata for all cells and save to file
organodata <- allorg@meta.data
write.table(organodata, "Organoid_meta_data.txt", row.names = TRUE, sep = "\t")

# Save h5 seurat object to file
SaveH5Seurat(allorg, filename = "allorg", overwrite = TRUE)
