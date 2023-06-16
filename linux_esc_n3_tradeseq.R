# Load libraries
library(Seurat)
library(sctransform)
library(ggplot2)
library(glmGamPoi)
library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(grDevices)
library(RColorBrewer)
library(BiocParallel)

# Load data matrix files
D2 <- Read10X_h5("/ddn/gs1/home/delkerda/single_cell/Hu/D2/outs/filtered_feature_bc_matrix.h5")
D4 <- Read10X_h5("/ddn/gs1/home/delkerda/single_cell/Hu/D4/outs/filtered_feature_bc_matrix.h5")
Epi <- Read10X_h5("/ddn/gs1/home/delkerda/single_cell/Hu/Epi/outs/filtered_feature_bc_matrix.h5")

# Create and combine Seurat objects and determine mito RNA content
d2_esc <- CreateSeuratObject(counts = D2$"Gene Expression", assay = "RNA", project = "d2")
d4_esc <- CreateSeuratObject(counts = D4$"Gene Expression", assay = "RNA", project = "d4")
ep_esc <- CreateSeuratObject(counts = Epi$"Gene Expression", assay = "RNA", project = "epi")
allesc <- merge(ep_esc, y=c(d2_esc,d4_esc), add.cell.ids = c("epi","d2","d4"), 
                 project = "mouse_esc")
allesc <- PercentageFeatureSet(allesc, pattern = "^mt-", col.name = "percent.mt")

# Filter out low quality cells and SCT normalize data
subesc <- subset(allesc, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & 
                   nCount_RNA > 1500 & nCount_RNA < 20000 & percent.mt < 20)
subesc <- SCTransform(subesc, method = "glmGamPoi", vars.to.regress = "percent.mt")

# Run dimensional reduction and find clusters
subesc <- RunPCA(subesc, verbose = FALSE)
subesc <- RunUMAP(subesc, dims = 1:30, verbose = FALSE)
subesc <- FindNeighbors(subesc, dims = 1:30, verbose = FALSE)
subesc <- FindClusters(subesc, verbose = FALSE)

# Convert seurat object to single cell experiment and filter cells
sce <- as.SingleCellExperiment(subesc, assay = "SCT")
geneFilter <- apply(assays(sce)$counts,1,function(x){sum(x >= 3) >= 10})
sce <- sce[geneFilter, ]

# Add slingshot data and make plot of inferred cell trajectories
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP',
                 approx_points = 150, start.clus = c("7","4"))
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

# Plot initial trajectories estimates
plot(reducedDims(sce)$UMAP, col = brewer.pal(9,'Set1')[sce$seurat_clusters], 
     pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

# Construct smooth curves and order cells
curves <- getCurves(sce, approx_points = 150, thresh = 0.01, stretch = 0.8, 
                    allow.breaks = FALSE, shrink = 0.99)
curves <- as.SlingshotDataSet(curves)
plot(reducedDims(sce)$UMAP, col = brewer.pal(9,"Set1")[sce$seurat_clusters], 
     asp = 1, pch = 16)
lines(curves, lwd = 3, col = 'black')

# Fit negative binomial GAM for differential expression analysis with tradeseq
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 8
sce <- fitGAM(sce, nknots = 11, parallel = TRUE, BPPARAM = BPPARAM)

# Save single cell experiment to file
sce_new <- readRDS("ESC_lineages_filtered.rds")

# Test for dynamic expression and identify progenitor marker genes
ATres <- associationTest(sce)
startRes <- startVsEndTest(sce)
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, counts(sce), gene = sigGeneStart)
plotGeneCount(curves, counts(sce), gene = sigGeneStart)

# Identify genes with different expression patterns
patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat])
plotSmoothers(sce, counts(sce), gene = rownames(patternRes)[oPat][4])
plotGeneCount(curves, counts(sce), gene = rownames(patternRes)[oPat][4])

# Discover differentiated cell type markers
endRes <- diffEndTest(sce)
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[1]]
plotSmoothers(sce, counts(sce), sigGene)
plotGeneCount(curves, counts(sce), gene = sigGene)

# Find early drivers of differentiation
plotGeneCount(curves, counts(sce), clusters = clusters, models = sce)
earlyDERes <- earlyDETest(sce, knots = c(1, 2))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly])
plotSmoothers(sce, counts(sce), gene = rownames(earlyDERes)[oEarly][2])

# Save single cell experiment to file and lineage specific DEG lists
saveRDS(sce, "ESC_lineages_filtered.rds")
write.table(startRes, "Progenitor_genes_filtered.txt", sep = "\t", row.names = TRUE)
write.table(ATres, "Association_genes_filtered.txt", sep = "\t", row.names = TRUE)
write.table(patternRes, "Pattern_genes_filtered.txt", sep = "\t", row.names = TRUE)
write.table(earlyDERes, "Driver_genes_filtered.txt", sep = "\t", row.names = TRUE)
write.table(endRes, "Marker_genes_filtered.txt", sep = "\t", row.names = TRUE)

