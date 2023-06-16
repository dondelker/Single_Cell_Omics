# Load libraries
library(Seurat)
library(sctransform)
library(ggplot2)
library(glmGamPoi)

# Load data matrix files
day09 <- Read10X(data.dir = "C://Users/delkerda/Desktop/Single_cell/mouse_placenta/2022/Day9/filtered_feature_bc_matrix")
day10 <- Read10X(data.dir = "C://Users/delkerda/Desktop/Single_cell/mouse_placenta/2022/Day10/filtered_feature_bc_matrix")
day12 <- Read10X(data.dir = "C://Users/delkerda/Desktop/Single_cell/mouse_placenta/2022/Day12/filtered_feature_bc_matrix")
day14 <- Read10X(data.dir = "C://Users/delkerda/Desktop/Single_cell/mouse_placenta/2022/Day14/filtered_feature_bc_matrix")

# Create and combine Seurat objects and determine mito RNA content
plac09 <- CreateSeuratObject(counts = day09, project = "d09")
plac10 <- CreateSeuratObject(counts = day10, project = "d10")
plac12 <- CreateSeuratObject(counts = day12, project = "d12")
plac14 <- CreateSeuratObject(counts = day14, project = "d14")
allplac <- merge(plac14, y=c(plac09,plac10,plac12), add.cell.ids = c("d14","d09","d10","d12"), 
                 project = "placenta")
allplac <- PercentageFeatureSet(allplac, pattern = "^mt-", col.name = "percent.mt")

# Subset Seurat object by cell barcode
cells.use <- read.delim("final_barcodes2.txt", header = FALSE)
subset_plac <- subset(allplac, cells = cells.use$V1)

# SCT normalize data
subset_plac <- SCTransform(subset_plac, method = "glmGamPoi", vars.to.regress = "percent.mt")

# Run dimensional reduction and find clusters
subset_plac <- RunPCA(subset_plac, verbose = FALSE)
subset_plac <- RunUMAP(subset_plac, dims = 1:30, verbose = FALSE)
subset_plac <- FindNeighbors(subset_plac, dims = 1:30, verbose = FALSE)
subset_plac <- FindClusters(subset_plac, verbose = FALSE)

# Create UMAP
DimPlot(subset_plac, label = TRUE)

# Retrieve metadata for all cells and save to file
trophodata <- subset_plac@meta.data
write.table(trophodata, "final_trophoblast_meta_data2.txt", row.names = TRUE, sep = "\t")

# Find differentially expressed genes between cell clusters and write to file
markers <- FindAllMarkers(subset_plac, min.pct = 0.25, min.diff.pct = 0.25)
write.table(markers, "final_tropho_markers2.txt", row.names = TRUE, col.names = NA, sep = "\t")

# Add custom cluster annotations
cluster_anno <- c("SynTI","GC","JZP1","SpT","JZP1","LaTP","S-TGC","S-TGC precursor",
                     "SynTI precursor","SpT precursor","SynTII","SynTII precursor",
                     "LaTP2","SpT precursor","JZP2")

names(cluster_anno) <- levels(subset_plac)
subset_tropho <- RenameIdents(subset_plac, cluster_anno)
DimPlot(subset_tropho, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
