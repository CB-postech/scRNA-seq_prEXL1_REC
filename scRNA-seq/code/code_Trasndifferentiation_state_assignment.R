# Load required packages
library(scater)
library(scran)
library(Seurat)
library(pheatmap)
library(RColorBrewer)

# Convert to SingleCellExperiment object
REC.sce <- as.SingleCellExperiment(REC)

# Cluster cells and compute size factors
clusters <- quickCluster(REC.sce)
REC.sce <- computeSumFactors(REC.sce, clusters = clusters)

# Log-normalize counts
REC.sce.norm <- logNormCounts(REC.sce, pseudo_count = 1)

# Model gene variance
gene_variance <- modelGeneVar(REC.sce.norm)

# Plot mean-variance relationship
plot(gene_variance$mean, gene_variance$total,
     xlab = "mean_log-expression",
     ylab = "variance")
curve(metadata(gene_variance)$trend(x),
      col = "blue",
      add = T)

# Identify highly variable genes (HVGs)
hvg.norm <- getTopHVGs(gene_variance, var.field = 'bio', var.threshold = 0.1)
length(hvg.norm)

# Convert to Seurat object
REC <- as.Seurat(REC.sce.norm, counts = "counts", data = "logcounts")
VariableFeatures(REC) <- hvg.norm

# Reset PCA and UMAP reductions
REC@reductions$PCA <- NULL
REC@reductions$UMAP <- NULL

# Scale data and run PCA
REC <- ScaleData(REC, features = VariableFeatures(REC))
REC <- RunPCA(REC, features = VariableFeatures(REC), npcs = 300)

# Determine the number of PCs to use
ElbowPlot(REC, ndims = 20)
PCs <- 200

# Cluster cells
REC <- FindNeighbors(REC, dims = 1:PCs)
REC <- FindClusters(REC, resolution = 0.4)

# Run UMAP
REC <- RunUMAP(REC, dims = 1:PCs)

# Visualize clusters
DimPlot(REC, label = F, cols = c('#E7BF7E','#362935','#D07046','#974238','#A4AB47'), pt.size = 1.3, reduction = 'palantir') + coord_fixed()

# Calculate average expression per cluster
REC_aver <- AverageExpression(REC, group.by = 'seurat_clusters')
REC_aver <- REC_aver$originalexp
REC_aver <- REC_aver[VariableFeatures(REC_aver),]

# Calculate correlation matrix
cormat <- round(cor(REC_aver),2)

# Visualize the correlation matrix as a heatmap
state_heatmap <- pheatmap(cormat, cutree_rows = 3, cutree_cols = 3, color = colorRampPalette(brewer.pal(9,'YlOrRd'))(100), fontsize = 10)

# Assign states based on clusters
REC$states <- as.character(REC$seurat_clusters)
REC$states[REC$states %in% c('1','2','3')] <- 'State 1'
REC$states[REC$states %in% c('0')] <- 'State 2'
REC$states[REC$states %in% c('4')] <- 'State 3'

# Visualize states
DimPlot(REC, group.by = 'states', cols = c('#2B7A0B','#7DCE13','#EAE509'), pt.size = 1.3) + coord_fixed()