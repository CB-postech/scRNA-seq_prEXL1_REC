# Load required packages
library(scran)
library(scater)
library(Seurat)
library(viridis)

# Cluster cells and compute size factors
clusters <- quickCluster(pEXL1_GFP_QC2)
pEXL1_GFP_QC2 <- computeSumFactors(pEXL1_GFP_QC2, clusters = clusters)

# Log-normalize counts
pEXL1_GFP_QC2.norm <- logNormCounts(pEXL1_GFP_QC2, pseudo_count = 1)

# Model gene variance
gene_variance <- modelGeneVar(pEXL1_GFP_QC2.norm)

# Plot mean-variance relationship
plot(gene_variance$mean, gene_variance$total,
     xlab = "mean_log-expression",
     ylab = "variance")
curve(metadata(gene_variance)$trend(x),
      col = "blue",
      add = T)

# Identify highly variable genes (HVGs)
hvg.norm <- getTopHVGs(gene_variance, fdr.threshold = 0.1, var.field = 'bio', var.threshold = 0.01)
length(hvg.norm)

# Convert to Seurat object
pEXL1_GFP_QC2.seurat <- as.Seurat(pEXL1_GFP_QC2.norm,
                         counts = "counts",
                         data = "logcounts")
VariableFeatures(pEXL1_GFP_QC2.seurat) <- hvg.norm

# Scale data and run PCA
pEXL1_GFP_QC2.seurat <- ScaleData(pEXL1_GFP_QC2.seurat, features = VariableFeatures(pEXL1_GFP_QC2.seurat))
pEXL1_GFP_QC2.seurat <- RunPCA(pEXL1_GFP_QC2.seurat, features = VariableFeatures(pEXL1_GFP_QC2.seurat))

# Determine the number of PCs to use
ElbowPlot(pEXL1_GFP_QC2.seurat, ndims = 50)
PCs <- 7

# Cluster cells
pEXL1_GFP_QC2.seurat <- FindNeighbors(pEXL1_GFP_QC2.seurat, dims = 1:PCs)
pEXL1_GFP_QC2.seurat <- FindClusters(pEXL1_GFP_QC2.seurat, resolution = 2)

# Run UMAP
pEXL1_GFP_QC2.seurat <- RunUMAP(pEXL1_GFP_QC2.seurat, dims = 1:PCs, min.dist = 0.15)

# Visualize clusters
DimPlot(pEXL1_GFP_QC2.seurat, label = F, cols = c('darkgoldenrod','darkgoldenrod1','darkolivegreen','darkolivegreen1','forestgreen',
                                         'gold','gold4','goldenrod','green4','khaki1',
                                         'khaki3','khaki4','lemonchiffon','lightsalmon','lightsalmon3',
                                         'olivedrab','olivedrab3','palegreen','palegreen3','peachpuff2',
                                         'seagreen','yellowgreen','tan2'), pt.size = 1.5) + coord_fixed()

# Load root atlas cell type marker genes
celltype_marker <- read.csv('D:/project/DGIST_REC/plant_celltype_marker.csv', header = T)

# Calculate SEC module score
pEXL1_GFP_QC2.seurat <- AddModuleScore(pEXL1_GFP_QC2.seurat, features = list(celltype_marker$gene[celltype_marker$celltype %in% c('SEC')]), name = 'SEC')

# Visualize SEC module score
FeaturePlot(pEXL1_GFP_QC2.seurat, features = 'SEC1', order = T, pt.size = 1.5, cols = rev(viridis::rocket(50))) + coord_fixed()

# Calculate guard cell module score
pEXL1_GFP_QC2.seurat <- AddModuleScore(pEXL1_GFP_QC2.seurat, features = list(celltype_marker$gene[celltype_marker$celltype %in% c('guard_cell')]), name = 'GC')

# Visualize guard cell module score
FeaturePlot(pEXL1_GFP_QC2.seurat, features = 'GC1', order = T, pt.size = 3, cols = rev(viridis::rocket(50))) + coord_fixed()

# Calculate phloem module score
pEXL1_GFP_QC2.seurat <- AddModuleScore(pEXL1_GFP_QC2.seurat, features = list(celltype_marker$gene[celltype_marker$celltype %in% c('phloem')]),
                              name = 'phloem')
# Visualize phloem module score
FeaturePlot(pEXL1_GFP_QC2.seurat, features = 'phloem1', order = T, pt.size = 5, cols = rev(viridis::rocket(50))) + coord_fixed()

# Calculate xylem module score
pEXL1_GFP_QC2.seurat <- AddModuleScore(pEXL1_GFP_QC2.seurat, features = list(celltype_marker$gene[celltype_marker$celltype %in% c('xylem')]),
                              name = 'xylem')
# Visualize xylem module score
FeaturePlot(pEXL1_GFP_QC2.seurat, features = 'xylem1', order = T, pt.size = 5, cols = rev(viridis::rocket(50))) + coord_fixed()

# Calculate procambium module score
pEXL1_GFP_QC2.seurat <- AddModuleScore(pEXL1_GFP_QC2.seurat, features = list(celltype_marker$gene[celltype_marker$celltype %in% c('Procambium')]),
                              name = 'Procambium')
# Visualize procambium module score
FeaturePlot(pEXL1_GFP_QC2.seurat, features = 'Procambium1', order = T, pt.size = 5, cols = rev(viridis::rocket(50))) + coord_fixed()

# Subset REC
REC <- pEXL1_GFP_QC2.seurat[,!pEXL1_GFP_QC2.seurat$seurat_clusters %in% c('5','8','19','22','21','10')]