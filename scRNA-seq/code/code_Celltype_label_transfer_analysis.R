# Load required libraries
library(Seurat)
library(RColorBrewer)

# Read data
leaf <- readRDS('D:/project/DGIST_REC/public/stomatal_Leaf.rds')

# Find intersection of variable features between REC and leaf
hvg <- intersect(VariableFeatures(REC), VariableFeatures(leaf))

# Find transfer anchors between reference (leaf) and query (REC)
transfer.anchors <- FindTransferAnchors(
  reference = leaf,
  query = REC,
  features = hvg,
  reference.assay = "RNA",
  query.assay = "originalexp",
  reduction = "cca",
  npcs = 10
)

# Transfer cell type labels from reference to query
celltype.predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = leaf$celltype,
  weight.reduction = 'cca',
  dims = 1:5
)

# Add predicted cell type labels to REC metadata
REC <- AddMetaData(
  REC,
  metadata = celltype.predictions,
  col.name = paste0('leaf.celltype', colnames(celltype.predictions))
)

# Plot predicted cell types on the palantir reduction
DimPlot(
  REC,
  group.by = 'leaf.celltypepredicted.id',
  reduction = 'palantir',
  pt.size = 1.3,
  cols = brewer.pal(9, 'Paired')
) +
  coord_fixed()