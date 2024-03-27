# Load required libraries
library(Seurat)
library(Hmisc)
library(reshape2)
library(ggplot2)
library(viridis)
library(RColorBrewer)

# Read in the data
leaf <- readRDS('D:/project/DGIST_REC/public/stomatal_Leaf.rds')

# Calculate average expression for each dataset
REC_aver <- AverageExpression(REC, group.by = 'states')
leaf_aver <- AverageExpression(leaf, group.by = 'celltype')

# Extract the average expression matrices
REC_aver <- REC_aver$originalexp
leaf_aver <- leaf_aver$RNA

# Find highly variable genes (HVGs) that are common between REC and leaf
hvg <- intersect(VariableFeatures(REC), VariableFeatures(leaf))

# Subset the average expression matrices to only include HVGs
REC_aver <- as.matrix(REC_aver)[rownames(REC_aver) %in% hvg,]
leaf_aver <- as.matrix(leaf_aver)[rownames(leaf_aver) %in% hvg,]

# Check dimensions of the matrices
dim(leaf_aver)
dim(REC_aver)

# Transpose the matrices
REC_aver <- t(REC_aver)
leaf_aver <- t(leaf_aver)

# Calculate correlation matrix between leaf_aver and REC_aver
celltype.corrmat <- cor(t(leaf_aver), t(REC_aver), method = 'pearson')
celltype.corrmat <- round(celltype.corrmat, 2)

# Scale the correlation matrix
data <- scale(celltype.corrmat)

# Perform hierarchical clustering
ord <- hclust(dist(data, method = "euclidean"), method = "ward.D")$order

# Melt the correlation matrix for plotting
melted_celltype.corrmat <- melt(celltype.corrmat, na.rm = TRUE)

# Set correlation values less than or equal to 0 to NA
melted_celltype.corrmat$value[melted_celltype.corrmat$value <= 0] <- NA

# Order the rows based on hierarchical clustering
melted_celltype.corrmat$Var1 <- factor(melted_celltype.corrmat$Var1, levels = names(table(melted_celltype.corrmat$Var1))[ord])

# Create a heatmap using ggplot2
celltype.corrmat.heatmap <- ggplot(melted_celltype.corrmat, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  theme_bw() +
  coord_fixed() +
  scale_fill_stepsn(n.breaks = 100, colours = colorRampPalette(rev(brewer.pal(9,'RdBu')))(100), na.value = 'lightgrey', limits = c(-0.5, 0.5)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 15)) +
  labs(fill = "Pearson's\ncorrelation") + 
  theme(legend.position = "right", axis.text.x = element_text(angle = 90)) +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5)

# Print the heatmap
celltype.corrmat.heatmap