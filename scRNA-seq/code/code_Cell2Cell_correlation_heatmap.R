# Load required libraries
library(zoo)
library(reshape2)
library(viridis)
library(ggplot2)

# Get top 100 markers for each stage
S1.top100.marker <- head(REC_markerByState$gene[REC_markerByState$cluster %in% 'State 1'], 100)
S2.top100.marker <- head(REC_markerByState$gene[REC_markerByState$cluster %in% 'State 2'], 100)
S3.top100.marker <- head(REC_markerByState$gene[REC_markerByState$cluster %in% 'State 3'], 100)

# Get expression data for selected markers
REC.expr <- REC@assays$originalexp@data[c(S1.top100.marker, S2.top100.marker, S3.top100.marker),]

# Order expression data by pseudotime
REC.expr <- REC.expr[,order(REC$palantir_pseudotime)]

# Apply rolling window to calculate average expression
f_avg_rolling_win <- rollapply(
  data = as.matrix(t(REC.initial.expr)),
  width = 10,  # window width
  FUN = function(w) mean(w),
  align = "right", by = 10,
  partial = TRUE
)

# Calculate correlation matrix
cormat <- round(cor(t(f_avg_rolling_win), method = 'spearman'), 3)
dim(cormat)

# Melt the correlation matrix
melted_cormat <- melt(cormat, na.rm = TRUE)

# Set values below 0.5 to NA
melted_cormat$value[melted_cormat$value < 0.5] <- NA

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  theme_minimal() +  # minimal theme
  coord_fixed() +
  scale_fill_stepsn(n.breaks = 50, colours = viridis::inferno(9), na.value = 'white') +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 15)) +
  labs(fill = "Pearson's\ncorrelation") + theme(legend.position = "right")

# Print the heatmap
print(ggheatmap)