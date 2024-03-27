## code_QC_part_1

# Load required library
library(DropletUtils)

# Read 10x data
pEXL1_GFP <- read10xCounts('D:/project/pEXL1_GFP/raw_gene_bc_matrices/')

# Calculate barcode ranks
br.out <- barcodeRanks(counts(pEXL1_GFP))

# Save barcode rank plot
png(paste0(' ', rawsce_list[i], '.png'))

# Plot barcode ranks
o <- order(br.out$rank)
plot(br.out$rank[o], br.out$total, log = "xy", xlab = "Rank", ylab = "Total", main = rawsce_list[i])
lines(br.out$rank[o], br.out$fitted[o], col = "red")

# Perform emptyDrops analysis (identify empty droplets)
set.seed(2022)
e.out <- emptyDrops(counts(rawsce)) 

# Summarize emptyDrops results
table(Sig=e.out$FDR <= 0.05, Limited=e.out$Limited)

# Select barcodes identified as cells
is.cell <- e.out$FDR <= 0.05
print(sum(is.cell, na.rm=TRUE))
print(table(br.out$rank == sum(is.cell, na.rm=TRUE)))

# Draw lines on the barcode rank plot
abline(h=min(br.out$fitted[o], na.rm=TRUE), col="red", lty=2)
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), legend=c("knee", "inflection", "FDR_0.05"))

# Save the plot
saved_plot <- recordPlot()
dev.off()

# Update rawsce object by selecting only barcodes identified as cells
colnames(rawsce) = colData(rawsce)$Barcode
pEXL1_GFP_QC1 <- rawsce[,which(e.out$FDR <= 0.05)]