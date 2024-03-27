# Load required packages
library(scater)
library(ggplot2)
library(viridis)

preprocess <- function(sce){
  # Ensure unique feature names
  # rownames(sce) = uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
  
  # Identify mitochondrial and chloroplast genes
  mtgenes = rowData(sce)[grep("ATMG", rowData(sce)$ID),]$ID
  cpgenes = rowData(sce)[grep("ATCG", rowData(sce)$ID),]$ID
  
  # Add per-cell QC metrics
  sce <- addPerCellQC(
    sce,
    subsets = list(MT=mtgenes, CP=cpgenes, MTCP=c(mtgenes, cpgenes)),
    percent_top = c(50, 100, 200, 500),
    detection_limit = 5
  )
  
  # Calculate log10 sum and detected counts
  sce$log10_sum = log10(sce$sum + 1)
  sce$log10_detected = log10(sce$detected + 1)
  
  # Remove cells with zero counts
  sce <- sce[,sce$sum!=0]
  
  return(sce)
}

filtering <- function(sce, umi, mtcppct){
  # Filter cells based on total counts and MTCP percentage
  filter_by_total_counts = sce$sum > umi
  filter_by_mtcp_percent = sce$subsets_MTCP_percent < mtcppct
  
  # Run PCA on QC metrics
  sce <- runColDataPCA(sce, variables = list("sum", "detected", "subsets_MTCP_percent", "percent.top_500"))
  
  # Combine filters
  sce$use <- (
    filter_by_total_counts &
      filter_by_mtcp_percent
  )
  
  # Plot PCA colored by cell use
  plotReducedDim(sce, dimred="PCA_coldata", colour_by="use")
  
  # Subset to filtered cells
  sce = sce[,sce$use]
  
  return(sce)
}

# Preprocess the data
pEXL1_GFP_QC1.sce <- preprocess(pEXL1_GFP_QC1)

# Plot QC histograms
hist(pEXL1_GFP_QC1.sce$log10_sum, breaks = 100)
hist(pEXL1_GFP_QC1.sce$log10_detected, breaks = 100, ylim=c(0,500))
hist(pEXL1_GFP_QC1.sce$subsets_MT_percent, breaks=100)
hist(pEXL1_GFP_QC1.sce$subsets_CP_percent, breaks=100)
hist(pEXL1_GFP_QC1.sce$subsets_MTCP_percent, breaks=100)


# Set filtering parameters
umi = 1000
mtcppct = 10

# Filter the data
pEXL1_GFP_QC2 <- filtering(pEXL1_GFP_QC1.sce, umi, mtcppct)

