# Load necessary libraries
if (!requireNamespace("NOISeq", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("NOISeq")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
if (!requireNamespace("ShrinkBayes", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("ShrinkBayes")
}
if (!requireNamespace("ShrinkBayes", quietly = TRUE)) {
      install.packages("BiocManager")
      BiocManager::install("ShrinkBayes")
}
if (!requireNamespace("edgeR", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("edgeR")
}
if (!requireNamespace("voom", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("voom")
}
if (!requireNamespace("DESingle", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("DESingle")
}

library(tximport)

# Packages for DE analysis:
library(NOISeq)
library(DESeq2)
library(ShrinkBayes)
library(edgeR)
library(voom)
library(DESingle)

  


files <- c("sample1.sf", "sample2.sf, "sample3.sf", "sample4.sf", "sample5.sf")
tx2gene <- read.csv("tx2gene.csv") # Contains transcript-to-gene mapping

# Import counts
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Access gene-level counts
counts <- txi$counts

# Metadata
meta <- data.frame(
  SampleID = colnames(counts),
  Condition = c("Control", "Control", "Starvation", "Starvation")
)

rownames(meta) <- meta$SampleID

# Create the NOISeq object
# `counts` is the count matrix, and `meta` is the sample metadata
data <- readData(data = counts, conditions = meta$Condition)

