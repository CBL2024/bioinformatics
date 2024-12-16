# Load necessary libraries
if (!requireNamespace("NOISeq", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("NOISeq")
}
library(NOISeq)
library(tximport)

files <- c("sample1/quant.sf", "sample2/quant.sf")
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

