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

files <- c("quants.sf") # Path to quant file. (from salmon)

tx2gene <- read.csv("tx2gene.csv") #transcript-to-gene mapping

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

counts <- txi$counts # Extract counts

# Metadata
meta <- data.frame(
  SampleID = colnames(counts),
  Condition = c("Control", "Starvation")
)

rownames(meta) <- meta$SampleID

# NOISeq 

noi_obj <- readData(data = counts, conditions = meta$Condition)


# DESeq2 

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Condition)
dds
dds <- DESeq(dds)
res <- results(dds)
head(res)
# Plot
plotMA(res, main = "MA Plot")
rld <- rlog(dds, blind = TRUE)
plotPCA(rld, intgroup = "Condition")
DE_genes <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
write.csv(DE_genes, "DE_genes_DESeq2.csv")


# Shrink Bayes

shrink_result <- shrinkbayes(counts = counts, design = meta$Condition, method = "bayes")
shrink_result$results
# Plot
ma_plot_data <- data.frame(
  logFC = shrink_result$results$logFC,
  pvalue = shrink_result$results$pvalue
)
ggplot(ma_plot_data, aes(x = logFC, y = -log10(pvalue))) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-log10(p-value)", title = "MA-Plot")




