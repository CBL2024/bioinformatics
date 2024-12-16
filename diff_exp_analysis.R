# Load necessary libraries
install("devtools")

if (!require("BiocManager", quietly = TRUE, ))
  install.packages("BiocManager")
BiocManager::install("tximport")

if (!requireNamespace("NOISeq", quietly = TRUE)) {
  BiocManager::install("NOISeq")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
if (!requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install("edgeR")
}
if (!requireNamespace("voom", quietly = TRUE)) {
  BiocManager::install("voom")
}
if (!requireNamespace("DESingle", quietly = TRUE)) {
  BiocManager::install("DESingle")
}

install.packages(c("sp","pixmap", "snowfall", "VGAM", "mclust", "logcondens", "Iso","XML","rgl"), repos="http://cran.r-project.org")
if (!requireNamespace("ShrinkBayes", quietly = TRUE)) {
  BiocManager::install("ShrinkBayes")
}

library("devtools")
install.packages("INLA", repos = c("https://inla.r-inla-download.org/R/stable", "https://cran.r-project.org"))
install_github("markvdwiel/ShrinkBayes")

setwd("C:/Users/ThinkPad/Documents/bioinformatics") # CHANGE THIS
library("tximport")
library("NOISeq")
library("DESeq2")
library("ShrinkBayes")
library("edgeR")
library("voom")
library("DESingle")

files <- c("quants.sf","quants2.sf","quants3.sf","quants4.sf") # Path to quant file. (from salmon)
print(files) #PRINT
tx2gene <- read.csv("tx2gene.csv") #transcript-to-gene mapping

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
counts <- txi$counts # Extract counts

head(counts) #PRINT
colnames(counts) <- c("Sample1", "Sample2", "Sample3", "Sample4")

# Metadata
meta <- data.frame(
  SampleID = colnames(counts),
  Condition = c("Control","Control","Starvation","Starvation") # Add Starvation also as condition
)

rownames(meta) <- meta$SampleID

head(meta)

# NOISeq 

noi_obj <- readData(data = counts, conditions = meta$Condition)

# DESeq2 

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Condition)
dds
dds <- DESeq(dds)
DESeq2_table <- results(dds)
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


# DGEList
dge <- DGEList(counts = counts, group = meta$Condition)
dge
dge <- calcNormFactors(dge) # Normalise
dge$samples # View normalised data
dge_result <- exactTest(dge)
topTags(dge_result)

#Plot
plotMA(dge_result, main = "MA-Plot")
rld <- rlog(dge)  # rlog transformation for PCA 
plotPCA(rld, intgroup = "Condition") # Plot PCA



# Voom

design <- model.matrix(~ 0 + meta$Condition)
colnames(design) <- levels(meta$Condition)
design # View design matrix

v <- voom(counts, design, plot = TRUE) # Voom transformation
head(v$E)  # E = transformed expression data, view the voom transformation

fit <- lmFit(v, design) # Fit with linear model
fit2 <- contrasts.fit(fit, contrasts = c(-1, 1))  # DE analysis between conditions
fit2 <- eBayes(fit2) 
topTable(fit2) # View table of results

# Plot
plotMA(fit2, main = "MA-Plot", ylim = c(-5, 5))
voom_table <- topTable(fit2, adjust = "fdr", p.value = 0.05)


# Create a combined table for all the DE analysis pipelines.

combined_results <- Reduce(function(x, y) merge(x, y, by = "GeneID", all = TRUE),
                           list(DESeq2_table, edgeR_table, voom_table))
