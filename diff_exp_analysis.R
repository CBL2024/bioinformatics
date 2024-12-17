setwd("C:/Users/Thinkpad/Documents/bioinformatics") # CHANGE THIS

if (!require("BiocManager", quietly = TRUE, )) {
  install.packages("BiocManager") }
if (!require("devtools", quietly = TRUE, )) {
  install.packages("devtools") }
if (!require("tximport", quietly = TRUE, )) {
  BiocManager::install("tximport") }
if (!requireNamespace("NOISeqBIO ", quietly = TRUE)) {
  BiocManager::install("NOISeq") }
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2") }
if (!requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install("edgeR") }
if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma") }
if (!requireNamespace("DEsingle", quietly = TRUE)) {
  BiocManager::install("DEsingle") }
if (!requireNamespace("sleuth", quietly = TRUE)) {
  BiocManager::install("pachterlab/sleuth") }
if (!requireNamespace("baySeq", quietly = TRUE)) {
  BiocManager::install("baySeq") }
if (!requireNamespace("samr", quietly = TRUE)) {
  BiocManager::install("samr") }
if (!requireNamespace("parallel", quietly = TRUE)) {
  BiocManager::install("parallel") }
library(parallel)
library(tximport)

# Import aligned data
files <- list.files(path = "quant_files", pattern = ".sf", full.names = TRUE, 
                    recursive = TRUE)
tx2gene <- read.csv("tx2gene.csv") #transcript-to-gene mapping
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
counts <- txi$counts # Extract counts
colnames(counts) <- c("Sample1", "Sample2", "Sample3")
counts
gene_lengths <- rowMeans(txi$length, na.rm = TRUE) # Keep average gene lengths

# Metadata
meta <- data.frame(
  sample = colnames(counts),
  Condition = c("Control", "Starvation","Starvation"), # Add Starvation also as condition
  path = files
)

rownames(meta) <- meta$sample
head(meta)

# NOISeq 
library(NOISeq)

noi_obj <- readData(data = counts, factor = meta)
noi_res <- noiseq(noi_obj, 
                  factor = "Condition", 
                  norm = "tmm", 
                  replicates = "biological")
# *****This uses tmm normalisation*****
up <- degenes(noi_res, q = 0.9, M = "up")
down <- degenes(noi_res, q = 0.9, M = "down")
noi_degs <- rbind(up, down)
noi_degs


# OUTPUT: noi_degs

# DESeq2 
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Condition)
dds
dds <- DESeq(dds)
res <- results(dds)
deseq_degs <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
deseq_degs

# OUTPUT: DESeq2_degs


library(edgeR)

dge <- DGEList(counts = counts, group = meta$Condition)
dge
dge <- calcNormFactors(dge) # TMM - Trimmed Mean of M-values
dge <- estimateDisp(dge)
dge_result <- exactTest(dge) # the exact test (for 2-group comparison)
edger <- topTags(dge_result)
edger_degs <- edger[edger$PValue < 0.05, ] # Only DEs < 0.05 
edger_degs

# Voom
library(limma)

design <- model.matrix(~ 0 + meta$Condition)
colnames(design) <- levels(meta$Condition)
design # View design matrix
v <- voom(counts, design, plot = FALSE) # Voom transformation
fit <- lmFit(v, design) # Fit with linear model
fit2 <- contrasts.fit(fit, contrasts = c(-1, 1))  # DE analysis between conditions
fit2 <- eBayes(fit2) 
limma <- topTable(fit2) # View table of results
limma_degs <- limma[limma$adj.P.Val < 0.05, ] # Only DEs < 0.05 
limma_degs


#Sleuth
library(sleuth)
so <- sleuth_prep(meta, ~ condition)
so <- sleuth_fit(so, ~ condition)
sleuth_fit_summary(so)
sleuth_degs <- sleuth_wald(so, hypothesis = "conditionTreatment")
sleuth_degs

# bayseq
library(baySeq)
cl <- makeCluster(4)

CD <- new("countData", 
          data = counts, 
          replicates = as.factor(meta$Condition),  
          groups = list(NDE = rep(1, length(meta$Condition)), # NDE = No differential expression = Null hyp
                        DE = as.numeric(meta$Condition))) # DE = Differential expression = Alt hyp


CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl=cl)

CD <- getLikelihoods(CD, bootStraps = 3, verbose = FALSE)

CD@estProps # Estimate psoteria
CD@posteriors[1:10,]
CD@posteriors[101:110,]
CD@estProps[2]

topCounts(CD, group = "DE")
seglens <- mobAnnotation$end - mobAnnotation$start + 1

cD <- new("countData", data = mobData, seglens = seglens, annotation = mobAnnotation)


# SAMseq
library(samr)
samdata <- list(x = as.matrix(counts), y = as.factor(meta$Condition), 
                geneid = rownames(counts), logged2 = FALSE)
sam_res <- SAMseq(samdata$x, samdata$y, resp.type = "Quantitative", nperms = 1000)
sam_degs <- sam_res$siggenes.table
sam_degs
colnames(sam_degs)[1] <- "Gene"

# combined table for all the DE analysis pipelines.
# outputs to combined table: noi_degs, deseq_degs, edger_degs, limma_degs, sleuth_degs, sam_degs
dim(noi_degs)
dim(edger_degs)

combined_results <- merge(deseq_degs, noi_degs, edger_degs, by = col(deseq_degs)[0])
combined_results

gene_ids <- rownames(combined_results)
combined_results$Gene_Length <- gene_lengths[match(gene_ids, names(gene_lengths))]
print(combined_results)
