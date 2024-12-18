if (!require("BiocManager", quietly = TRUE, )) {
  install.packages("BiocManager") }
if (!require("devtools", quietly = TRUE, )) {
  install.packages("devtools") }
if (!require("tximport", quietly = TRUE, )) {
  BiocManager::install("tximport") }
if (!requireNamespace("NOISeqBIO ", quietly = TRUE)) {
  BiocManager::install("noiseqbio") }
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
if (!requireNamespace("parallel", quietly = TRUE)) {
  BiocManager::install("parallel") }
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  BiocManager::install("ggplot2") }
library(parallel)
library(tximport)
library(DESeq2)

# Define control and starvation file paths
base_path <- "~/bioinformatics/quants"
prefix <- "SRR63156"
suffix_control <- c("74", "76", "78", "80")
suffix_starvation <- c("82", "84", "86", "88")

# Define vectors of file paths and conditions and combine into DF
paths <- c()
conditions <- c()
for (suffix in suffix_control) {
  file_path <- paste(base_path, "/", prefix, suffix, "_1/", "quant.sf", sep="")
  paths <- c(paths, file_path)
  conditions <- c(conditions, "control")
}
for (suffix in suffix_starvation) {
  file_path <- paste(base_path, "/", prefix, suffix, "_1/", "quant.sf", sep="")
  paths <- c(paths, file_path)
  conditions <- c(conditions, "starvation")
}
paths_df <- data.frame(conditions = conditions, paths = paths)

# reference file
gtf_file_path <- "~/bioinformatics/Drosophila_melanogaster.BDGP6.46.113.gtf"

# Read the GTF file using read.table, skipping comments (lines starting with "#")
gtf_data <- read.table(gtf_file_path, sep = "\t", header = FALSE, comment.char = "#", stringsAsFactors = FALSE)

# Assign column names based on the GTF format
colnames(gtf_data) <- c("seqname", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attribute")

# Extract the 'transcript_id' and 'gene_id' from the 'attribute' column
extract_attributes <- function(attribute_column, key) {
  sapply(strsplit(attribute_column, "; "), function(x) {
    match <- grep(paste0("^", key), x, value = TRUE)
    if (length(match) > 0) {
      gsub(paste0(key, " "), "", match)
    } else {
      NA
    }
  })
}

# Extract transcript_id and gene_id
gtf_data$transcript_id <- extract_attributes(gtf_data$attribute, "transcript_id")
gtf_data$gene_id <- extract_attributes(gtf_data$attribute, "gene_id")

# Create the tx2gene data frame
tx2gene <- data.frame(
  transcript_id = gtf_data$transcript_id,
  gene_id = gtf_data$gene_id
)

# Remove rows with NA values (optional, but can be useful)
tx2gene <- na.omit(tx2gene)

# View the first few rows of the tx2gene data frame

# Get count data
txi <- tximport(paths_df$paths, type = "salmon", tx2gene = tx2gene)
counts <- txi$counts # Extract counts
txi$length
gene_lengths <- data.frame(row.names = rownames(counts), rowMeans(txi$length, na.rm = TRUE)) 
gene_lengths
write.csv(gene_lengths, "~/bioinformatics/gene_lengths.csv")

# Metadata
meta <- data.frame(
  sample = paths_df$paths, # colnames(counts),
  Condition = paths_df$conditions,
  path = paths_df$paths
)

rownames(meta) <- meta$sample
head(meta)

# NOISeq 
library(NOISeq)

noi_obj <- readData(data = counts, factor = meta)
noi_res <- noiseqbio(noi_obj, 
                  factor = "Condition", 
                  norm = "tmm")
# *****This uses tmm normalisation*****
up <- degenes(noi_res, q = 0.9, M = "up")
down <- degenes(noi_res, q = 0.9, M = "down")
noi_degs <- rbind(up, down)
noi_degs
write.csv(noi_res@results,"~/bioinformatics/noi_res.csv")
DE.plot(noi_res, q = 0.9, graphic = "expr", log.scale = TRUE)

# OUTPUT: noi_degs

# DESeq2 
library(DESeq2)
counts <- round(counts)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Condition)
dds <- DESeq(dds)
deseq_res <- results(dds)
topTable(dds)

# OUTPUT: DESeq2_degs


library(edgeR)

dge <- DGEList(counts = counts, group = meta$Condition)
dge <- calcNormFactors(dge) # TMM - Trimmed Mean of M-values
dge <- estimateDisp(dge)
dge_result <- exactTest(dge) # the exact test (for 2-group comparison)
edger <- topTags(dge_result)
edger_degs <- edger[edger$PValue < 0.05, ] # Only DEs < 0.05 

# Voom
library(limma)

design <- model.matrix(~ 0 + meta$Condition)  
colnames(design) <- c("Control", "Starvation")
design # View design matrix
v <- voom(counts, design, plot = T) # Voom transformation
fit <- lmFit(v, design) # Fit with linear model

contrast.matrix <- makeContrasts(
  ConditionB_vs_ConditionA = Starvation - Control, 
  levels = design
)
contrast.matrix

fit <- contrasts.fit(fit, contrasts = contrast.matrix)  # DE analysis between conditions
fit <- eBayes(fit) 
limma_res <- topTable(fit, adjust = "BH", sort.by = "P", number = Inf)  


library(samr)
samdata <- list(x = as.matrix(counts), y = as.factor(meta$Condition), 
                geneid = rownames(counts), genenames = rownames(counts), logged2=FALSE)

samr.obj <- samr(samdata, resp.type="Two class unpaired", nperms=1000)
pv=samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)

write.csv(pv,"~/bioinformatics/samr_res.csv")
write.csv(samr.obj$foldchange,"~/bioinformatics/samr_logfc.csv")



#Sleuth
library(sleuth)
so <- sleuth_prep(meta, ~ condition)
so <- sleuth_fit(so, ~ condition)
sleuth_fit_summary(so)
sleuth_degs <- sleuth_wald(so, hypothesis = "conditionTreatment")
sleuth_degs

# bayseq
library(baySeq)
Condition_factor <- as.factor(meta$Condition)

cl <- makeCluster(16)

groups <- list(
  NDE = rep(1, length(Condition_factor)),  # Null hypothesis: All samples in one group
  DE  = as.numeric(Condition_factor)       # Alternative hypothesis: 1 = Control, 2 = Starvation
)
replicates <- as.character(meta$Condition)

CD <- new("countData", 
          data = counts, 
          replicates = replicates,
          groups = groups,
          densityFunction = bbDensity)

libsizes(CD) <- getLibsizes(CD)
libsizes(CD)

CD@annotation <- data.frame(name = row.names(counts))

CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl=cl)

CD@groups
sapply(names(CD@groups), function(group) lapply(CD@priors$priors[[group]], head, 5))


CD2 <- getLikelihoods(CD, bootStraps = 10, verbose = TRUE)

CD2@estProps
CD2@posteriors[1:10,]
CD2@posteriors[101:110,]
CD2@estProps[2]

write.csv(cbind(CD2@annotation, CD2@posteriors),"~/bioinformatics/bayseq_beta.csv")

stopCluster(cl)


plotMA.CD(CD2, samplesA = "control", samplesB = "starvation",
          col = c(rep("red", 100), rep("black", 900)))

topCounts(CD2, group = "DE", normaliseData = TRUE)  













bayseq <- read.csv("~/bioinformatics/bayseq_res.csv")
merged_df <- merge(bayseq, gene_lengths, by.x = "name", by.y = "row.names", all.x = TRUE)
hist(merged_df$rowMeans.txi.length..na.rm...TRUE., main="BaySeq", xlim=c(0,20000),
     col = rgb(1, 0, 0, 0.5), breaks = 50, xlab = "Gene Length")

bayes_de <- subset(merged_df, DE > -0.05 )
hist(bayes_de$rowMeans.txi.length..na.rm...TRUE., breaks = 50, 
     col = rgb(0, 0, 1, 0.5), add = TRUE )

legend("topright", legend = c("All genes", "Differentially expressed"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))

write.csv(merged_df, "~/bioinformatics/DE_results_with_length/bay_lengths.csv")


limma <- read.csv("~/bioinformatics/limma_res.csv")
merged_df <- merge(limma, gene_lengths, by.x = "name", by.y = "row.names", all.x = TRUE)
hist(merged_df$rowMeans.txi.length..na.rm...TRUE., main="Limma", xlim=c(0,20000),
     col = rgb(1, 0, 0, 0.5), breaks = 50, xlab = "Gene Length")

limma_de <- subset(merged_df, adj.P.Val < 0.05 )
hist(limma_de$rowMeans.txi.length..na.rm...TRUE., breaks = 50, 
     col = rgb(0, 0, 1, 0.5), add = TRUE )

legend("topright", legend = c("All genes", "Differentially expressed"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))

write.csv(merged_df, "~/bioinformatics/DE_results_with_length/limma_lengths.csv")


edger <- read.csv("~/bioinformatics/edger_res.csv")
merged_df <- merge(edger, gene_lengths, by.x = "name", by.y = "row.names", all.x = TRUE)
hist(merged_df$rowMeans.txi.length..na.rm...TRUE., main="EdgeR", xlim=c(0,20000),
     col = rgb(1, 0, 0, 0.5), breaks = 50, xlab = "Gene Length")

edger_de <- subset(merged_df, PValue < 0.05 )
hist(edger_de$rowMeans.txi.length..na.rm...TRUE., breaks = 50, 
     col = rgb(0, 0, 1, 0.5), add = TRUE )

legend("topright", legend = c("All genes", "Differentially expressed"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))

write.csv(merged_df, "~/bioinformatics/DE_results_with_length/edger_lengths.csv")


noiseq <- read.csv("~/bioinformatics/noi_res.csv")
merged_df <- merge(noiseq, gene_lengths, by.x = "name", by.y = "row.names", all.x = TRUE)
hist(merged_df$rowMeans.txi.length..na.rm...TRUE., main="NOISeq", xlim=c(0,20000),
     col = rgb(1, 0, 0, 0.5), breaks = 50, xlab = "Gene Length")

noiseq_de <- subset(merged_df, prob > 0.95 )
hist(noiseq_de$rowMeans.txi.length..na.rm...TRUE., breaks = 50, 
     col = rgb(0, 0, 1, 0.5), add = TRUE )

legend("topright", legend = c("All genes", "Differentially expressed"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))

svg("~/bioinformatics/NOISeq_hist.svg")

write.csv(merged_df, "~/bioinformatics/DE_results_with_length/noiseq_lengths.csv")


deseq <- read.csv("~/bioinformatics/deseq_res.csv")
merged_df <- merge(deseq, gene_lengths, by.x = "name", by.y = "row.names", all.x = TRUE)
hist(merged_df$rowMeans.txi.length..na.rm...TRUE., main="DESeq2", xlim=c(0,20000),
     col = rgb(1, 0, 0, 0.5), breaks = 50, xlab = "Gene Length")

deseq_de <- subset(merged_df, padj < 0.05 )
hist(deseq_de$rowMeans.txi.length..na.rm...TRUE., breaks = 50, 
     col = rgb(0, 0, 1, 0.5), add = TRUE )

legend("topright", legend = c("All genes", "Differentially expressed"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))

write.csv(merged_df, "~/bioinformatics/DE_results_with_length/deseq_lengths.csv")

samr <- read.csv("~/bioinformatics/samr_res.csv")

merged_df <- merge(samr, gene_lengths, by.x = "name", by.y = "row.names", all.x = TRUE)
merged_df
hist(merged_df$rowMeans.txi.length..na.rm...TRUE., main="samr", xlim=c(0,20000),
     col = rgb(1, 0, 0, 0.5), breaks = 50, xlab = "Gene Length")

samr_de <- subset(merged_df, x < 0.05 )
hist(samr_de$rowMeans.txi.length..na.rm...TRUE., breaks = 50, 
     col = rgb(0, 0, 1, 0.5), add = TRUE )

legend("topright", legend = c("All genes", "Differentially expressed"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))

write.csv(merged_df, "~/bioinformatics/DE_results_with_length/samr_lengths.csv")



# Find the common genes across methods


common_genes <- intersect(
  samr_de$name,    # Extract gene identifiers from SAMR result
  intersect(
    deseq_de$name, # Extract gene identifiers from DESeq2 result
    intersect(
      limma_de$name, # Extract gene identifiers from limma result
      bayes_de$name   # Extract gene identifiers from bayes result
    )
  )
)

common_genes

counts_meta <- counts
colnames(counts_meta) <- meta$Condition
counts_meta

# Extract expression data for the common genes
expression_data <- counts_meta[common_genes, ]  # Assuming 'expression_matrix' has gene names as row names

# Create a heatmap using the 'pheatmap' package
library(pheatmap)

# Plot the heatmap
expression <- counts_meta[edger_de$name, ]

hmap <- pheatmap(expression, scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", main="edgeR heatmap",
         show_colnames = TRUE)

expression <- counts_meta[noiseq_de$name, ]

hmap <- pheatmap(expression, scale = "row", clustering_distance_rows = "euclidean", 
                 clustering_distance_cols = "euclidean", main="NOISeq heatmap",
                 show_colnames = TRUE)

expression <- counts_meta[limma_de$name, ]

hmap <- pheatmap(expression, scale = "row", clustering_distance_rows = "euclidean", 
                 clustering_distance_cols = "euclidean", main="limma heatmap",
                 show_colnames = TRUE)

expression <- counts_meta[bayes_de$name, ]

hmap <- pheatmap(expression, scale = "row", clustering_distance_rows = "euclidean", 
                 clustering_distance_cols = "euclidean", main="baySeq heatmap",
                 show_colnames = TRUE)

expression <- counts_meta[deseq_de$name, ]

hmap <- pheatmap(expression, scale = "row", clustering_distance_rows = "euclidean", 
                 clustering_distance_cols = "euclidean", main="DESeq2 heatmap",
                 show_colnames = TRUE)

expression <- counts_meta[edger_de$name, ]

hmap <- pheatmap(expression, scale = "row", clustering_distance_rows = "euclidean", 
                 clustering_distance_cols = "euclidean", main="edgeR heatmap",
                 show_colnames = TRUE)

expression <- counts_meta[samr_de$name, ]

hmap <- pheatmap(expression, scale = "row", clustering_distance_rows = "euclidean", 
                 clustering_distance_cols = "euclidean", main="samR heatmap",
                 show_colnames = TRUE)



library(ggplot2)
deseq$neg_log10_pvalue <- -log10(deseq$padj)
# Create a volcano plot
ggplot(deseq, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
  geom_point(aes(color = pvalue < 0.05), size = 1) +  # Color points based on p-value significance
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)", title = "DESeq2") +
  scale_color_manual(values = c("gray", "black"))  # Red for significant points

limma$neg_log10_pvalue <- -log10(limma$adj.P.Val)
# Create a volcano plot
ggplot(limma, aes(x = logFC, y = neg_log10_pvalue)) +
  geom_point(aes(color = adj.P.Val < 0.05), size = 1) +  # Color points based on p-value significance
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)", title = "limma") +
  scale_color_manual(values = c("gray", "black"))  # Red for significant points

edger$neg_log10_pvalue <- -log10(edger$PValue)
# Create a volcano plot
ggplot(edger, aes(x = logFC, y = neg_log10_pvalue)) +
  geom_point(aes(color = PValue < 0.05), size = 1) +  # Color points based on p-value significance
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)", title = "edgeR") +
  scale_color_manual(values = c("gray", "black"))  # Red for significant points


noiseq <- na.omit(noiseq)
noiseq$neg_log10_pvalue <- -log10(1-noiseq$prob)
# Create a volcano plot
ggplot(noiseq, aes(x = log2FC, y = neg_log10_pvalue)) +
  geom_point(aes(color = 1-prob < 0.05), size = 1) +  # Color points based on p-value significance
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)", title = "noiseq") +
  scale_color_manual(values = c("gray", "black"))  # Red for significant points


samr <- read.csv("~/bioinformatics/samr_res.csv")

samr$neg_log10_pvalue <- -log10(samr$pvalue)
# Create a volcano plot
ggplot(samr, aes(x = logFC, y = neg_log10_pvalue)) +
  geom_point(aes(color = pvalue < 0.05), size = 1) +  # Color points based on p-value significance
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)", title = "samr") +
  scale_color_manual(values = c("gray", "black"))  # Red for significant points


