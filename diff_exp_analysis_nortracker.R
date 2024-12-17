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
  BiocManager::install("pachterlab/sleuth") }"quant_files"
if (!requireNamespace("baySeq", quietly = TRUE)) {
  BiocManager::install("baySeq") }
if (!requireNamespace("samr", quietly = TRUE)) {
  BiocManager::install("samr") }
if (!requireNamespace("parallel", quietly = TRUE)) {
  BiocManager::install("parallel") }
BiocManager::install("rtracklayer")
library(parallel)
library(tximport)
library(tracklayer)
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

# Import aligned data
# files <- list.files(path = "quant_files", pattern = ".sf", full.names = TRUE, 
#                     recursive = TRUE)
# tx2gene <- read.csv("tx2gene.csv") #transcript-to-gene mapping

# reference file
gtf_file_path <- "~/bioinformatics/Drosophila_melanogastr.BDGP6.46.113.gtf"
gtf_data <- tracklayer::import(gtf_file_path)
tx2gene <- data.frame(
  transcript_id = mcols(dtf_data)$transcript_id,
  gene_id = mcols(dtf_data)$gene_id
)

# load reference data file #####################################################

# Path to your GTF file
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
head(tx2gene)

################################################################################

# Get count data
txi <- tximport(paths_df$paths, type = "salmon", tx2gene = tx2gene)
counts <- txi$counts # Extract counts
gene_lengths <- rowMeans(txi$length, na.rm = TRUE) # Keep average gene lengths
gene_names <- rownames(counts)

# save this data
write.csv(gene_lengths, "~/bioinformatics/gene_lengths.csv", row.names = FALSE)
write.csv(gene_names, "~/bioinformatics/gene_names.csv", row.names = FALSE)

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
print(CD)

CD@estProps # Estimate psoteria
CD@posteriors[1:10,]
CD@posteriors[101:110,]
CD@estProps[2]

topCounts(CD, group = "DE")

results <- getDE(CD)


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
