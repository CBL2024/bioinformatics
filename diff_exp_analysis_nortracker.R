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
if (!requireNamespace("parallel", quietly = TRUE)) {
  BiocManager::install("parallel") }
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
counts <- round(counts)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Condition)
dds
dds <- DESeq(dds)
deseq_res <- results(dds)
deseq_res


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
colnames(design) <- c("Control", "Starvation")
design # View design matrix
v <- voom(counts, design, plot = FALSE) # Voom transformation
fit <- lmFit(v, design) # Fit with linear model

contrast.matrix <- makeContrasts(
  ConditionB_vs_ConditionA = Starvation - Control, 
  levels = design
)
contrast.matrix

fit <- contrasts.fit(fit, contrasts = contrast.matrix)  # DE analysis between conditions
fit <- eBayes(fit) 
limma_res <- topTable(fit, adjust = "BH", sort.by = "P", number = Inf)  


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
groups

CD <- new("countData", 
          data = counts, 
          replicates = as.factor(meta$Condition),  
          groups = groups)

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

write.csv(topCounts(CD2, group = "DE"),"~/bioinformatics/posteriors.csv")

stopCluster(cl)

topCounts(CD2, group = "DE")

plotMA.CD(CD2, samplesA = "control", samplesB = "starvation",
          col = c(rep("red", 100), rep("black", 900)))
