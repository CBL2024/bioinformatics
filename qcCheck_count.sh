### Bioinformatics
# Group 7
# Topic: Identification of the best RNAseq differential expression pipeline to tackle transcript length bias

# Background
# We are looking at how genes are activated or inactivated depending on the composition of the diet. Specifically, 
# we determine the transcriptome of the brain of fruit flies that have been fed on different diets (normal food, 
# fully starved, food with no sugar, food with no amino acids, etc) using RNA-seq. This type of analysis is called 
# Differential Expression analysis, and RNA-seq is a powerful method to do this as it can identify all genes in the 
# genome that are activated or inactivated in response to some external stimulus in a single experiment.

# Conditions
# 2 groups: Fully starved vs. control (not starved)
# Analysis methods (SESeq, edgeR, baySeq, SAMseq)

# Data
# Raw data: all raw data are in .fastq.gz files, and all the .fastq.gz files for each condition per fly has two files (e.g., 1 forward (f)and 1 backward (b) sequences)

# Raw data list
# Complete Starvation
# Fly1: SRR6315688 (f), SRR6315689 (b)
# Fly2: SRR6315686 (f), SRR6315687 (b)
# Fly3: SRR6315684 (f), SRR6315685 (b)
# Fly4: SRR6315682 (f), SRR6315683 (b)

# Normal feeding
# Fly1: SRR6315680 (f), SRR6315681 (b)
# Fly2: SRR6315678 (f), SRR6315679 (b)
# Fly3: SRR6315676 (f), SRR6315677 (b)
# Fly4: SRR6315674 (f), SRR6315675 (b)

# Note:
# This is a script for aligning and assembling the sequence of neuropiptide sequences in fruit flies (e.g., control vs. fully starved conditions)
# Unless specified, all of the following scripts are in BASH

# Step 1: Manage data
# Making new directories
mkdir ~/Desktop/bioinformatics/fly_assembly
mkdir ~/Desktop/bioinformatics/fly_assembly/data
cd ~/Desktop/bioinformatics/fly_assembly/data # Check if you are now in the new folder

# Put all your data from Download to your folder
cp ~/Download/*fastq.gz ~/Desktop/fly_assembly/data/ # Check again if all your files are present and without any duplicate files

# Check the files
ls -lh # -l: long, -h: human-readable - checking file size

# Decompress the files and check the head command lines
zcat ~/Desktop/bioinformatics/fly_assembly/data/*.fastq.gz | head -n 20

# Step 2: Running QC checks, trimming and filtering data
cd ../ # cd pack to the fly_assembly folder from the fly_assembly/data folder
mkdir script
vi ~Desktop/bioinformatics/fly_assembly/script/trim_data.sh

# Paste the following in trim_data.sh

# ---------------------------------------- #
#!/bin/bash

#A simple script to trim multiple read files

# Specify the directory containing your data
DATA_DIR=./data

for fn in SRR6315688 SRR6315686 SRR6315684 SRR6315682
do
    echo "Now trimming $fn"
    fastp \
        -i ${DATA_DIR}/${fn}_1.fastq.gz -I ${DATA_DIR}/${fn}_2.fastq.gz \
        -o ${DATA_DIR}/f${fn}_1.fastq.gz -O ${DATA_DIR}/f${fn}_2.fastq.gz \
        --qualified_quality_phred 15 \
        --unqualified_percent_limit 40 \
        --length_required 80 \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --thread 8 \
        --html ${DATA_DIR}/${fn}.html \
        --json ${DATA_DIR}/${fn}.json
done
# ---------------------------------------- #

# Save the script by pressing ESC, then type :wq + ENTER

# Allow editing, reading and executing the script
cd ~/Desktop/bioinformatics/fly_assembly/script
chmod u+rwx trim_data.sh

# Run script
bash ~/Desktop/bioinformatics/fly_assembly/script/trim_data.sh

# Step 4: Count tables using salmon
# Note: we want to find out how these collections of short reads related to gene expression
# Here, we will use salmon (pseudo-aligner: which is faster and less resource-intensive than aligners), given that Drophila has a very well-annotated genome

# Download Drosophila transcriptome
# Go to the ensembl ftp site: https://www.ensembl.org/info/data/ftp/index.html
# Search for Drosophila, but this time click on FASTA under cDNA.
# Get the download link for Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz
# Download this file using wget

# Make new directory for output
cd ~/Desktop/bioinformatics/fly_assembly
mkdir ./data/quants

# Generate an index file with the following command
vi ~/Desktop/bioinformatics/fly_assembly/script/salmon.sh

# ---------------------------------------- #
#!/bin/bash
# A simple script to pseudo-align reads

# Specify the directories for input and output
DATA_DIR=./data
OUTPUT_DIR=./data/quants

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

for fn in SRR6315688 SRR6315686 SRR6315684 SRR6315682
do
    echo "Now aligning $fn"
    salmon quant \
        -i Dmel6.32 \                         # Reference transcriptome
        -l A \                                # Library type
        -1 ${DATA_DIR}/f${fn}_1.fastq.gz \    # Input file 1 from data folder
        -2 ${DATA_DIR}/f${fn}_2.fastq.gz \    # Input file 2 from data folder
        -p 8 \                                # Number of threads
        -o ${OUTPUT_DIR}/${fn}                # Output directory for each sample
done
# ---------------------------------------- #

# Save the script by pressing ESC, then type :wq + ENTER

# Allow editing, reading and executing the script
cd ~/Desktop/bioinformatics/fly_assembly/script
chmod u+rwx salmon.sh

# Run script
cd ~/Desktop/bioinformatics/fly_assembly
bash ~/Desktop/bioinformatics/fly_assembly/script/salmon.sh
