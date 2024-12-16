### Bioinformatics
# Group 7
# Topic: Identification of the best RNAseq differential expression pipeline to tackle transcript length bias

# Background
# We are looking at how genes are activated or inactivated depending on the composition of the diet. Specifically, we determine the transcriptome of the brain of fruit flies that have been fed on different diets (normal food, fully starved, food with no sugar, food with no amino acids, etc) using RNA-seq. This type of analysis is called Differential Expression analysis, and RNA-seq is a powerful method to do this as it can identify all genes in the genome that are activated or inactivated in response to some external stimulus in a single experimen

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

# Download all data
mkdir ../script
vi ~/Desktop/bioinformatics/fly_assembly/script/load_data.sh
# -------------------------------------------------#
#!/bin/bash

#A simple script to extract multiple sra files

for fn in SRR6315674 SRR6315676 SRR6315680 SRR6315682 SRR6315684 SRR6315686

do

echo now extracting $fn

fasterq-dump $fn -v -e 8

done
# -------------------------------------------------#

# Save the script and exit vi by ESC, then type :wq + ENTER

# Allow editing, reading and executing the script
cd ~/Desktop/bioinformatics/fly_assembly/script
chmod u+rwx load_data.sh

# Run script
cd ~/Desktop/bioinformatics/fly_assembly
bash ~/Desktop/bioinformatics/fly_assembly/script/load_data.sh

# Check the files
ls -lh # -l: long, -h: human-readable - checking file size

# Decompress the files and check the head command lines
zcat ~/Desktop/bioinformatics/fly_assembly/data/*.fastq.gz | head -n 20

# Step 2: QC, trimming and filtering data

vi ~/Desktop/bioinformatics/fly_assembly/script/trim_data.sh

# Paste the following in trim_data.sh

# ---------------------------------------- #
#!/bin/bash

# Specify the directory containing your data
DATA_DIR=./data

# For Complete Starvation
# Fly1: SRR6315688 (f), SRR6315689 (b)
# Fly2: SRR6315686 (f), SRR6315687 (b)
# Fly3: SRR6315684 (f), SRR6315685 (b)
# Fly4: SRR6315682 (f), SRR6315683 (b)

# Define pairs of forward and reverse read IDs
forward_reads=("SRR6315688" "SRR6315686" "SRR6315684" "SRR6315682")
reverse_reads=("SRR6315689" "SRR6315687" "SRR6315685" "SRR6315683")

# Loop through the pairs
for i in ${!forward_reads[@]}; do
    forward=${forward_reads[$i]}
    reverse=${reverse_reads[$i]}
   
    echo "Now trimming $forward and $reverse"

    fastp \
        -i ${DATA_DIR}/${forward}.fastq.gz -I ${DATA_DIR}/${reverse}.fastq.gz \
        -o ${DATA_DIR}/f${forward}.fastq.gz -O ${DATA_DIR}/f${reverse}.fastq.gz \
        --qualified_quality_phred 15 \
        --unqualified_percent_limit 40 \
        --length_required 80 \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --thread 8 \
        --html ${DATA_DIR}/${forward}.html \
        --json ${DATA_DIR}/${forward}.json
done

# For Normal Feeding
# Fly1: SRR6315680 (f), SRR6315681 (b)
# Fly2: SRR6315678 (f), SRR6315679 (b)
# Fly3: SRR6315676 (f), SRR6315677 (b)
# Fly4: SRR6315674 (f), SRR6315675 (b)

# Define pairs of forward and reverse read IDs
forward_reads=("SRR6315680" "SRR6315678" "SRR6315676" "SRR6315674")
reverse_reads=("SRR6315681" "SRR6315679" "SRR6315677" "SRR6315675")

# Loop through the pairs
for i in ${!forward_reads[@]}; do
    forward=${forward_reads[$i]}
    reverse=${reverse_reads[$i]}
   
    echo "Now trimming $forward and $reverse"

    fastp \
        -i ${DATA_DIR}/${forward}.fastq.gz -I ${DATA_DIR}/${reverse}.fastq.gz \
        -o ${DATA_DIR}/f${forward}.fastq.gz -O ${DATA_DIR}/f${reverse}.fastq.gz \
        --qualified_quality_phred 15 \
        --unqualified_percent_limit 40 \
        --length_required 80 \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --thread 8 \
        --html ${DATA_DIR}/${forward}.html \
        --json ${DATA_DIR}/${forward}.json
done
# ---------------------------------------- #

# Save the script and exit vi by ESC, then type :wq + ENTER

# Allow editing, reading and executing the script
cd ~/Desktop/bioinformatics/fly_assembly/script
chmod u+rwx trim_data.sh

# Run script
cd ~/Desktop/bioinformatics/fly_assembly
bash ~/Desktop/bioinformatics/fly_assembly/script/trim_data.sh

# Step 4: Count tables using salmon
# Note: we want to find out how these collections of short reads related to gene expression
# Here, we will use salmon (pseudo-aligner: which is faster and less resource-intensive than aligners), given that Drophila has a very well-annotated genome

# Download Drosophila transcriptome
wget https://ftp.ensembl.org/pub/release-106/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz -P ./data/

# Generate an index file with the following command
salmon index -t ./data/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz -i ./data/Dmel6.32

# Write a script for running salmon in all trimmed data
vi ~/Desktop/bioinformatics/fly_assembly/script/salmon.sh

# ---------------------------------------- #
#!/bin/bash

forward_reads=("SRR6315680" "SRR6315678" "SRR6315676" "SRR6315674" "SRR6315688" "SRR6315686" "SRR6315684" "SRR6315682")
reverse_reads=("SRR6315681" "SRR6315679" "SRR6315677" "SRR6315675" "SRR6315689" "SRR6315687" "SRR6315685" "SRR6315683")

# Ensure the quants directory exists
mkdir -p quants

for i in "${!forward_reads[@]}"
do
    fn_forward=${forward_reads[$i]}
    fn_reverse=${reverse_reads[$i]}
   
    echo "Now aligning $fn_forward with $fn_reverse"

    salmon quant -i ./data/Dmel6.32 -l A -1 ./data/f${fn_forward}.fastq.gz -2 ./data/f${fn_reverse}.fastq.gz -p 8 -o ./data/quants/${fn_forward} > ./data/quants/${fn_forward}_salmon.log 2>&1
done
# ---------------------------------------- #

# Save the script and exit vi by ESC, then type :wq + ENTER

# Allow editing, reading and executing the script
cd ~/Desktop/bioinformatics/fly_assembly/script
chmod u+rwx salmon.sh

# Run script
cd ~/Desktop/bioinformatics/fly_assembly
bash ~/Desktop/bioinformatics/fly_assembly/script/salmon.sh
