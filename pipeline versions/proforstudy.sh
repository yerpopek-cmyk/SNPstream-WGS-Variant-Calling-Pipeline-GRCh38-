#!/bin/bash
# Stop the script on any critical error
set -e  

# --- PIPELINE PARAMETERS ---
ACCESSION="SRR622461"
REFERENCE="GCA_000001405.29.fasta"
SAMPLE="sample"
THREADS=6

# --- CONDA ENVIRONMENTS ---
ENV_PREP="ngs_prep"       # Env for sra-tools, fastqc, fastp, bwa
ENV_BAM="ngs_bam"         # Env for samtools and bcftools
# -----------------

# [FUNCTION] Check if required tools are installed
check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo "[ERROR] Program $1 is not installed. Aborting."
        exit 1
    fi
}

echo "========================================="
echo "[INFO] Activating PREP environment: $ENV_PREP"
eval "$(conda shell.bash hook)"
conda activate $ENV_PREP
echo "========================================="

echo "[STEP 1/17] Downloading raw reads..."
prefetch $ACCESSION
fasterq-dump $ACCESSION --split-files --threads $THREADS
mv ${ACCESSION}_1.fastq ${SAMPLE}_R1.fastq
mv ${ACCESSION}_2.fastq ${SAMPLE}_R2.fastq

echo "[STEP 2-4/17] Quality Control and Trimming..."
mkdir -p fastqc_before fastp_reports fastqc_after
fastqc ${SAMPLE}_R1.fastq ${SAMPLE}_R2.fastq -t $THREADS -o fastqc_before/
fastp -i ${SAMPLE}_R1.fastq -I ${SAMPLE}_R2.fastq \
      -o ${SAMPLE}_R1_trimmed.fastq -O ${SAMPLE}_R2_trimmed.fastq \
      -h fastp_reports/fastp_report.html -j fastp_reports/fastp_report.json --thread $THREADS
fastqc ${SAMPLE}_R1_trimmed.fastq ${SAMPLE}_R2_trimmed.fastq -t $THREADS -o fastqc_after/

echo "[STEP 5/17] Indexing Reference..."
if [ ! -f "${REFERENCE}.bwt" ]; then bwa index $REFERENCE; fi

echo "[STEP 6/17] Alignment: Creating RAW SAM file..."
# Read Group for compatibility
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"

# Step A: Perform alignment and save as TEXT (SAM)
# We do NOT use the pipe '|' here so the .sam file stays on the disk
bwa mem -t $THREADS -R "$RG" $REFERENCE ${SAMPLE}_R1_trimmed.fastq ${SAMPLE}_R2_trimmed.fastq > ${SAMPLE}.sam
echo "[INFO] Raw alignment saved to ${SAMPLE}.sam (Open with 'less -S' to inspect)"

echo "========================================="
echo "[INFO] Switching to BAM environment: $ENV_BAM"
conda activate $ENV_BAM
echo "========================================="

echo "[STEP 6.1/17] Converting SAM to binary BAM..."
# Step B: Manual conversion using samtools view
# NO 'rm' command follows — the .sam file remains in the folder
samtools view -@ $THREADS -b ${SAMPLE}.sam -o ${SAMPLE}.bam

echo "[STEP 7/17] Sorting BAM..."
samtools sort -m 1G -@ $THREADS -o ${SAMPLE}_sorted.bam ${SAMPLE}.bam

echo "[STEP 8/17] Indexing sorted BAM..."
samtools index ${SAMPLE}_sorted.bam

echo "[STEP 9/17] Name-sorting for fixmate..."
samtools sort -n -m 1G -@ $THREADS -o ${SAMPLE}_namesorted.bam ${SAMPLE}_sorted.bam

echo "[STEP 10-11/17] Fixmate and Re-sorting..."
samtools fixmate -m ${SAMPLE}_namesorted.bam ${SAMPLE}_fixmate.bam
samtools sort -m 1G -@ $THREADS -o ${SAMPLE}_fixmate_sorted.bam ${SAMPLE}_fixmate.bam

echo "[STEP 12/17] Marking Duplicates..."
samtools markdup -@ $THREADS ${SAMPLE}_fixmate_sorted.bam ${SAMPLE}_marked.bam

echo "[STEP 13/17] Filtering: Removing Duplicates..."
samtools view -b -F 1024 -@ $THREADS ${SAMPLE}_marked.bam -o ${SAMPLE}_noduplicates.bam

echo "[STEP 14/17] Filtering: Removing Unmapped reads..."
samtools view -b -F 4 -@ $THREADS ${SAMPLE}_noduplicates.bam -o ${SAMPLE}_mapped.bam

echo "[STEP 15/17] Final BAM Indexing..."
samtools index ${SAMPLE}_mapped.bam

echo "[STEP 16/17] Variant Calling..."
mkdir -p variants
bcftools mpileup -f $REFERENCE ${SAMPLE}_mapped.bam | bcftools call -mv -Ov -o variants/${SAMPLE}_variants.vcf

echo "[STEP 17/17] VCF Statistics..."
bcftools stats variants/${SAMPLE}_variants.vcf > variants/${SAMPLE}_variants.stats

echo "========================================="
echo "EDUCATIONAL PIPELINE FINISHED"
echo "[NOTE] All intermediate files (SAM, temporary BAMs) are preserved for inspection."
echo "========================================="
