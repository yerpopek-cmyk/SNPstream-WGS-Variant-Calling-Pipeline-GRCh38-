#!/bin/bash
# Stop the script on any error
set -e  

# --- CONFIGURATION & PARAMETERS ---
ACCESSION="SRR622461"
REFERENCE="GCA_000001405.29.fasta"
SAMPLE="sample"
THREADS=6
ENV_PREP="ngs_prep"
ENV_BAM="ngs_bam"

# [FUNCTION] Initialize Conda for script usage
# This allows 'conda activate' to work inside the bash script
eval "$(conda shell.bash hook)"

# Trap for clean exit on interruption
trap "echo -e '\n[ERROR] Pipeline interrupted. Exiting...'; exit" INT TERM

echo "========================================="
echo "[INFO] Starting WGS Variant Calling Pipeline"
echo "[INFO] Target Accession: $ACCESSION"
echo "========================================="

# STEP 1-4: Data Acquisition and Quality Control
conda activate $ENV_PREP

echo "[STEP 1/17] Downloading raw reads from SRA..."
prefetch $ACCESSION && fasterq-dump $ACCESSION --split-files --threads $THREADS
mv ${ACCESSION}_1.fastq ${SAMPLE}_R1.fastq
mv ${ACCESSION}_2.fastq ${SAMPLE}_R2.fastq

echo "[STEP 2-4/17] Quality Control and Adapter Trimming..."
mkdir -p fastqc_before fastp_reports fastqc_after
fastqc ${SAMPLE}_R1.fastq ${SAMPLE}_R2.fastq -t $THREADS -o fastqc_before/
fastp -i ${SAMPLE}_R1.fastq -I ${SAMPLE}_R2.fastq \
      -o ${SAMPLE}_R1_trimmed.fastq -O ${SAMPLE}_R2_trimmed.fastq \
      --thread $THREADS \
      -h fastp_reports/fastp_report.html -j fastp_reports/fastp_report.json

# STEP 5: Reference Indexing
if [ ! -f "${REFERENCE}.bwt" ]; then 
    echo "[STEP 5/17] Indexing Reference Genome (This may take several minutes)..."
    bwa index $REFERENCE
fi

# STEP 6: Alignment (Educational: Creating a RAW SAM file)
echo "[STEP 6/17] Alignment: Mapping reads to reference (Creating SAM)..."
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"
bwa mem -t $THREADS -R "$RG" $REFERENCE ${SAMPLE}_R1_trimmed.fastq ${SAMPLE}_R2_trimmed.fastq > ${SAMPLE}.sam

echo "========================================="
echo "[INFO] Switching to Analysis Environment: $ENV_BAM"
conda activate $ENV_BAM
echo "========================================="

# STEP 6.1-15: BAM Processing (No file deletion for educational purposes)
echo "[STEP 6.1/17] Converting SAM (text) to BAM (binary)..."
samtools view -@ $THREADS -b ${SAMPLE}.sam -o ${SAMPLE}.bam

echo "[STEP 7-15/17] Post-alignment processing: Sorting, Marking Duplicates, and Filtering..."
# Sorting by coordinates
samtools sort -m 1G -@ $THREADS -o ${SAMPLE}_sorted.bam ${SAMPLE}.bam
samtools index ${SAMPLE}_sorted.bam

# Duplicate marking workflow
samtools sort -n -m 1G -@ $THREADS -o ${SAMPLE}_namesorted.bam ${SAMPLE}_sorted.bam
samtools fixmate -m ${SAMPLE}_namesorted.bam ${SAMPLE}_fixmate.bam
samtools sort -m 1G -@ $THREADS -o ${SAMPLE}_fixmate_sorted.bam ${SAMPLE}_fixmate.bam
samtools markdup -@ $THREADS ${SAMPLE}_fixmate_sorted.bam ${SAMPLE}_marked.bam

# Filtering: Remove duplicates (-F 1024) and unmapped reads (-F 4)
samtools view -b -F 1024 -@ $THREADS ${SAMPLE}_marked.bam -o ${SAMPLE}_noduplicates.bam
samtools view -b -F 4 -@ $THREADS ${SAMPLE}_noduplicates.bam -o ${SAMPLE}_mapped.bam
samtools index ${SAMPLE}_mapped.bam

# STEP 16-17: Variant Calling and Statistics
echo "[STEP 16/17] Variant Calling: Generating VCF..."
mkdir -p variants
# Note: ngs_bam environment includes openssl=1.0 to avoid libcrypto compatibility issues
bcftools mpileup -f $REFERENCE ${SAMPLE}_mapped.bam | bcftools call -mv -Ov -o variants/${SAMPLE}_variants.vcf

echo "[STEP 17/17] Generating VCF summary statistics..."
bcftools stats variants/${SAMPLE}_variants.vcf > variants/${SAMPLE}_variants.stats

echo "========================================="
echo "   PIPELINE FINISHED SUCCESSFULLY! 🎉    "
echo "========================================="
echo "[INFO] Final variants located in: /variants/"
echo "[INFO] Intermediate SAM/BAM files kept for review."
echo "========================================="
