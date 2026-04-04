#!/bin/bash
# Stop the script on any critical error
set -e  

# --- PIPELINE PARAMETERS ---
ACCESSION="SRR622461"
REFERENCE="GCA_000001405.29.fasta"
SAMPLE="sample"
THREADS=6

# --- CONDA ENVIRONMENTS ---
# Replace these with your actual environment names
ENV_PREP="ngs_prep"       # Env for sra-tools, fastqc, fastp, bwa
ENV_BAM="ngs_bam"         # Env for samtools and bcftools
# -----------------

# [FUNCTION] Check if required tools are installed
check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo "[ERROR] Program $1 is not installed or not in PATH. Aborting."
        exit 1
    fi
}

echo "========================================="
echo "[INFO] Initializing Conda and activating PREP environment: $ENV_PREP"
eval "$(conda shell.bash hook)"
conda activate $ENV_PREP
echo "========================================="

echo "[PREPARATION] Checking required software in $ENV_PREP..."
for tool in prefetch fasterq-dump fastqc fastp bwa; do
    check_tool $tool
done

echo "[STEP 1/17] Downloading raw reads from SRA (prefetch + fasterq-dump)..."
if prefetch $ACCESSION; then
    echo "[INFO] Prefetch completed successfully."
else
    echo "[WARNING] Prefetch failed. Retrying in 60 seconds..."
    sleep 60
    prefetch $ACCESSION
fi

# Convert .sra to paired fastq files
fasterq-dump $ACCESSION --split-files --threads $THREADS
mv ${ACCESSION}_1.fastq ${SAMPLE}_R1.fastq
mv ${ACCESSION}_2.fastq ${SAMPLE}_R2.fastq

# [CHECK] Verify if FASTQ files are not empty
if [ -s "${SAMPLE}_R1.fastq" ] && [ -s "${SAMPLE}_R2.fastq" ]; then
    echo "[INFO] FASTQ files created successfully and contain data."
else
    echo "[ERROR] FASTQ files are empty or were not created. Aborting."
    exit 1
fi

echo "[STEP 2/17] Quality Control: Raw reads (FastQC, pre-trimming)..."
mkdir -p fastqc_before
fastqc ${SAMPLE}_R1.fastq ${SAMPLE}_R2.fastq -t $THREADS -o fastqc_before/

echo "[STEP 3/17] Filtering and Adapter Trimming (fastp)..."
mkdir -p fastp_reports
fastp \
    -i ${SAMPLE}_R1.fastq -I ${SAMPLE}_R2.fastq \
    -o ${SAMPLE}_R1_trimmed.fastq -O ${SAMPLE}_R2_trimmed.fastq \
    -h fastp_reports/fastp_report.html \
    -j fastp_reports/fastp_report.json \
    --thread $THREADS

echo "[STEP 4/17] Quality Control: Trimmed reads (FastQC, post-trimming)..."
mkdir -p fastqc_after
fastqc ${SAMPLE}_R1_trimmed.fastq ${SAMPLE}_R2_trimmed.fastq -t $THREADS -o fastqc_after/

echo "[STEP 5/17] Indexing Reference Genome (BWA index)..."
# Skip if index already exists to save time
if [ ! -f "${REFERENCE}.bwt" ]; then
    bwa index $REFERENCE
else
    echo "[INFO] Genome index already exists. Skipping step."
fi

echo "[STEP 6/17] Alignment to Reference (BWA MEM + samtools view)..."
# Define Read Group for downstream compatibility
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"

# Attempt 1: Full threads
if bwa mem -t $THREADS -R "$RG" $REFERENCE ${SAMPLE}_R1_trimmed.fastq ${SAMPLE}_R2_trimmed.fastq | samtools view -@ $THREADS -b - > ${SAMPLE}.bam; then
    echo "[INFO] Alignment completed successfully."
else
    echo "[WARNING] Alignment failed. Retrying with 2 threads to reduce system load..."
    # Attempt 2: Safety mode
    if bwa mem -t 2 -R "$RG" $REFERENCE ${SAMPLE}_R1_trimmed.fastq ${SAMPLE}_R2_trimmed.fastq | samtools view -b - > ${SAMPLE}.bam; then
        echo "[INFO] Alignment successful on second attempt."
    else
        echo "[ERROR] Alignment failed again. Check memory and reference files."
        exit 1
    fi
fi

echo "========================================="
echo "[INFO] Switching Conda environment to: $ENV_BAM"
conda activate $ENV_BAM
echo "========================================="

echo "[STEP 7/17] Sorting BAM by genomic coordinates..."
samtools sort -m 1G -@ $THREADS -o ${SAMPLE}_sorted.bam ${SAMPLE}.bam
rm -f ${SAMPLE}.bam

echo "[STEP 8/17] Indexing coordinate-sorted BAM..."
samtools index ${SAMPLE}_sorted.bam

echo "[STEP 9/17] Sorting BAM by read names (prep for fixmate)..."
samtools sort -n -m 1G -@ $THREADS -o ${SAMPLE}_namesorted.bam ${SAMPLE}_sorted.bam

echo "[STEP 10/17] Adding technical tags for duplicates (samtools fixmate)..."
samtools fixmate -m ${SAMPLE}_namesorted.bam ${SAMPLE}_fixmate.bam
rm -f ${SAMPLE}_namesorted.bam

echo "[STEP 11/17] Re-sorting by coordinates after fixmate..."
samtools sort -m 1G -@ $THREADS -o ${SAMPLE}_fixmate_sorted.bam ${SAMPLE}_fixmate.bam
rm -f ${SAMPLE}_fixmate.bam

echo "[STEP 12/17] Marking PCR duplicates (samtools markdup)..."
samtools markdup -@ $THREADS ${SAMPLE}_fixmate_sorted.bam ${SAMPLE}_marked.bam
rm -f ${SAMPLE}_fixmate_sorted.bam

echo "[STEP 13/17] Removing marked duplicates..."
samtools view -b -F 1024 -@ $THREADS ${SAMPLE}_marked.bam -o ${SAMPLE}_noduplicates.bam
rm -f ${SAMPLE}_marked.bam

echo "[STEP 14/17] Removing unmapped reads..."
samtools view -b -F 4 -@ $THREADS ${SAMPLE}_noduplicates.bam -o ${SAMPLE}_mapped.bam
rm -f ${SAMPLE}_noduplicates.bam

echo "[STEP 15/17] Indexing final mapped BAM..."
samtools index ${SAMPLE}_mapped.bam

echo "[STEP 16/17] Variant Calling (bcftools mpileup | bcftools call)..."
mkdir -p variants
if bcftools mpileup -f $REFERENCE ${SAMPLE}_mapped.bam | bcftools call -mv -Ov -o variants/${SAMPLE}_variants.vcf; then
    echo "[INFO] Variant calling completed successfully (conda bcftools)."
else
    echo "[WARNING] Error detected (likely libcrypto conflict). Switching to system /usr/bin/bcftools..."
    /usr/bin/bcftools mpileup -f $REFERENCE ${SAMPLE}_mapped.bam | /usr/bin/bcftools call -mv -Ov -o variants/${SAMPLE}_variants.vcf
fi

echo "[STEP 17/17] Collecting VCF statistics (bcftools stats)..."
bcftools stats variants/${SAMPLE}_variants.vcf > variants/${SAMPLE}_variants.stats

echo "========================================="
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "Final alignment  : ${SAMPLE}_mapped.bam"
echo "Variant calls    : variants/${SAMPLE}_variants.vcf"
echo "VCF Statistics   : variants/${SAMPLE}_variants.stats"
echo "========================================="
