#!/bin/bash
# Stop the script on any critical error
set -e  

# --- PIPELINE PARAMETERS ---
ACCESSION="SRR622461"
REFERENCE="GCA_000001405.29.fasta"
SAMPLE="sample"
THREADS=6

# --- CONDA ENVIRONMENTS ---
# Enter the exact names of your created Conda environments here
ENV_PREP="ngs_prep"       # Environment for sra-tools, fastqc, fastp, bwa
ENV_BAM="ngs_bam"         # Separate clean environment for samtools and bcftools
# -----------------

echo "========================================="
echo "[INFO] Initializing Conda and activating prep environment: $ENV_PREP"
eval "$(conda shell.bash hook)"
conda activate $ENV_PREP
echo "========================================="

echo "[STEP 1/16] Downloading raw reads from SRA (prefetch + fasterq-dump)..."
# prefetch downloads the .sra file locally first — more reliable than direct fasterq-dump
if prefetch $ACCESSION; then
    echo "[INFO] Prefetch completed successfully."
else
    echo "[WARNING] Prefetch failed. Retrying in 60 seconds..."
    sleep 60
    prefetch $ACCESSION
fi

# Convert downloaded .sra file to paired fastq files
fasterq-dump $ACCESSION --split-files --threads $THREADS
mv ${ACCESSION}_1.fastq ${SAMPLE}_R1.fastq
mv ${ACCESSION}_2.fastq ${SAMPLE}_R2.fastq

echo "[STEP 2/16] Quality control of raw reads (FastQC, pre-trimming)..."
mkdir -p fastqc_before
fastqc ${SAMPLE}_R1.fastq ${SAMPLE}_R2.fastq -t $THREADS -o fastqc_before/

echo "[STEP 3/16] Filtering and adapter trimming (fastp)..."
mkdir -p fastp_reports
fastp \
    -i ${SAMPLE}_R1.fastq -I ${SAMPLE}_R2.fastq \
    -o ${SAMPLE}_R1_trimmed.fastq -O ${SAMPLE}_R2_trimmed.fastq \
    -h fastp_reports/fastp_report.html \
    -j fastp_reports/fastp_report.json \
    --thread $THREADS

echo "[STEP 4/16] Quality control of trimmed reads (FastQC, post-trimming)..."
mkdir -p fastqc_after
fastqc ${SAMPLE}_R1_trimmed.fastq ${SAMPLE}_R2_trimmed.fastq -t $THREADS -o fastqc_after/

echo "[STEP 5/16] Indexing reference genome (BWA index)..."
bwa index $REFERENCE

echo "[STEP 6/16] Aligning reads to reference and converting to BAM (BWA MEM + samtools view)..."
if bwa mem -t $THREADS $REFERENCE ${SAMPLE}_R1_trimmed.fastq ${SAMPLE}_R2_trimmed.fastq | samtools view -@ $THREADS -b - > ${SAMPLE}.bam; then
    echo "[INFO] Alignment completed successfully."
else
    echo "[WARNING] Alignment failed. Retrying with 2 threads to reduce memory load..."
    bwa mem -t 2 $REFERENCE ${SAMPLE}_R1_trimmed.fastq ${SAMPLE}_R2_trimmed.fastq | samtools view -b - > ${SAMPLE}.bam
fi

echo "========================================="
echo "[INFO] Switching Conda environment to: $ENV_BAM (for samtools/bcftools)..."
conda activate $ENV_BAM
echo "========================================="

echo "[STEP 7/16] Sorting BAM by genomic coordinates..."
samtools sort -m 1G -@ $THREADS -o ${SAMPLE}_sorted.bam ${SAMPLE}.bam
rm -f ${SAMPLE}.bam

echo "[STEP 8/16] Indexing coordinate-sorted BAM..."
samtools index ${SAMPLE}_sorted.bam

echo "[STEP 9/16] Sorting BAM by read names (prep for fixmate)..."
samtools sort -n -m 1G -@ $THREADS -o ${SAMPLE}_namesorted.bam ${SAMPLE}_sorted.bam

echo "[STEP 10/16] Adding technical tags for duplicate detection (samtools fixmate)..."
samtools fixmate -m ${SAMPLE}_namesorted.bam ${SAMPLE}_fixmate.bam
rm -f ${SAMPLE}_namesorted.bam

echo "[STEP 11/16] Re-sorting by coordinates after fixmate..."
samtools sort -m 1G -@ $THREADS -o ${SAMPLE}_fixmate_sorted.bam ${SAMPLE}_fixmate.bam
rm -f ${SAMPLE}_fixmate.bam

echo "[STEP 12/16] Marking PCR duplicates (samtools markdup)..."
samtools markdup -@ $THREADS ${SAMPLE}_fixmate_sorted.bam ${SAMPLE}_marked.bam
rm -f ${SAMPLE}_fixmate_sorted.bam

echo "[STEP 13/16] Removing marked duplicates from BAM..."
samtools view -b -F 1024 -@ $THREADS ${SAMPLE}_marked.bam -o ${SAMPLE}_noduplicates.bam
rm -f ${SAMPLE}_marked.bam

echo "[STEP 14/16] Removing unmapped reads from BAM..."
samtools view -b -F 4 -@ $THREADS ${SAMPLE}_noduplicates.bam -o ${SAMPLE}_mapped.bam
rm -f ${SAMPLE}_noduplicates.bam

echo "[STEP 15/16] Indexing final mapped BAM..."
samtools index ${SAMPLE}_mapped.bam

echo "[STEP 16/16] Variant calling (bcftools mpileup | bcftools call)..."
mkdir -p variants
if bcftools mpileup -f $REFERENCE ${SAMPLE}_mapped.bam | bcftools call -mv -Ov -o variants/${SAMPLE}_variants.vcf; then
    echo "[INFO] Variant calling completed successfully (conda bcftools)."
else
    echo "[WARNING] Error detected, likely a library conflict (libcrypto). Switching to system /usr/bin/bcftools..."
    /usr/bin/bcftools mpileup -f $REFERENCE ${SAMPLE}_mapped.bam | /usr/bin/bcftools call -mv -Ov -o variants/${SAMPLE}_variants.vcf
fi

echo "========================================="
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "Final alignment : ${SAMPLE}_mapped.bam"
echo "Variant calls   : variants/${SAMPLE}_variants.vcf"
echo "========================================="
