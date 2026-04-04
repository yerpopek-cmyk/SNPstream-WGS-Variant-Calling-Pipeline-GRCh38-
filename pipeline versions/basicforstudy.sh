#!/bin/bash
# Stop the script if any critical error occurs
set -e  

# --- PARAMETER SETTINGS ---
ACCESSION="SRR622461"
REFERENCE="GCA_000001405.29.fasta"
SAMPLE="sample"
THREADS=6

# --- CONDA ENVIRONMENT SETTINGS ---
ENV_PREP="ngs_prep"       # Environment: sra-tools, fastqc, fastp, bwa
ENV_BAM="ngs_bam"         # Environment: samtools, bcftools
# -----------------

echo "========================================="
echo "[INFO] Activating preparation environment: $ENV_PREP"
eval "$(conda shell.bash hook)"
conda activate $ENV_PREP
echo "========================================="

echo "[STEP 1/16] Downloading data from SRA..."
prefetch $ACCESSION
fasterq-dump $ACCESSION --split-files --threads $THREADS
mv ${ACCESSION}_1.fastq ${SAMPLE}_R1.fastq
mv ${ACCESSION}_2.fastq ${SAMPLE}_R2.fastq

echo "[STEP 2-4/16] Quality Control and Cleaning (fastp)..."
mkdir -p fastqc_before fastp_reports fastqc_after
fastqc ${SAMPLE}_R1.fastq ${SAMPLE}_R2.fastq -t $THREADS -o fastqc_before/
fastp -i ${SAMPLE}_R1.fastq -I ${SAMPLE}_R2.fastq \
      -o ${SAMPLE}_R1_trimmed.fastq -O ${SAMPLE}_R2_trimmed.fastq \
      -h fastp_reports/fastp_report.html -j fastp_reports/fastp_report.json --thread $THREADS
fastqc ${SAMPLE}_R1_trimmed.fastq ${SAMPLE}_R2_trimmed.fastq -t $THREADS -o fastqc_after/

echo "[STEP 5/16] Indexing Reference Genome..."
if [ ! -f "${REFERENCE}.bwt" ]; then bwa index $REFERENCE; fi

echo "[STEP 6/16] ALIGNMENT: Creating text-based SAM file..."
# Read Group (crucial for sample identification in downstream analysis)
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"

# We intentionally DO NOT use a pipe '|' to keep the large .sam on disk for educational review
bwa mem -t $THREADS -R "$RG" $REFERENCE ${SAMPLE}_R1_trimmed.fastq ${SAMPLE}_R2_trimmed.fastq > ${SAMPLE}.sam
echo "[INFO] Text-based alignment file created: ${SAMPLE}.sam"

echo "========================================="
echo "[INFO] Switching to BAM environment: $ENV_BAM"
conda activate $ENV_BAM
echo "========================================="

echo "[STEP 6.1/16] CONVERSION: Converting SAM to binary BAM..."
# Here the student can see how the text "giant" turns into a compact binary file
# We are NOT removing the original .sam file with 'rm'
samtools view -@ $THREADS -b ${SAMPLE}.sam -o ${SAMPLE}.bam
echo "[INFO] Binary file created: ${SAMPLE}.bam"

echo "[STEP 7/16] Sorting BAM by coordinates..."
samtools sort -m 1G -@ $THREADS -o ${SAMPLE}_sorted.bam ${SAMPLE}.bam

echo "[STEP 8/16] Indexing sorted BAM..."
samtools index ${SAMPLE}_sorted.bam

echo "[STEP 9/16] Sorting by name (required for duplicate marking)..."
samtools sort -n -m 1G -@ $THREADS -o ${SAMPLE}_namesorted.bam ${SAMPLE}_sorted.bam

echo "[STEP 10/16] Adding fixmate tags..."
samtools fixmate -m ${SAMPLE}_namesorted.bam ${SAMPLE}_fixmate.bam

echo "[STEP 11/16] Re-sorting by coordinates..."
samtools sort -m 1G -@ $THREADS -o ${SAMPLE}_fixmate_sorted.bam ${SAMPLE}_fixmate.bam

echo "[STEP 12/16] Marking duplicates (markdup)..."
samtools markdup -@ $THREADS ${SAMPLE}_fixmate_sorted.bam ${SAMPLE}_marked.bam

echo "[STEP 13-14/16] Filtering: Removing duplicates and unmapped reads..."
samtools view -b -F 1024 -@ $THREADS ${SAMPLE}_marked.bam -o ${SAMPLE}_noduplicates.bam
samtools view -b -F 4 -@ $THREADS ${SAMPLE}_noduplicates.bam -o ${SAMPLE}_mapped.bam

echo "[STEP 15/16] Final BAM indexing..."
samtools index ${SAMPLE}_mapped.bam

echo "[STEP 16/16] Variant Calling..."
mkdir -p variants
if bcftools mpileup -f $REFERENCE ${SAMPLE}_mapped.bam | bcftools call -mv -Ov -o variants/${SAMPLE}_variants.vcf; then
    echo "[INFO] Variant calling completed successfully."
else
    echo "[WARNING] Falling back to system /usr/bin/bcftools due to library error..."
    /usr/bin/bcftools mpileup -f $REFERENCE ${SAMPLE}_mapped.bam | /usr/bin/bcftools call -mv -Ov -o variants/${SAMPLE}_variants.vcf
fi

echo "========================================="
echo "EDUCATIONAL PIPELINE COMPLETED 🎉"
echo "[IMPORTANT] All intermediate files (SAM and temporary BAMs) have been saved."
echo "Check file size differences using: ls -lh ${SAMPLE}*"
echo "========================================="
