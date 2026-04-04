# Germline WGS Variant Calling Pipeline (GRCh38)

This is an automated bioinformatics pipeline for identifying germline genetic variants from Whole Genome Sequencing (WGS) data. Optimized for **GRCh38** and designed with an educational focus.

## 🚀 Quick Start (Installation)

To avoid library conflicts (especially the `libcrypto` error), follow these exact steps to set up your environments.

### 1. Clone the repository
```bash
git clone [https://github.com/yerpopek-cmyk/Germline-WGS-Variant-Calling-Pipeline-GRCh38-.git](https://github.com/yerpopek-cmyk/Germline-WGS-Variant-Calling-Pipeline-GRCh38-.git)
cd Germline-WGS-Variant-Calling-Pipeline-GRCh38-

# Environment for data prep and alignment
conda create -n ngs_prep -c bioconda -c conda-forge sra-tools fastqc fastp bwa -y

# Environment for BAM processing and Variant Calling (includes the libcrypto fix)
conda create -n ngs_bam -c bioconda -c conda-forge samtools bcftools openssl=1.0 -y

curl -L ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p13/GCA_000001405.29_GRCh38.p13_genomic.fna.gz | gunzip > GCA_000001405.29.fasta

chmod +x pipeline.sh
./pipeline.sh
