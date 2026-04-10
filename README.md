# SNPstream-WGS-Variant-Calling-Pipeline-GRCh38-
# SNPstream-WGS-Variant-Calling-Pipeline (GRCh38)

Этот проект представляет собой автоматизированный пайплайн для анализа данных полногеномного секвенирования (WGS) с целью поиска герминальных генетических вариантов.This project provides an automated bioinformatics pipeline for analyzing Whole Genome Sequencing (WGS) data to identify germline genetic variants.

## Описание процесса|Process Overview
Пайплайн выполняет полный цикл обработки данных: от сырых чтений до фильтрованного списка вариантов. В качестве референса используется сборка генома человека **GRCh38**.This pipeline implements a complete end-to-end data processing cycle, transitioning from raw sequencing reads to a filtered list of genetic variants. The workflow is optimized for high-throughput genomic data and utilizes the GRCh38 (Human Genome Build 38) as the reference assembly.

### Основные этапы|Key Analytical Stages
1. **Quality Control:** Оценка качества сырых данных (FastQC). Initial assessment of raw data integrity using FastQC.
2. **Preprocessing:** Тримминг адаптеров и фильтрация низкокачественных чтений (fastp). Automated adapter trimming and low-quality base filtering via fastp.
3. **Alignment:** Выравнивание чтений на референсный геном (BWA-MEM). High-sensitivity mapping of processed reads to the reference genome using the BWA-MEM algorithm.
4. **Post-alignment:** Сортировка, индексация и разметка дубликатов (Samtools, Sambamba). Coordinate sorting, indexing, and duplicate marking using Samtools and Sambamba to ensure high-confidence calls.
5. **Variant Calling:** Поиск SNP и Indels (HaplotypeCaller / BCFtools). Identification of Single Nucleotide Polymorphisms (SNPs) and short insertions/deletions (Indels) using BCFtools (or HaplotypeCaller).
6. **Annotation:** Аннотация найденных вариантов. Functional annotation of identified variants to determine biological significance. (is not included in my pipeline)

## Инструменты|Required Tools
To run the pipeline, ensure the following software suite is installed:
Для работы пайплайна требуются установленные инструменты:
* `bwa`
* `samtools`
* `bcftools`
* `fastp` / `fastqc`

## Technical Standards & Features
Data Integrity: To ensure high-confidence results, the pipeline implements rigorous Phred Quality Score filtering. Only high-quality reads contribute to the final variant calls, significantly reducing the false discovery rate (FDR).

Performance & Scalability: The workflow is built with full multi-threading support for all computationally intensive steps (BWA, Samtools, fastp). This architecture significantly reduces the time required for WGS analysis on modern multi-core workstations and HPC environments.

Standardization: All generated output follows the VCF v4.2 standard. This ensures full compatibility with the IGV(Integrative Genomics Viewer), UCSC Genome Browser, and other downstream functional annotation tools.

## Как запустить|Getting Started
1. Clone the repository:
  ``bash
  git clone https://github.com/yerpopek-cmyk/SNPstream-WGS-Variant-Calling-Pipeline-GRCh38-.git
cd SNPstream-WGS-Variant-Calling-Pipeline-GRCh38-
2. Prepare your environment:
Ensure your reference genome (.fasta) and input SRA accessions are correctly configured in the script.
3. Execute the pipeline:
chmod +x pipeline.sh
./pipeline.sh                              ./"name".sh   |   bash "name".sh
