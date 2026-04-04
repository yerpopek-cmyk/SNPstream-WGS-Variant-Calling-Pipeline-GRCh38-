# Scientific Rationale of the Pipeline (Student Guide)

Hello! If you are reading this, it means you are ready to tear down the "engine" of NGS (Next-Generation Sequencing) analysis. Our pipeline consists of 17 steps. Here is why each one of them is essential:

### Block 1: Data Preparation
* **Step 1 (SRA Download):** We download the "raw" sequencing data. `prefetch` ensures a reliable download, while `fasterq-dump` converts the data into FASTQ files, which are the standard format for sequence data.
* **Steps 2-4 (QC & Trimming):** We assess data quality using `FastQC`. If there is "junk" or sequencer adapters at the ends of the reads, `fastp` trims them off. Think of this as washing your vegetables before cooking.

### Block 2: Alignment (Mapping)
* **Step 5 (Index):** We create a genome "index." This allows the software to find the correct location for a read in fractions of a second rather than years.
* **Step 6 (Alignment):** The `bwa mem` program searches for the best fit for each read on the reference genome.
    * **IMPORTANT:** This results in a **.SAM** file (Sequence Alignment Map). It is a massive text file that you can actually read with your own eyes using the `less -S` command.

### Block 3: BAM Processing (Clean-up)
* **Step 6.1 (Convert):** We convert the SAM into a **BAM** file (Binary). It contains the same information but is compressed to save disk space.
* **Steps 7-15 (Sorting & Markdup):**
    * We **Sort** the reads in order (like numbered pages in a book).
    * We identify **duplicates** (PCR clones) and mark them. If we don't ignore these, they can lead to "false positive" mutation calls.
    * We remove "garbage" reads that failed to map anywhere on the genome.

### Block 4: Final Stage (Variant Calling)
* **Step 16 (BCFtools):** The program looks at every "column" of aligned reads. If the reference genome has an "A" but our data consistently shows a "T" — we have found a mutation (variant)!
* **Step 17 (Stats):** We calculate the total number of findings to ensure the results are biologically logical and consistent.

**Remember:** In this educational pipeline, we **DO NOT delete** intermediate files. Compare their sizes using the `ls -lh` command to witness the magic of bioinformatics data compression!
