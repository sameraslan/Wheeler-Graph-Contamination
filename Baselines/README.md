# Simulating and Analyzing Mycoplasma Contamination in Sequencing Data

# Overview

This repository includes a script for simulating single-end reads from human and Mycoplasma genomes, mixing them, and then analyzing the mixed reads for Mycoplasma contamination. The process involves simulating reads using the `art_illumina` tool, mixing these reads with a custom Python script, and then using Bowtie2 and Samtools for alignment and analysis. This script is ideal for studying cross-contamination scenarios in sequencing experiments.

---

# Prerequisites

- ART Illumina (a read simulator for next-generation sequencing), installable from [ART Illuminas website](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). More details:

  1. **Visit the ART Illumina Website:** Go to [the ART webpage](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm).
  2. **Download:** Find the appropriate version for your operating system (Linux or macOS) and download the compressed file.
  3. **Extract:** Once downloaded, extract the files. You can use a graphical file manager or use a terminal command like `tar -xvf art_illumina_*.tar.gz`, replacing `*` with the version number.
  4. **Access ART Illumina:** The extracted folder should contain the executable files for ART Illumina. You may need to give execute permissions to the binaries using a command like `chmod +x art_illumina`.
- Bowtie2, a tool for aligning sequencing reads, installable via `sudo apt-get install bowtie2` or through the [Bowtie 2 website](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
- Samtools, a suite of programs for interacting with high-throughput sequencing data, installable via `sudo apt-get install samtools`.

---

# Running the Script

Execute the script in your terminal:

```bash
bash ContaminationAnalysisPipeline.sh
```

# The Script

What the script does:

1. Simulating reads for both the Human and Mycoplasma genomes using the ART Illumina tool.
2. Mixing these simulated reads using MixingTheReads.py.
3. Building Bowtie2 indices for both genomes.
4. Aligning the mixed reads to the Mycoplasma genome using Bowtie2.
5. Converting the alignment results from SAM to BAM format using Samtools.
6. Counting the number of reads aligned to Mycoplasma and the total number of reads in the mixed FASTQ file for analysis.
