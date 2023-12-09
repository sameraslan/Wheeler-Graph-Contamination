# Mock Dataset Generation for Genomic Analysis

---

## Overview

This folder contains a Python script designed to generate a mock dataset for genomic analysis. The script creates random DNA sequences with specified GC content and motifs, saving them in FASTA format. It then simulates sequencing reads from these sequences, producing a FASTQ file. This is particularly useful for testing genomic algorithm/data pipelines such as ours.

---

## Prerequisites

If you are in the wheeler_graph_env then you don't need to worry about this. Otherwise, you need:

- Biopython library, installable via `pip install biopython`.

---

## Running the Script

To generate the mock dataset, run the script with the following command:

```bash
python generate_mock_genomic_data.py
```

The script does not require any command-line arguments as the parameters are predefined within the script itself.

---

## Script Functions

The script includes several key functions:

* `generate_sequence(length, gc_content, motifs)`: Generates a random DNA sequence with specified GC content and optional motifs.
* `generate_fasta(file_name, sequence_length, gc_content, motifs, num_sequences)`: Creates a FASTA file with a specified number of generated DNA sequences.
* `generate_fastq(fasta_files, labels, output_file, read_length, num_reads_per_file)`: Simulates sequencing reads from the given FASTA files, producing a FASTQ file.

---

## Output Files

The script generates the following files:

* FASTA files for each DNA source type (e.g., `host.fasta`, `contaminant1.fasta`, `contaminant2.fasta`).
* A FASTQ file (`mock_data.fastq`) containing simulated sequencing reads from the generated FASTA files.

---

## Customization

You can customize the script by modifying the motifs, GC contents, sequence lengths, and other parameters within the script to fit your specific requirements.
