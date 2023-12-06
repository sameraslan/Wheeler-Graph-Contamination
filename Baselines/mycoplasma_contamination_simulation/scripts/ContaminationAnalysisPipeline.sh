#!/bin/bash
data_path="/home/ak511/Desktop/WG-contamination/mycoplasma_contamination_simulation/Data/"
scripts_path="/home/ak511/Desktop/WG-contamination/mycoplasma_contamination_simulation/scripts"

# Script to simulate reads, mix them, and analyze for Mycoplasma contamination

# Simulate single-end reads from Human genome
# -ss HS25: Use the HiSeq 2500 sequencing system
# -i $data_path/Human.fasta: Input FASTA file for human genome
# -l 100: Generate reads of 100 bp
# -f 20: 20-fold coverage
# -o human_single_end_reads: Output prefix for human reads
art_illumina -ss HS25 -i "$data_path/Human.fasta" -l 100 -f 20 -o human_single_end_reads

# Simulate single-end reads from Mycoplasma genome
# -i $data_path/mycoplasma.fasta: Input FASTA file for Mycoplasma genome
# -o mycoplasma_reads: Output prefix for Mycoplasma reads
art_illumina -ss HS25 -i "$data_path/mycoplasma.fasta" -l 100 -f 20 -o mycoplasma_reads

# Mix the reads using the custom Python script
python "$scripts_path/MixingTheReads.py"

# Build Bowtie2 index for the Human genome
# $data_path/Human.fasta: Input FASTA file for human genome
# Human_index: Output prefix for the index
bowtie2-build "$data_path/Human.fasta" Human_index

# Build Bowtie2 index for the Mycoplasma genome
# $data_path/mycoplasma.fasta: Input FASTA file for Mycoplasma genome
# mycoplasma_index: Output prefix for the index
bowtie2-build "$data_path/mycoplasma.fasta" mycoplasma_index

# Align the mixed reads to the Mycoplasma genome
# -x mycoplasma_index: Use the Mycoplasma genome index
# -U mixed_reads.fq: Input FASTQ file of mixed reads
# -S aligned_to_mycoplasma.sam: Output SAM file
bowtie2 -x mycoplasma_index -U mixed_reads.fq -S aligned_to_mycoplasma.sam

# Convert SAM to BAM
samtools view -bS aligned_to_mycoplasma.sam > aligned_to_mycoplasma.bam

# Count the number of reads aligned to Mycoplasma
# -F 4 flag excludes unaligned reads
samtools view -c -F 4 aligned_to_mycoplasma.bam

# Count the total number of reads in the mixed FASTQ file
grep -c "^@" mixed_reads.fq
