import random
from Bio import SeqIO

def generate_fasta(file_name, sequence_length, num_sequences=1):
    with open(file_name, 'w') as f:
        for i in range(num_sequences):
            sequence = ''.join(random.choices(['A', 'C', 'G', 'T'], k=sequence_length))
            f.write(f'>seq{i}\n{sequence}\n')

def generate_fastq(fasta_files, labels, output_file, read_length, num_reads_per_file):
    with open(output_file, 'w') as of:
        for fasta_file, label in zip(fasta_files, labels):
            for record in SeqIO.parse(fasta_file, "fasta"):
                full_seq = str(record.seq)
                for _ in range(num_reads_per_file):
                    start = random.randint(0, len(full_seq) - read_length)
                    seq = full_seq[start:start + read_length]  # Random subsequence
                    qual = ''.join(random.choices(['!', '"', '#', '$', '%', '&', "'"], k=read_length))  # Mock quality score
                    of.write(f'@{label}|{record.id}|pos{start}\n{seq}\n+\n{qual}\n')

fasta_files = ['host.fasta', 'contaminant1.fasta', 'contaminant2.fasta']
sequence_lengths = [500, 500, 500]

for file_name, seq_len in zip(fasta_files, sequence_lengths):
    generate_fasta(file_name, seq_len, 1)

labels = ['host', 'contaminant1', 'contaminant2']
generate_fastq(fasta_files, labels, 'mock_data.fastq', 100, 50)  # Read length 100, 50 reads per fasta file