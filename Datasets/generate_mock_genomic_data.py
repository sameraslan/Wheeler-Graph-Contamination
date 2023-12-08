import random
from Bio import SeqIO

def generate_sequence(length, gc_content, motifs=None):
    gc_pairs = int(length * gc_content)
    at_pairs = length - gc_pairs
    sequence = random.choices(['G', 'C'], k=gc_pairs) + random.choices(['A', 'T'], k=at_pairs)
    random.shuffle(sequence)

    if motifs:
        for motif in motifs:
            insert_position = random.randint(0, len(sequence) - len(motif))
            sequence[insert_position:insert_position+len(motif)] = motif

    return ''.join(sequence)

def generate_fasta(file_name, sequence_length, gc_content, motifs=None, num_sequences=1):
    with open(file_name, 'w') as f:
        for i in range(num_sequences):
            sequence = generate_sequence(sequence_length, gc_content, motifs)
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

# Define motifs and GC content for each type of sequence
motifs = {
    "host": ["ATGC", "GATTACA"],
    "contaminant1": ["CGT", "TGCA"],
    "contaminant2": ["AAAA", "TTTT"]
}
gc_contents = {
    "host": 0.4,  # 40% GC
    "contaminant1": 0.6,  # 60% GC
    "contaminant2": 0.3  # 30% GC
}

fasta_files = ['host.fasta', 'contaminant1.fasta', 'contaminant2.fasta']
sequence_lengths = [10000, 100, 100]
labels = ['host', 'contaminant1', 'contaminant2']

# Generate FASTA files
for file_name, seq_len, label in zip(fasta_files, sequence_lengths, labels):
    generate_fasta(file_name, seq_len, gc_contents[label], motifs[label], 1)

# Generate FASTQ file
generate_fastq(fasta_files, labels, 'mock_data.fastq', 50, 500)  # Read length 50, 500 reads per fasta file