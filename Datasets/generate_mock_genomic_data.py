import random
from Bio import SeqIO

def generate_sequence(length, gc_content, motifs=None):
    """
    Generates a random DNA sequence with a specified GC content and optional motifs.
    
    :param length: Length of the DNA sequence to generate.
    :param gc_content: Proportion of the sequence that should be G or C nucleotides.
    :param motifs: List of specific motifs (subsequences) to be included in the sequence.
    :return: A string representing the generated DNA sequence.
    """
    # Calculate the number of GC and AT pairs based on the desired GC content
    gc_pairs = int(length * gc_content)
    at_pairs = length - gc_pairs

    # Randomly generate the GC and AT parts of the sequence
    sequence = random.choices(['G', 'C'], k=gc_pairs) + random.choices(['A', 'T'], k=at_pairs)
    random.shuffle(sequence)

    # Insert motifs into the sequence at random positions, if provided
    if motifs:
        for motif in motifs:
            insert_position = random.randint(0, len(sequence) - len(motif))
            sequence[insert_position:insert_position + len(motif)] = motif

    return ''.join(sequence)

def generate_fasta(file_name, sequence_length, gc_content, motifs=None, num_sequences=1):
    """
    Generates a FASTA file containing specified number of DNA sequences.

    :param file_name: Name of the output FASTA file.
    :param sequence_length: Length of each DNA sequence.
    :param gc_content: Proportion of G or C nucleotides in each sequence.
    :param motifs: Optional motifs to include in the sequences.
    :param num_sequences: Number of sequences to generate.
    """
    with open(file_name, 'w') as f:
        for i in range(num_sequences):
            sequence = generate_sequence(sequence_length, gc_content, motifs)
            f.write(f'>seq{i}\n{sequence}\n')

def generate_fastq(fasta_files, labels, output_file, read_length, num_reads_per_file):
    """
    Generates a FASTQ file from a list of FASTA files, simulating sequencing reads.

    :param fasta_files: List of input FASTA files.
    :param labels: Labels to prepend to the read IDs, indicating the source file.
    :param output_file: Name of the output FASTQ file.
    :param read_length: Length of each simulated read.
    :param num_reads_per_file: Number of reads to simulate from each FASTA file.
    """
    with open(output_file, 'w') as of:
        for fasta_file, label in zip(fasta_files, labels):
            for record in SeqIO.parse(fasta_file, "fasta"):
                full_seq = str(record.seq)
                for _ in range(num_reads_per_file):
                    # Simulate a sequencing read of specified length from the sequence
                    start = random.randint(0, len(full_seq) - read_length)
                    seq = full_seq[start:start + read_length]

                    # Generate a mock quality score string
                    qual = ''.join(random.choices(['!', '"', '#', '$', '%', '&', "'"], k=read_length))
                    of.write(f'@{label}|{record.id}|pos{start}\n{seq}\n+\n{qual}\n')

# Define unique motifs and different GC contents for each type of sequence (host and contaminants)
motifs = {
    "host": ["ATGCGTAC", "GATTACAGAT"],
    "contaminant1": ["GGGGGG", "CCCCCC"],
    "contaminant2": ["AAAAAAAAAA", "TTTTTTTTTT"]
}
gc_contents = {
    "host": 0.5,
    "contaminant1": 0.9,
    "contaminant2": 0.1
}

fasta_files = ['host.fasta', 'contaminant1.fasta', 'contaminant2.fasta']
sequence_lengths = [10000, 100, 100]
labels = ['host', 'contaminant1', 'contaminant2']

# Generate FASTA files for host and contaminants with specified parameters
for file_name, seq_len, label in zip(fasta_files, sequence_lengths, labels):
    generate_fasta(file_name, seq_len, gc_contents[label], motifs[label], 1)

# Generate a FASTQ file simulating sequencing reads from the generated FASTA file(s)
generate_fastq(fasta_files, labels, 'mock_data.fastq', 50, 500)  # Read length 50, 500 reads per fasta file