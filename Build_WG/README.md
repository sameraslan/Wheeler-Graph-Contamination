# Building a Wheeler Graph for Contamination

---

Overview
This folder contains the files for building a wheeler graph, specifically a De Bruijn graph from one or multiple FASTA files. In our specific contamination use case, a common scenario would be to have an input file set of one target genome FASTA and one or multiple contaminant genome FASTA files. The resultant wheeler graph is constructed from these multiple target and contaminant FASTA files, where each edge in the graph is labeled with a nucleotide, along with a set of corresponding organisms that the label came from. For example, the notation 'S10 -> S1 [ label = "A (0,1,3)" ];' indicates that the edge from Node 10 to Node 1 has the label 'A' and is associated with organisms 0, 1, and 3.

---

Prerequisites

- Python installed on your system.
- Biopython library, installable via `pip install biopython`.
- Graphviz library (for optional visualization), installable via `pip install graphviz`.

---

Running the Script

Use the following command structure in your terminal:

python generate_debruijn.py [OPTIONS] FASTA_FILES...

Replace 'FASTA_FILES...' with the paths to your FASTA files.

Command-Line Arguments

- -h / --help: Display usage information. Use this for help with options.
- -o / --ofile FILE: Specify the output file (in DOT format) for the Wheeler graph. Replace 'FILE' with your output filename.
- -k / --kmer LENGTH: Set the k-mer length for graph construction. K-mers are substrings of length 'LENGTH'.
- -l / --seqLen LENGTH: Define the sequence length to be considered from each FASTA entry. If 'LENGTH' is specified, only the first 'LENGTH' nucleotides of each sequence are used. If not specified or set to '-1', the full sequence length is used.
- -a / --alnNum NUMBER: Specify the number of sequences to process from the FASTA file. Only the first 'NUMBER' sequences will be considered.

---

Example Usage:
`python generate_debruijn.py -o output.dot -k 5 -l 100 -a 3 target.fasta contaminant1.fasta contaminant2.fasta`

In this command:

- 'output.dot' is the file where the Wheeler graph will be saved.
- K-mers of length 5 will be used.
- The first 100 nucleotides of each sequence will be considered.
- The script will process the first 3 sequences from each FASTA file.

---

Visualization

After running the script, use Graphviz to view and analyze the constructed Wheeler graph. This visualization will help in understanding the relationships and overlaps between the target and contaminant genomes.

Example Usage:
`dot -Tpng output.dot -o graph.png`
