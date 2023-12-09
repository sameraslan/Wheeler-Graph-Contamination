"""
This file takes in a list of FASTA files and outputs a De Bruijn graph in DOT format.
Adapted from `https://github.com/Kuanhao-Chao/Wheeler_Graph_Toolkit`.
Modified by Samer Aslan, Owen Bianchi, Anirudh Kashyap, and Senapati Diwangkara.
"""

import sys
import getopt
import os
from Bio import AlignIO
import node as n

queue = []  # Initialize a queue for breadth-first search (BFS)

# Usage message for command-line execution
USAGE = '''usage: DeBruijnGraph_generator.py [-h / --help] [-v/ --version] [-o / --ofile FILE] [-k / --kmer k-mer length] [-l / --seqLen sequence length] [-a / --alnNum alignment number] sequence FASTA files'''

# Main function to process command-line arguments and generate the graph
def main(argv):
    # Default values for the command-line arguments
    k = 3
    verbose = False
    outfile = "tmp.dot"
    seqLen = -1
    alnNum = 3

    # Parsing command-line arguments
    try:
        opts, args = getopt.getopt(argv,"hvo:k:l:a:",["ofile=", "kmer=", "seqLen=", "alnNum="])
    except getopt.GetoptError:
        print(USAGE)
        sys.exit(2)
    
    # Processing each option and argument
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(USAGE)
            sys.exit()
        if opt in ("-v", "--version"):
            verbose = True
        elif opt in ("-o", "--ofile"):
            outfile = arg
        elif opt in ("-k", "--kmer"):
            k = int(arg)
        elif opt in ("-l", "--seqLen"):
            seqLen = int(arg)
        elif opt in ("-a", "--alnNum"):
            alnNum = int(arg)

    # Check for at least one input FASTA file
    if len(args) < 1:
        print(USAGE)
        print("Please input at least one FASTA file")
        sys.exit(2)

    # Adjust k-mer length and initialize data structures
    k = k - 1
    kmer_dic = {}
    organism_id = 0  # Starting organism ID

    # Processing each FASTA file
    for fasta_file in args:
        process_fasta_file(fasta_file, k, seqLen, alnNum, organism_id, kmer_dic)
        organism_id += 1  # Increment organism ID for next file

    # Write the graph to the output file
    print("outfile: ", outfile)
    try:    
        os.remove(outfile)
    except OSError:
        pass
    fw = open(outfile, "a")
    fw.write("strict digraph  {\n")
    visited = []  # List for visited nodes for BFS
    bfs(visited, kmer_dic["$"], fw)
    fw.write("}\n")
    fw.close()

# Function to process a single FASTA file and populate the kmer dictionary
def process_fasta_file(fasta_file, k, seqLen, alnNum, organism_id, kmer_dic):
    # Reading alignment data from FASTA file
    alignment = AlignIO.read(fasta_file, "fasta")
    alignment = alignment[:alnNum]
    alignment_len = alignment.get_alignment_length()
    seq_number = len(alignment)
    kmer_set = set()  # Set to store unique k-mers

    # Processing each sequence in the alignment
    for row_idx in range(seq_number):
        # Removing gaps and slicing sequence if necessary
        row_seq = str(alignment[row_idx].seq)
        row_seq_parsed = row_seq.replace("-", "")
        if seqLen == -1 or seqLen > len(row_seq_parsed):
            lcl_seqLen = len(row_seq_parsed)
        else:
            lcl_seqLen = seqLen
        row_seq_parsed = row_seq_parsed[:lcl_seqLen] + "$"

        # Generating k-mers from the sequence
        for i in range(lcl_seqLen + 1):
            node_kmer = row_seq_parsed[i:i + k]
            kmer_set.add(node_kmer)

    # Creating nodes for each unique k-mer and adding to the dictionary
    node_idx = 0
    for node_kmer in kmer_set:
        if node_kmer in kmer_dic:
            kmer_dic[node_kmer].add_organism_id(organism_id)  # Add organism ID if node exists
        else:
            new_node = n.node(node_idx, node_kmer, alignment_len, {organism_id})
            node_idx += 1
            kmer_dic[node_kmer] = new_node

    # Building the De Bruijn graph by connecting nodes
    for row_idx in range(seq_number):
        row_seq = str(alignment[row_idx].seq)
        row_seq_parsed = row_seq.replace("-", "")
        if seqLen == -1 or seqLen > len(row_seq_parsed):
            lcl_seqLen = len(row_seq_parsed)
        else:
            lcl_seqLen = seqLen
        row_seq_parsed = row_seq_parsed[:lcl_seqLen] + "$"
        prev_node = ""
        for i in range(lcl_seqLen + 1):
            node_kmer = row_seq_parsed[i:i + k]
            curr_node = kmer_dic[node_kmer]
            if prev_node != "":
                curr_node.add_child(prev_node)
                prev_node.add_parent(curr_node)
            prev_node = curr_node

# Breadth-first search function to traverse the graph and write to file
def bfs(visited, node, fw):
    visited.append(node)
    queue.append(node)
    while queue:
        m = queue.pop(0)
        for m_child in m.children:
            # Writing graph edges to the output file
            organism_id_label = ','.join(map(str, m_child.organismID))
            edge_label = f"{m_child.nodelabel[0]} ({organism_id_label})"
            fw.write(f"\tS{m.nodeid} -> S{m_child.nodeid} [ label = \"{edge_label}\" ];\n")
        for child in m.children:
            if child not in visited:
                visited.append(child)
                queue.append(child)

if __name__ == "__main__":
    main(sys.argv[1:])
