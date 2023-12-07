"""Base code from KuanHao Chao, modified for contamination"""

"This file takes in a list of FASTA files and outputs a De Bruijn graph in DOT format"
import sys
import getopt
import os
from Bio import AlignIO
import node as n

queue = []  # Initialize a queue

USAGE = '''usage: DeBruijnGraph_generator.py [-h / --help] [-v/ --version] [-o / --ofile FILE] [-k / --kmer k-mer length] [-l / --seqLen sequence length] [-a / --alnNum alignment number] sequence FASTA files'''

# We loop through all the FASTA files and construct the graph
def main(argv):
    k = 3
    verbose = False
    outfile = "tmp.dot"
    seqLen = -1
    alnNum = 3
    try:
        opts, args = getopt.getopt(argv,"hvo:k:l:a:",["ofile=", "kmer=", "seqLen=", "alnNum="])
    except getopt.GetoptError:
        print(USAGE)
        sys.exit(2)
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

    if len(args) < 1:
        print(USAGE)
        print("Please input at least one FASTA file")
        sys.exit(2)

    k = k - 1
    kmer_dic = {}
    organism_id = 0  # Starting organism ID
    for fasta_file in args:
        process_fasta_file(fasta_file, k, seqLen, alnNum, organism_id, kmer_dic)
        organism_id += 1  # Increment for the next file

    print("outfile: ", outfile)
    try:    
        os.remove(outfile)
    except OSError:
        pass
    fw = open(outfile, "a")
    fw.write("strict digraph  {\n")
    visited = []  # List for visited nodes.
    bfs(visited, kmer_dic["$"], fw)
    fw.write("}\n")
    fw.close()

def process_fasta_file(fasta_file, k, seqLen, alnNum, organism_id, kmer_dic):
    alignment = AlignIO.read(fasta_file, "fasta")
    alignment = alignment[:alnNum]
    alignment_len = alignment.get_alignment_length()
    seq_number = len(alignment)
    kmer_set = set()

    for row_idx in range(seq_number):
        row_seq = str(alignment[row_idx].seq)
        row_seq_parsed = row_seq.replace("-", "")
        if seqLen == -1 or seqLen > len(row_seq_parsed):
            lcl_seqLen = len(row_seq_parsed)
        else:
            lcl_seqLen = seqLen
        row_seq_parsed = row_seq_parsed[:lcl_seqLen] + "$"

        for i in range(lcl_seqLen + 1):
            node_kmer = row_seq_parsed[i:i + k]
            kmer_set.add(node_kmer)

    # Constructing and updating nodes with organism ID.
    node_idx = 0
    for node_kmer in kmer_set:
        if node_kmer in kmer_dic:
            kmer_dic[node_kmer].add_organism_id(organism_id)  # If the node already exists, we just add the organism ID (for quick query later)
        else:
            new_node = n.node(node_idx, node_kmer, alignment_len, {organism_id})
            node_idx += 1
            kmer_dic[node_kmer] = new_node

    # Constructing the graph.
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

def bfs(visited, node, fw):
    visited.append(node)
    queue.append(node)
    while queue:
        m = queue.pop(0)
        for m_child in m.children:
            organism_id_label = ','.join(map(str, m_child.organismID))
            edge_label = f"{m_child.nodelabel[0]} ({organism_id_label})"
            fw.write(f"\tS{m.nodeid} -> S{m_child.nodeid} [ label = \"{edge_label}\" ];\n")
        for child in m.children:
            if child not in visited:
                visited.append(child)
                queue.append(child)

if __name__ == "__main__":
    main(sys.argv[1:])