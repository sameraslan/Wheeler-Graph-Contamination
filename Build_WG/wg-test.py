import sys
import random
import argparse

def create_graph_from_seq(sequence, filename):
    num_nodes = len(sequence) + 1
    node_names = [f'S{i}' for i in range(1, num_nodes+1)]
    
    edges = []
    for i, c in enumerate(sequence):
        edges.append((0, i+1, c))
        
    with open(filename, 'w') as f:
        f.write('strict digraph {\n')
        for src, dst, label in edges:
            src_name = node_names[src] 
            dst_name = node_names[dst]
            f.write(f'\t{src_name} -> {dst_name} [ label = "{label}" ];\n')
        f.write('}')

def main():
    parser = argparse.ArgumentParser(description='Generate Wheeler Graph from Sequence')
    parser.add_argument('-s', '--sequence', type=str, required=True, help='Input sequence')
    parser.add_argument('-o', '--outfile', type=str, default='out.dot', help='Output filename (default: out.dot)')

    args = parser.parse_args()
    seq = args.sequence
    outfile = args.outfile
    
    create_graph_from_seq(seq, outfile)


if __name__ == '__main__':
    main()