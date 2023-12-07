import networkx as nx
# import matplotlib.pyplot as plt

DOTFILE = "./wheel.dot"
DOTFILE = "./graphs/test-nonorg.dot"
FASTQ_FILE = "./hw2q1_A.fastq"

def query_graph(G, pattern):
    """
    Given a Wheeler Graph G and a pattern P, return the organism IDs that match the pattern.
    Doesn't really use the consecutive property of the pattern, but it's a start.
    
    Args:
        G: a wheeler graph with labels
        pattern: string of pattern to be exactly matched

    Returns:
        List of path graphs

    """

    matching_edges = []
    paths = []

    def matching_func(edge, char, prev_matching_edges):
        if edge[2]['label'] == char:
            if len(prev_matching_edges) == 0:
                return True
            else:
                return edge[0] in [e[1] for e in prev_matching_edges]
        return False
    

    for pattern_idx in range(len(pattern) - 1, 0, -1):
        new_matching_edges = list(filter(lambda e: matching_func(e, pattern[pattern_idx], matching_edges), G.edges.data()))
        
        # If no matching edges, return empty list
        if len(new_matching_edges) == 0:
            return []
    
        
        if len(paths) == 0:
            # If path is empty, append all matching edges
            paths = [(e[0], e[1]) for e in new_matching_edges]
        else:
            # Append paths
            for edge in new_matching_edges:
                for path in paths:
                    if path[-1][1] == edge[0]:
                        path.append(edge)
                        

    return paths

if __name__ == "__main__":
    G = nx.nx_pydot.read_dot(DOTFILE)
    pattern = "CTAG"
    paths = query_graph(G, pattern)
    print(paths)

# def parse_fastq(fh):
#     """ Parse reads from a FASTQ filehandle.  For each read, we
#         return a name, nucleotide-string, quality-string triple. """
#     reads = []
#     while True:
#         first_line = fh.readline()
#         if len(first_line) == 0:
#             break  # end of file
#         name = first_line[1:].rstrip()
#         seq = fh.readline().rstrip()
#         fh.readline()  # ignore line starting with +
#         qual = fh.readline().rstrip()
#         reads.append((name, seq, qual))
#     return reads

# if __name__ == "__main__":
#     G = nx.nx_pydot.read_dot(DOTFILE)

#     with open(FASTQ_FILE) as fh:
#         reads = parse_fastq(fh)

#     for read in reads:
#         print(read[1])
#         a = query_graph(G, read[1])
#         print(a)
