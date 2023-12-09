import networkx as nx

# Function to extract the nucleotide character from a graph edge label
def extract_nucleotide(label):
    label = label.strip('\"')
    return label.split()[0]  # ex. returns 'A'

# Function to extract organism IDs from a graph edge label
def extract_organism_ids(label):
    label = label.strip('\"')
    ids_part = label.split('(')[-1].strip(')')
    return set([int(id) for id in ids_part.split(',')])  # ex. returns {0, 1, 3}

# Function to query the De Bruijn graph for a given pattern
def query_graph(G, pattern):
    paths = []

    # Look for starting nodes that correspond to the first char of the pattern
    for node in G.nodes:
        for edge in G.edges(node, data=True):
            if extract_nucleotide(edge[2]['label']) == pattern[0]:
                paths.append(([node, edge[1]], extract_organism_ids(edge[2]['label'])))

    # For each path, extend it by adding edges that match the next char of the pattern
    for char in pattern[1:]:
        new_paths = []
        for path, ids in paths:
            last_node = path[-1]
            for edge in G.edges(last_node, data=True):
                if extract_nucleotide(edge[2]['label']) == char:
                    target_node = edge[1]
                    new_path = path + [target_node]
                    new_ids = ids & extract_organism_ids(edge[2]['label'])  # Intersecting organism IDs
                    if new_ids:  # Only consider the path if there are common organism IDs
                        new_paths.append((new_path, new_ids))

        paths = new_paths
        if not paths:
            return [], set()

    # If at least one valid path is found, return the organism IDs from any of the paths
    if paths:
        #print(paths)
        return paths[0][1]  # Only need to return the IDs from one path as all paths will have the same IDs

    return set()

# Ensure valid fatsq file
def is_valid_input_file(input):
    # Every valid FASTQ record should have exactly 4 lines
    if len(input) % 4 != 0:
        return False, "Error: The FASTQ file does not have a multiple of 4 lines."

    for i in range(0, len(input), 4):
        # The first line of each record should start with the "@" symbol
        if not input[i].startswith("@"):
            return False, f"Error: Line {i + 1} does not start with '@'."

        # The third line of each record should start with the "+" symbol
        if not input[i + 2].startswith("+"):
            return False, f"Error: Line {i + 3} does not start with '+'."

        # The sequence and quality strings should have the same length
        if len(input[i + 1].strip()) != len(input[i + 3].strip()):
            return False, f"Error: Length mismatch between sequence and quality on lines {i + 2} and {i + 4}."

    return True, ""

# Function to read input file
def read_input_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Make sure we have a valid input file
    is_valid, error_msg = is_valid_input_file(lines)
    if not is_valid:
        print(error_msg)
        sys.exit(1)

    reads = [line.strip() for line in lines[1::4]]  # Starting on line 1 every four lines
    quality_chars = lines[3::4]  # Starting on line 3 every fout lines
    quality_scores = [[ord(char) - 33 for char in line.strip()] for line in quality_chars]

    return reads, quality_scores

# Below main is to test on individual patterns
if __name__ == "__main__":
    DOTFILE = "./graphs/example.dot"
    G = nx.nx_pydot.read_dot(DOTFILE)
    pattern = "CTAGG"
    organism_ids = query_graph(G, pattern)
    print(f"Organism(s) with this Pattern: {organism_ids}")