import networkx as nx
import sys
# import matplotlib.pyplot as plt

FASTQ_FILE = "./hw2q1_A.fastq"
DOTFILE = "./graphs/test.dot"

def extract_nucleotide(label):
    label = label.strip('\"')
    return label.split()[0] # ex. returns 'A'

def extract_organism_ids(label):
    label = label.strip('\"')
    ids_part = label.split('(')[-1].strip(')')
    return [int(id) for id in ids_part.split(',')]  # ex. returns [0, 1, 3]

def query_graph(G, pattern):
    paths = []
    organism_ids = []

    # Look for starting nodes that correspond to the first char of the pattern
    for node in G.nodes:
        for edge in G.edges(node, data=True):
            if extract_nucleotide(edge[2]['label']) == pattern[0]: # Check if edge matches first character
                # Start a new path with this edge
                paths.append([edge[1]])  # Return the target node of the starting edges

    # For each path (valid starting node), extend it by adding edges that match the next char of the pattern
    for char in pattern[1:]:
        new_paths = []
        print("\nchar:", char)
        print("paths:", paths)
        for path in paths:
            last_node = path[-1]  # Get the last node of the path
            print("last node:", last_node)
            
            for edge in G.edges(last_node, data=True):
                print("edge:", edge)
                if extract_nucleotide(edge[2]['label']) == char:
                    target_node = edge[1]
                    print("target:", target_node)
                    new_path = [last_node].append(target_node)  # extends current path by adding target node of the edge
                    print("new path:", new_path)
                    new_paths.append(new_path) 
                    print("new paths:", new_paths)
                    organism_ids.extend(extract_organism_ids(edge[2]['label']))

        paths = new_paths

        if not paths:
            return []

    return set(organism_ids)

if __name__ == "__main__":
    G = nx.nx_pydot.read_dot(DOTFILE)
    pattern = "CTAG"
    organism_ids = query_graph(G, pattern)
    print(organism_ids)


# # Ensure valid fatsq file
# def is_valid_input_file(input):
#     # Every valid FASTQ record should have exactly 4 lines
#     if len(input) % 4 != 0:
#         return False, "Error: The FASTQ file does not have a multiple of 4 lines."

#     for i in range(0, len(input), 4):
#         # The first line of each record should start with the "@" symbol
#         if not input[i].startswith("@"):
#             return False, f"Error: Line {i + 1} does not start with '@'."
        
#         # The third line of each record should start with the "+" symbol
#         if not input[i + 2].startswith("+"):
#             return False, f"Error: Line {i + 3} does not start with '+'."
        
#         # The sequence and quality strings should have the same length
#         if len(input[i + 1].strip()) != len(input[i + 3].strip()):
#             return False, f"Error: Length mismatch between sequence and quality on lines {i + 2} and {i + 4}."
    
#     return True, ""

# # Function to read input file
# def read_input_file(filename):
#     with open(filename, 'r') as file:
#         lines = file.readlines()

#     # Make sure we have a valid input file
#     is_valid, error_msg = is_valid_input_file(lines)
#     if not is_valid:
#         print(error_msg)
#         sys.exit(1)

#     reads = [line.strip() for line in lines[1::4]]  # Starting on line 1 every four lines
#     quality_chars = lines[3::4]  # Starting on line 3 every fout lines
#     quality_scores = [[ord(char) - 33 for char in line.strip()] for line in quality_chars]

#     return reads, quality_scores

# if __name__ == "__main__":
#     G = nx.nx_pydot.read_dot(DOTFILE)

#     reads, quality_scores = read_input_file(FASTQ_FILE)

#     for read in reads:
#         a = query_graph(G, read)
#         print(a)
