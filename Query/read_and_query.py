import networkx as nx

DOTFILE = "./graphs/test.dot"

def extract_nucleotide(label):
    label = label.strip('\"')
    return label.split()[0]  # ex. returns 'A'

def extract_organism_ids(label):
    label = label.strip('\"')
    ids_part = label.split('(')[-1].strip(')')
    return set([int(id) for id in ids_part.split(',')])  # ex. returns {0, 1, 3}

def query_graph(G, pattern):
    paths = []
    organism_ids_sets = []

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
            return []

    print(paths)
    # Collect all unique organism IDs from valid paths
    valid_organism_ids = set()
    for _, ids in paths:
        valid_organism_ids.update(ids)

    return valid_organism_ids

if __name__ == "__main__":
    G = nx.nx_pydot.read_dot(DOTFILE)
    pattern = "CTAG"
    organism_ids = query_graph(G, pattern)
    print(organism_ids)