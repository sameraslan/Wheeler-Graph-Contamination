import networkx as nx
# import matplotlib.pyplot as plt

DOTFILE = "../Kuanhao_WG_Toolkit/data/graph/DeBruijnGraph/DeBruijn_k_2_a1.dot"
PATTERN = "TCA"
CURRENT_PATTERN_IDX = 0
LAST_EXPLORED_NODE = None
ORGANISM_IDS = []

# Read the graph
G = nx.nx_pydot.read_dot(DOTFILE)

for node in G.nodes():
    # TODO: this doesn't explore new branch properly 
    # e.g. (S1, S2), (S2, S3), (S1, S4) doesn't explore S1 -> S4
    # We should cache something in the branching node
    CURRENT_PATTERN_IDX = 0
    ORGANISM_IDS = []

    # Traverse the graph using DFS
    for edge in nx.dfs_edges(G, node):
        # If it is a new branch, reset the pattern matching
        if LAST_EXPLORED_NODE and LAST_EXPLORED_NODE != edge[0]:
            # Explored new branch, reseting
            # TODO: cache the CURRENT_PATTERN_IDX here?
            break

        # Match the pattern
        if G.get_edge_data(edge[0], edge[1])[0]["label"] == PATTERN[CURRENT_PATTERN_IDX]:
            # Match

            # Take note of the organism id
            # ORGANISM_IDS.append(G.get_edge_data(edge[0], edge[1])[0]["organismID"])
            ORGANISM_IDS.append("a1")

            # If the pattern is fully matched
            if CURRENT_PATTERN_IDX == len(PATTERN) - 1:
                print("Pattern fully matched, starts at {} ends at {}".format(node, edge[1]))
                print("Organism IDs: {}".format(ORGANISM_IDS))
                # TODO: take a majority vote of the organism ids
                break

            # Increment the pattern index
            CURRENT_PATTERN_IDX += 1
        else:
            # Doesn't match, break
            break

        # Take note of the last node to determine if we're exploring new branch
        LAST_EXPLORED_NODE = edge[1]


