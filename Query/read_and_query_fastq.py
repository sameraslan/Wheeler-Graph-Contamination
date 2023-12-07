import networkx as nx
# import matplotlib.pyplot as plt

DOTFILE = "../Kuanhao_WG_Toolkit/data/graph/DeBruijnGraph/DeBruijn_k_2_a1.dot"
PATTERN = "GATCCA" # S1 -> S3 -> S6 -> S4 -> S8 -> S7

"""
Given a Wheeler Graph G and a pattern P, return the organism IDs that match the pattern.
Assumption 1: no loop
Limitation 1: returns only 1 pattern location
Limitation 2: exact matching

Args:
    G: a wheeler graph with labels
    pattern: string of pattern to be exactly matched

Returns:
    List of path graphs

"""
def query_graph(G, pattern):
    dfs_edges = list(nx.edge_dfs(G))

    starting_edge = 0 # index of dfs_edges to start to, 
    # Either
    # 1. the first tuple in dfs_edges (init)
    # 2. earliest tuple whose 'label' matches the first char of the pattern
    # 3. Next tuple after the last explored node
    next_starting_edge = 0 # index of dfs_edges to start to,
    current_pattern_idx = 0

    # Variable for branching logic
    last_explored_node = None
    pattern_idx_cache = {}
    # organism_ids_cache = {}
    path_graph_cache = {}

    # Variable for results
    # organism_ids = []
    current_path_graph = nx.DiGraph()
    path_graphs = []

    # for edge in dfs_edges[possible_starting_edge:]:
    while starting_edge + current_pattern_idx <= len(dfs_edges):
        edge = dfs_edges[starting_edge + current_pattern_idx]

        # If new branch, reset from the cache
        if last_explored_node and last_explored_node != edge[0]:
            current_pattern_idx = pattern_idx_cache[edge[0]]
            current_path_graph = path_graph_cache[edge[0]]
            # This is a hack, might collide with next_starting_edge?
            starting_edge = starting_edge - current_pattern_idx + 1
            next_starting_edge = starting_edge

        # If edge[0] node branches, take note of the current_pattern_idx
        if G.out_degree(edge[0]) > 1:
            pattern_idx_cache[edge[0]] = current_pattern_idx
            path_graph_cache[edge[0]] = current_path_graph

        # If a non-starting edge equals the first char of the pattern, make a note of it
        if G.get_edge_data(edge[0], edge[1])[0]["label"] == pattern[0]:
            if next_starting_edge == starting_edge: # Ignore the rest if found
                next_starting_edge = starting_edge + current_pattern_idx
        
        # ACTUAL MATCHING LOGIC STARTS
        # If pattern matches, append the organism id and increment the pattern index
        if G.get_edge_data(edge[0], edge[1])[0]["label"] == pattern[current_pattern_idx]:
            nx.add_path(current_path_graph, [edge[0], edge[1]])
            current_pattern_idx += 1

            # If the pattern is fully matched, return the organism ids
            if current_pattern_idx == len(pattern):
                return [current_path_graph]
                # path_graphs.append(current_path_graph)
            

        # else, if pattern doesn't match...
        else:
            # If no next_starting_edge, just continue from next
            if next_starting_edge == starting_edge:
                starting_edge = starting_edge + current_pattern_idx + 1
                next_starting_edge = starting_edge
            # Else, if next_starting_edge is found, match pattern from there
            else:
                starting_edge = next_starting_edge
            
            current_pattern_idx = 0
            current_path_graph = nx.DiGraph()
        # ACTUAL MATCHING LOGIC ENDS

        # Take note of the last explored node to detect branch
        last_explored_node = edge[1]

    # Not found
    print("Not found")
    return []

if __name__ == "__main__":
    G = nx.nx_pydot.read_dot(DOTFILE)
    # Manually remove loop
    G.remove_edges_from([('S7', 'S6'), ('S2', 'S4')])
    # print(G.edges.data())
    a = query_graph(G, PATTERN)
    print("Found {} match".format(len(a)))
    if(len(a) > 0):
        print(a[0].edges)
        # [('S0', 'S1'), ('S1', 'S3'), ('S3', 'S6'), ('S6', 'S4'), ('S4', 'S8'), ('S8', 'S7')]
