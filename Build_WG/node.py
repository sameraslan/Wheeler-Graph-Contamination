"""
Node class for representing nodes in a De Bruijn graph.
Adapted from `https://github.com/Kuanhao-Chao/Wheeler_Graph_Toolkit`.
Modified by Samer Aslan, Owen Bianchi, Anirudh Kashyap, and Senapati Diwangkara.
"""

class node:
    # Constructor to initialize a node
    def __init__(self, nodeID, nodeLabel, nodeColid, organismID=None):
        # Node attributes: ID, label, column ID, set of organism IDs, parents and children
        self.nodeid = nodeID  # Unique identifier for the node
        self.nodelabel = nodeLabel  # Label of the node (usually a k-mer)
        self.nodecolid = nodeColid  # Column ID, used for aligning nodes in visualization
        self.organismID = organismID if organismID else set()  # Set of organism IDs that this node is associated with
        self.parents = []  # List of parent nodes
        self.children = []  # List of child nodes
        print("Initializing a node: ", self.nodeid, self.nodelabel)

    # Method to add a parent node if it's not already a parent
    def add_parent(self, parent):
        if parent not in self.parents:
            self.parents.append(parent)
            
    # Method to add a child node if it's not already a child
    def add_child(self, child):
        if child not in self.children:
            self.children.append(child)

    # Method to add an organism ID to this node
    def add_organism_id(self, organismID):
        self.organismID.add(organismID)  # Adding a new organism ID to the set

    # Setter method to update the node ID
    def set_nodeid(self, nodeID):
        self.nodeid = nodeID  # Updating the node's ID
