class node:
    def __init__(self, nodeID, nodeLabel, nodeColid, organismID=None):
        self.nodeid = nodeID
        self.nodelabel = nodeLabel
        self.nodecolid = nodeColid
        self.organismID = organismID if organismID else set()  # set of organism IDs for this node
        self.parents = []
        self.children = []
        print("Initializing a node: ", self.nodeid, self.nodelabel)

    def add_parent(self, parent):
        if parent not in self.parents:
            self.parents.append(parent)
        # else:
        #     print("The parent was added before!")

    def add_child(self, child):
        if child not in self.children:
            self.children.append(child)
        # else:
        #     print("The child was added before!")

    # Add organism label to the node
    def add_organism_id(self, organismID):
        self.organismID.add(organismID)

    def set_nodeid(self, nodeID):
        self.nodeid = nodeID

    # def merge_node(self, node):
    #     # check children
    #     for c in self.children:
    #         print(c)
    #     # check parents
    #     for c in self.children:
    #         print(c)
        