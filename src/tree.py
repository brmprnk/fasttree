class Tree:

    def __init__(self, nodes, sequences):
        self.nodes = nodes
        self.N = len(nodes)
        self.sequences = sequences

    def newick(self):
        """Returns newick string representation of Tree
        """