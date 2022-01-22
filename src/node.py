class Node:

    def __init__(self, name, index, sequence, identical_sequences):
        self.name = name
        self.sequence = sequence
        self.branchlength = 0
        self.leftchild = None
        self.rightchild = None
        self.parent = None
        self.active = True
        self.leaf = False
        self.index = index
        self.identical_sequences = identical_sequences  # Default is [self.name]

        # Top-hits is uninitialized
        self.tophits = None

        # FastNJ best join heuristic
        self.best_join = None

        if sequence != 'nosequence':
            self.profile = self.create_profile(sequence)

    def create_profile(self, sequence: str) -> list: 
        """ Create a profile of an input sequence

        Args:
            sequence (str) : sequence of leaf node

        Returns:
            profile (list) : profile matrix of leaf node
        """
        profile = []
        nt_values = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        for nt in sequence:
            column = [0 for i in range(4)]
            column[nt_values[nt]] = 1.0
            profile.append(column)
        return profile

    def set_inactive(self): # inactivates nodes 
        self.active = False

    def __str__(self):
        return "Node {} with children ({}, {}) and parent {}".format(self.index, self.leftchild, self.rightchild, self.parent)
