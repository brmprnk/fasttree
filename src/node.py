class Node:

    def __init__(self, name, index, sequence):
        self.name = name
        self.sequence = sequence
        self.branchlength = 0
        self.leftchild = None
        self.rightchild = None
        self.active = True
        self.index = index
        if sequence != 'nosequence':
            self.profile = self.create_profile(sequence)

    def create_profile(self, sequence):
        profile = []
        nt_values = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        for nt in sequence:
            column = [0 for i in range(4)]
            column[nt_values[nt]] = 1.0
            profile.append(column)
        return profile

    def set_inactive(self):
        self.active = False