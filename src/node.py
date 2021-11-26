class Node:

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.branchlength = 0
        self.leftchild = None
        self.rightchild = None
        self.active = True

    def set_inactive(self):
        self.active = False