import math


class Tree:
    """

    """

    def __init__(self, nodes, args):
        self.nodes = nodes  # List object holding all nodes in the Tree

        # Constants as described in the paper
        self.N = len(nodes)
        self.m = int(math.sqrt(self.N))

        self.lambda1 = 0.5  # lambda1 (float): lambda value BIONJ

        # Step 1 of algorithm : Create total profile T
        self.T = self.update_T()

        # And set variables dictated by user input flags
        self.no_top_hits = args.no_top_hits
        self.verbose = args.verbose

    def update_T(self) -> list:
        """Updating total profile T with current list of nodes.

        The total profile T is the average profile over all active nodes.
        When we do a join, we also need to update the total profile.

        Returns:
            T (list): Total profile T
        """
        all_profiles = []
        for node in self.nodes:
            if node.active:
                all_profiles.append(node.profile)

        T = []  # Total profile
        for n in range(len(all_profiles[0])):  # loop over length of string
            Tx = [0, 0, 0, 0]  # frequency of total frequency vector is 0 for each base
            for x in range(len(all_profiles)):  # loop over #active nodes in list
                for y in range(4):
                    Tx[y] += all_profiles[x][n][y]  # add values of frequency vector to total frequency vector (Tx)
                    if x == len(all_profiles) - 1:  # if frequency vector of last node is added divide by #active nodes
                        Tx[y] = Tx[y] / len(all_profiles)
            T.append(Tx)  # append total frequency vectors to total profile

        self.T = T
        return T

    def newick_str(self):
        """ Generates a Newick string based on the list of nodes.

            The Newick list will initially contain the index of each node.
            Cycle through all nodes to build the structure.
            The index of internal nodes (= has children) will be replaced by the index of their children.
            The index of leaf nodes (= has no children) will be replaced by the name of the leaf node.
            This results in a list with all the names of the nodes with the right hierarchy of brackets.

        Returns:
            Newick string (str): the desired format of the tree structure (newick string)
        """

        print("Newick_str, ", len(self.nodes))
        # Initiate newick string with root node
        for node in self.nodes:
            print("Node", node.index)
            if node.parent is None:
                newick_list = ['(', node.leftchild, ',', node.rightchild, ')']
                queue = [node.leftchild, node.rightchild]  # indices of nodes that are waiting to be replaced
        # print('Newick string initialized')

        print(len(queue))
        while len(queue) > 0:
            # print('New iteration while loop')
            add_to_queue = []
            for ii in queue:
                replace = newick_list.index(
                    ii)  # Find where this node was located in the newick list, this entry should be replaced with the index of the children or name of the node
                if self.nodes[ii].leaf:  # skip leaf nodes
                    newick_list[replace] = self.nodes[ii].name  # Replace node index by it's name (this is a leaf node)
                else:  # If not a leaf node,
                    newick_list[replace:replace + 1] = ('(', self.nodes[ii].leftchild, ',', self.nodes[ii].rightchild,
                                                        ')')  # Replace node index by the index of it's children
                    add_to_queue.extend([self.nodes[ii].leftchild, self.nodes[ii].rightchild])
                queue.remove(ii)  # Node is checked, remove from queue
            queue.extend(add_to_queue)  # Add all new nodes in the newick list to the queue
            # print('Newick list at end of iteration', newick_list)

        newick_str = "".join(newick_list) + ';'
        print('\nNewick string:', newick_str)

        return newick_str