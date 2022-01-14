import math

import src.util as util


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
            These can be joined together to form the Newick string.

            Duplicate sequences are taken into consideration by allowing a node to have more than two children,
            as described in Price et. al. [1]: "It constructs a tree for the unique subset of sequences and then creates
            multifurcating nodes, without support values, as parents of the redundant sequences."

        Returns:
            Newick string (str): the desired format of the tree structure (newick string)
        """
        # Initiate newick string with root node
        for node in self.nodes:
            if node.parent is None:
                newick_list = ['(', node.leftchild, ',', node.rightchild, ')']
                queue = [node.leftchild, node.rightchild]  # indices of nodes that are waiting to be replaced
        # print('Newick string initialized')

        while len(queue) > 0:
            add_to_queue = []
            for ii in queue:
                replace = newick_list.index(ii)  # Find where this node was located in the newick list, this entry should be replaced with the index of the children or name of the node
                if self.nodes[ii].leaf:  # skip leaf nodes
                    # if only one copy of this sequence was in the alignment
                    if len(self.nodes[ii].identical_sequences) == 1:
                        newick_list[replace] = str(self.nodes[ii].name) + ":" + str(round(self.nodes[ii].branchlength,3))  # Replace node index by it's name (this is a leaf node), and corresponding branch length
                    # if multiple copies of this sequence were in the alignment file
                    else:
                        copy_names = self.nodes[ii].identical_sequences     # names of all the identical sequences
                        node_duplicates = str(copy_names[0]) + ":" + str(round(self.nodes[ii].branchlength, 3)) # initialize by first name (no comma in front)
                        for name in copy_names[1:]:                         # add all other names (preceding by comma)
                            node_duplicates = node_duplicates + "," + str(name) + ":" + str(round(self.nodes[ii].branchlength, 3))
                        newick_list[replace] = node_duplicates              # replace node index by all names
                else:  # If not a leaf node
                    newick_list[replace:replace + 1] = ('(', self.nodes[ii].leftchild, ',', self.nodes[ii].rightchild,')')  # Replace node index by the index of it's children
                    add_to_queue.extend([self.nodes[ii].leftchild, self.nodes[ii].rightchild])  # add the children to the queue, they still need to be checked
                queue.remove(ii)  # Node is checked, remove from queue
            queue.extend(add_to_queue)  # Add all new nodes in the newick list to the queue

        # Join all names from the newick list to form one string
        newick_str = "".join(newick_list) + ';'
        print('Newick string:', newick_str)

        return newick_str

    def BranchLength(self):
        """Compute Branch lengths for each node

        Args:
            ft (Tree) : The Tree object
            minimized_join (list): containing the just joined nodes
        """
        nr_leafs = len(self.nodes)

        for n1 in self.nodes:
            for n2 in self.nodes:
                if n1 == n2: #if nodes are the same
                    continue
                elif n1.leaf == False and n2.leaf == False: #if internal node
                    continue
                else:
                    # n1 = minimized_join[0].index
                    # n2 = minimized_join[1].index
                    # connect single leaf with other single leaf
                    if n1.leaf == True and n2.leaf == True:
                        fraction = util.uncorrected_distance(self, [self.nodes[n1.index], self.nodes[n2.index]])
                        self.nodes[n1.index].branchlength = util.JC_distance(fraction)
                        self.nodes[n2.index].branchlength = util.JC_distance(fraction)
                    # connect single leaf with other branch
                    elif n1.leaf == True:
                        d12 = util.uncorrected_distance(self, [self.nodes[n1.index], self.nodes[n2.leftchild]])
                        d13 = util.uncorrected_distance(self, [self.nodes[n1.index], self.nodes[n2.rightchild]])
                        d23 = util.uncorrected_distance(self, [self.nodes[n2.leftchild], self.nodes[n2.rightchild]])
                        self.nodes[n1.index].branchlength = (util.JC_distance(d12) + util.JC_distance(d13) - util.JC_distance(d23)) / 2
                        self.nodes[n2.index].branchlength = (util.JC_distance(d12) + util.JC_distance(d13) - util.JC_distance(d23)) / 2
                    # connect single leaf with other branch
                    elif n2.leaf == True:
                        d12 = util.uncorrected_distance(self, [self.nodes[n2.index], self.nodes[n1.leftchild]])
                        d13 = util.uncorrected_distance(self, [self.nodes[n2.index], self.nodes[n1.rightchild]])
                        d23 = util.uncorrected_distance(self, [self.nodes[n1.leftchild], self.nodes[n1.rightchild]])
                        self.nodes[n1.index].branchlength = (util.JC_distance(d12) + util.JC_distance(d13) - util.JC_distance(d23)) / 2
                        self.nodes[n2.index].branchlength = (util.JC_distance(d12) + util.JC_distance(d13) - util.JC_distance(d23)) / 2
                    # connect two branches
                    # else:
                    #     d13 = util.uncorrected_distance(self, [self.nodes[self.nodes[n1].leftchild], self.nodes[self.nodes[n2].leftchild]])
                    #     d14 = util.uncorrected_distance(self, [self.nodes[self.nodes[n1].leftchild], self.nodes[self.nodes[n2].rightchild]])
                    #     d23 = util.uncorrected_distance(self, [self.nodes[self.nodes[n1].rightchild], self.nodes[self.nodes[n2].leftchild]])
                    #     d24 = util.uncorrected_distance(self, [self.nodes[self.nodes[n1].rightchild], self.nodes[self.nodes[n2].rightchild]])
                    #     d12 = util.uncorrected_distance(self, [self.nodes[self.nodes[n1].leftchild], self.nodes[self.nodes[n1].rightchild]])
                    #     d34 = util.uncorrected_distance(self, [self.nodes[self.nodes[n2].leftchild], self.nodes[self.nodes[n2].rightchild]])
                    #     self.nodes[n1].branchlength = (util.JC_distance(d13) + util.JC_distance(d14) + util.JC_distance(d23) + util.JC_distance(
                    #         d24)) / 4 - (
                    #                                         util.JC_distance(d12) + util.JC_distance(d34)) / 2
                    #     self.nodes[n2].branchlength = (util.JC_distance(d13) + util.JC_distance(d14) + util.JC_distance(d23) + util.JC_distance(
                    #         d24)) / 4 - (
                    #                                         util.JC_distance(d12) + util.JC_distance(d34)) / 2