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
                    newick_list[replace:replace + 1] = ('(', self.nodes[ii].leftchild, ',', self.nodes[ii].rightchild,')', ":" , str(round(self.nodes[ii].branchlength, 3)) ) # Replace node index by the index of it's children
                    add_to_queue.extend([self.nodes[ii].leftchild, self.nodes[ii].rightchild])  # add the children to the queue, they still need to be checked
                queue.remove(ii)  # Node is checked, remove from queue
            queue.extend(add_to_queue)  # Add all new nodes in the newick list to the queue

        # Join all names from the newick list to form one string
        newick_str = "".join(newick_list) + ';'
        print('Newick string:', newick_str)

        return newick_str

    def BranchLength(self):
        """Compute Branch lengths for each leaf node

        Args:
            self : The Tree object
            
        """
        for n1 in self.nodes: 
            for n2 in self.nodes:
                if n1 == n2:
                    continue
                if n1.parent == n2.parent:
                    # connect single leaf with other single leaf
                    if n1.leaf == True  and n2.leaf == True:
                        d1 = util.uncorrected_distance(self, [self.nodes[n1.index], self.nodes[n1.parent]]) # uncorrected distance between leaf and parent
                        self.nodes[n1.index].branchlength = util.JC_distance(d1) # append branchlength to node
                        d2 = util.uncorrected_distance(self, [self.nodes[n2.index], self.nodes[n1.parent]]) # uncorrected distance between second leaf and parent
                        self.nodes[n2.index].branchlength = util.JC_distance(d2) # append branchlength to node
                    # connect single leaf with internal node
                    elif n1.leaf == True and n2.leaf == False:   
                        # distance of internal node d(1,23) = (d12+d13-d23)/2 
                        d12 = util.uncorrected_distance(self, [self.nodes[n1.index], self.nodes[n2.leftchild]])
                        d13 = util.uncorrected_distance(self, [self.nodes[n1.index], self.nodes[n2.rightchild]])
                        d23 = util.uncorrected_distance(self, [self.nodes[n2.leftchild], self.nodes[n2.rightchild]])
                        self.nodes[n2.index].branchlength = (util.JC_distance(d12) + util.JC_distance(d13) - util.JC_distance(d23)) / 2
                        # distance of leaf
                        d1 = util.uncorrected_distance(self, [self.nodes[n1.index], self.nodes[n1.parent]]) # uncorrected distance between leaf and parent
                        self.nodes[n1.index].branchlength = util.JC_distance(d1) # append branchlength to node
                        
                    # connect single leaf with internal node
                    elif n2.leaf == True and n1.leaf == False:  
                        # distance of internal node d(1,23) = (d12+d13-d23)/2 
                        d12 = util.uncorrected_distance(self, [self.nodes[n2.index], self.nodes[n1.leftchild]])
                        d13 = util.uncorrected_distance(self, [self.nodes[n2.index], self.nodes[n1.rightchild]])
                        d23 = util.uncorrected_distance(self, [self.nodes[n1.leftchild], self.nodes[n1.rightchild]])
                        self.nodes[n1.index].branchlength = (util.JC_distance(d12) + util.JC_distance(d13) - util.JC_distance(d23)) / 2
                        # distance of leaf
                        d2 = util.uncorrected_distance(self, [self.nodes[n2.index], self.nodes[n1.parent]]) # uncorrected distance between leaf and parent
                        self.nodes[n2.index].branchlength = util.JC_distance(d2) # append branchlength to node

                    # connect two internal nodes    
                    elif n1.leaf == False and n2.leaf == False:
                        # distance of two internal nodes d(12,34) = (d13+d14+d23+d24)/4 - (d12+d34)/2 
                        d13 = util.uncorrected_distance(self, [self.nodes[n1.leftchild], self.nodes[n2.leftchild]])
                        d14 = util.uncorrected_distance(self, [self.nodes[n1.leftchild], self.nodes[n2.rightchild]])
                        d23 = util.uncorrected_distance(self, [self.nodes[n1.rightchild], self.nodes[n2.leftchild]])
                        d24 = util.uncorrected_distance(self, [self.nodes[n1.rightchild], self.nodes[n2.rightchild]])
                        d12 = util.uncorrected_distance(self, [self.nodes[n1.leftchild], self.nodes[n1.rightchild]])
                        d34 = util.uncorrected_distance(self, [self.nodes[n2.leftchild], self.nodes[n2.rightchild]])
                        self.nodes[n1.index].branchlength = (util.JC_distance(d13) + util.JC_distance(d14) + util.JC_distance(d23) + util.JC_distance(d24)) / 4 - (
                                util.JC_distance(d12) + util.JC_distance(d34)) / 2

    def variance_correcion(self, ij):
        ''' calculate variance correction with formula's:
                v(ij) ≡ ∆(i, j)/2
                v(k) has kids so look at leftchild and rightchild so becomes v(k_leftchild, k_rightchild)
                v(k) = o for leaves

        Args:
            tree object
            list with nodes for which the variance correction should be calculated (could have a length of 1 or 2 depending on v(ij) or v(k))
        returns:
            variance correction v(ij) or v(k)
        '''
        if len(ij) > 1:
            return util.profile_distance([ij[0].profile, ij[1].profile])
        elif ij[0].leaf:
            return 0
        else:
            return util.profile_distance(
                [self.nodes[ij[0].rightchild].profile, self.nodes[ij[0].leftchild].profile])


    def variance(self, join: list):
        """Variance between nodes is given by V (i, j) = ∆(i, j) − ν(i) − ν(j)

        args:
            tree object
            list containing nodes which are joined

        returns:
            Variance V_ij
        """
        # if leafs
        if len(join) == 2:
            V_ij = util.profile_distance([join[0].profile, join[1].profile])
            return V_ij
        # if internal node
        if len(join) == 3:
            indice = [join[0].index, join[1].index, join[2].index]
            del_ij = util.profile_distance_nodes(self, indice)
            u_i = self.variance_correcion([join[0], join[1]])
            u_j = self.variance_correcion([join[2]])
            V_ij = del_ij - u_i - u_j
            return V_ij      


    def update_lambda(self, join: list):
        ''' calculate a new lambda value to minimize the variance of the distance estimates for the new node ij,
            using the formula λ =1/2 + SUM(V (j, k) − V (i, k))/(2(n − 2)V (i, j))

        Args:
            tree object
            list with nodes which are just joined

        '''
        #Check if joined nodes have children
        i = join[0]
        j = join[1]
        if i.leaf and j.leaf:
            join = [i, j]
        elif i.leaf:
            join = [self.nodes[j.leftchild], self.nodes[j.rightchild], i]
        elif j.leaf:
            join = [self.nodes[i.leftchild], self.nodes[i.rightchild], j]
        else:
            join = [self.nodes[i.leftchild], self.nodes[i.rightchild], j]

        #calculate variance of nodes including internal nodes (children)
        V_ij = self.variance(join)

        #count number of active nodes
        N_active = 1
        for j in self.nodes:
            if j.name == join[0] or j.name == join[1]:
                continue
            if not j.active:
                continue
            N_active += 1

        # Given these variances, BIONJ weights the join of i, j so as to minimize the variance of the distance estimates for the new node ij, using the formula λ =1/2 + SUM(V (j, k) − V (i, k))/(2(n − 2)V (i, j))
        # FT computes the numerator with (n − 2)(ν(i) − ν(j)) + (j, k) − (i, k) see outdistance for calculation of sums using T
        sumV = (N_active - 2) * (self.variance_correcion([join[1]]) - self.variance_correcion([join[0]])) + N_active * util.profile_distance([join[0].profile, self.T]) - N_active * util.profile_distance([join[1].profile, self.T])
        new_lambda = 0.5 + (sumV) / (2 * (N_active - 2) * V_ij)
        
        # check if lambda is within boundaries [0,1]
        if new_lambda > 1:
            new_lambda = 1
        if new_lambda < 0:
            new_lambda = 0

        #update lambda
        self.lambda1 = new_lambda


    
