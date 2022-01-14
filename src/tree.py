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
                replace = newick_list.index(
                    ii)  # Find where this node was located in the newick list, this entry should be replaced with the index of the children or name of the node
                if self.nodes[ii].leaf:  # skip leaf nodes
                    newick_list[replace] = str(self.nodes[ii].name) + ":" + str(round(self.nodes[ii].branchlength,3))  # Replace node index by it's name (this is a leaf node), and corresponding branch length
                else:  # If not a leaf node,
                    newick_list[replace:replace + 1] = ('(', self.nodes[ii].leftchild, ',', self.nodes[ii].rightchild,
                                                        ')')  # Replace node index by the index of it's children
                    add_to_queue.extend([self.nodes[ii].leftchild, self.nodes[ii].rightchild])
                queue.remove(ii)  # Node is checked, remove from queue
            queue.extend(add_to_queue)  # Add all new nodes in the newick list to the queue
            # print('Newick list at end of iteration', newick_list)

        newick_str = "".join(newick_list) + ';'
        print('Newick string:', newick_str)

        return newick_str

    def BranchLength(self):
        """Compute Branch lengths for each leaf node

        Args:
            self : The Tree object
            
        """
        for n1 in self.nodes: 
            if n1.leaf == True:        
                d = util.uncorrected_distance(self, [self.nodes[n1.index], self.nodes[n1.parent]]) # uncorrected distance between node and parent
                self.nodes[n1.index].branchlength = util.JC_distance(d) # append branchlength to node

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

        if len(join) == 2:
            # indices = [join[0].index, join[1].index]
            V_ij = util.profile_distance([join[0].profile, join[1].profile])
            return V_ij

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
        # print('sumV', sumV)
        # print('deelstreep', (2 * (N_active - 2) * V_ij))
        # print('V_ij', V_ij)
        new_lambda = abs(0.5 + (sumV) / (2 * (N_active - 2) * V_ij))
        print('new', abs(new_lambda))
        if new_lambda > 0.999:
            new_lambda = 0.5
        #update lambda
        self.lambda1 = new_lambda


    
