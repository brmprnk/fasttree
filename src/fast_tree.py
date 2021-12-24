"""
Maybe some comments about implementation here.

References to page numbers in this code are referring to the paper:
[1] Price at al. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix.
    Molecular Biology and Evolution, vol 26 (7). 2009.
    Paper can be found on https://pubmed.ncbi.nlm.nih.gov/19377059/
"""
import math
import pprint as pp  # For pretty printing (Replace with own code before submission)
import sys
from itertools import combinations

from src.node import Node
# from src.NNI import NNI

def fast_tree(sequences) -> str:
    """FastTree Algorithm.

    Args:
        sequences (dict): Mapping of sequences to their names, as was provided in the program's input

    Returns:
        (str): A phylogenetic tree in Newick format.
    """
    print("The sequences entered into the program : ")
    pp.pprint(sequences)
    nodes = []
    for ii, seq in enumerate(sequences.keys()):
        node = Node(seq, ii, sequences[seq])
        node.leaf = True
        nodes.append(node)

    # Actual first step : Unique sequences ( page 1646, do later )

    # Step 1 of algorithm : Create total profile T
    T = Profile(nodes)  # <----------------------------------------------------------------- why are we doing this??
    # print("Total profile T")
    # pp.pprint(T)

    # Step 2 : Top hits sequence
    # Skip for now

    # Step 3 : Create initial topology
    CreateInitialTopology(nodes, T)
    print("Initial topology is created:")
    # PrintNewick(nodes)

    # Step 4 : Nearest Neighbor Interchanges
    NNI(nodes)
    print("NNI is finished:")

    # Final step: print tree topology as Newick string
    PrintNewick(nodes)

def Profile(nodes: list) -> list:
    """Calculate Profile of internal nodes
    
        internalNodes (list): The sequences of internal nodes
    Args:
    Returns:
        Profile (list): the profile matrix containing ratios 
    """
    sequences = []
    for node in nodes:
        sequences.append(node.sequence)
    columns = [''.join(seq) for seq in zip(*sequences)]
    return [[float(col.count(base)) / float(len(col)) for base in 'ACGT'] for col in columns]

def averageProfile(nodes: list) ->  list:
    '''Calculates the average of profiles of internal nodes

        profiles of internalNodes (list): The profiles of internal nodes
    Returns:
        Average profile (list): the profile matrix containing average of two profiles
    '''
    lambda1 = 0.5
    p1 = nodes[0].profile
    p2 = nodes[1].profile
    ptotal = []
    for i in range(len(p1)):
        pbase = []
        for j  in range((4)):
            pbase.append((p1[i][j] * lambda1 + p2[i][j] * (1 - lambda1)))
        ptotal.append(pbase)
    return ptotal


def profileDistanceNodes(nodes: list, indices: list):
    '''indices [i, j, k]
    ∆(ij, k) = λ∆(i, k) + (1 − λ)∆(j, k)
    profile of parent after joining i and j

    args:
        list of all nodes 
        index of nodes [i, j, k]
    returns:
        Profile distance between joined nodes and other nodes
    '''
    lamb = 0.5
    profDist = lamb * (profileDistanceNew([nodes[indices[0]].profile, nodes[indices[2]].profile ])) + (1-lamb) * (profileDistanceNew([nodes[indices[1]].profile, nodes[indices[2]].profile ]))
    return profDist

def profileDistanceNew(profiles:list) -> float: 
    ''' calculates “profile distance” that is the average distance between profile characters over all positions
    
    input: 
        list of profiles of which you want to calculate profile distance (i,j)
    
    returns:
        profile distance ∆(i, j)
    '''
    value = 0
    for L in range(len(profiles)):
        for a in range(4):

            for b in range(4):
                
                if profiles[0][L][a] > 0 and profiles[1][L][b] > 0:
                    P = 0
                else:
                    P = 1
                value += profiles[0][L][a] * profiles[1][L][b] * P
                
    
    return value / len(profiles)
            

def findAllKids(nodes:list, nodeIndex: int) -> list:
    '''Finds all children of input nodeIndex
        
        Input:
            All nodes (list): list of all nodes
            Index of node (int): index of node of which you want to know the children
        
        Returns:
            All children (list): list with indexes of all children of specific node
    '''
    kids = [] #list of all kids of node with index nodeIndex
    tempKids = [nodeIndex] #queue of list with all kids
    while True:
        acting = tempKids[0]
        tempKids.remove(acting)
        if nodes[acting].rightchild != None:
            right = nodes[acting].rightchild
            tempKids.append(right)
            kids.append(right)
        if nodes[acting].leftchild != None:
            left = nodes[acting].leftchild
            tempKids.append(left)
            kids.append(left)           
        if len(tempKids) == 0: #if no queue left
            # print(kids)
            return kids
        
def uncorDistance(nodes: list, join: list):
    '''uncorrected distance between joined nodes (node_ij) and other nodes (node_k)
    du(ij, k) = ∆(ij, k) − u(ij) − u(k)
    args;
        list with input nodes
        indices of nodes which are joined
    returns:
        uncorrected Distance 
    '''
    if len(join) == 2:
        indices = [join[0].index, join[1].index]
        #du(i, j) = ∆(i, j)−u(i)−u(j), where u(i) = 0 and u(j) = 0 as i and j are leaves
        del_ij = profileDistanceNew([nodes[indices[0]].profile, nodes[indices[1]].profile])
        return del_ij
    if len(join) == 3:
        indices = [join[0], join[1], join[2].index]
        del_ijk = profileDistanceNodes(nodes, indices)
        u_ij = updistance(nodes, [indices[0], indices[1]])
        u_k = updistance(nodes, [indices[2]])
        du_ijk = del_ijk - u_ij - u_k
        return du_ijk   

def updistance(nodes: list, ijk: list):
    ''' calculate updistance with formula's:
            u(ij) ≡ ∆(i, j)/2 
            u(k) has kids so look at leftchild and rightchild so becomes u(k_leftchild, k_rightchild)
            u(k) = o for leaves

    Args:
        list with all nodes
        list with nodes for which the updistance should be calculated (could have a length of 1 or 2 depending on u(ij) or u(k))
    returns:
        updistance u(ij) or u(k) 
    '''
    if len(ijk) > 1:
        return profileDistanceNew([nodes[ijk[0]].profile, nodes[ijk[1]].profile ]) / 2
    elif nodes[ijk[0]].leaf == True: 
        return 0
    else:
        return profileDistanceNew([nodes[nodes[ijk[0]].rightchild].profile, nodes[nodes[ijk[0]].leftchild].profile])  / 2


def outDistanceNew(i: Node, nodes: list, T: list) -> float:
    '''The average profile distance between a node and all other
       nodes can be inferred from the total profile T: r(i) = (n∆(i, T) − ∆(i, i) − (n − 1)u(i) + u(i) − sum u(j))/(n-2)
        
    Args:
        active nodes; list of nodes; T total profile of current topology
    returns:
        out distance of one node 
    '''

    N_active_nodes = 1  # i is always an active node; is 
    sumJ = 0
    for j in nodes:
        if j.name == i:
            continue
        if not j.active:
            continue
        N_active_nodes += 1
        sumJ += updistance(nodes, [j.index])
    # !!!!!!!!!! ∆(i, i) snap ik niet
    sum_du_ij = N_active_nodes * profileDistanceNew([i.profile, T]) - profileDistanceNew([i.profile, i.profile]) - (
        N_active_nodes - 1) * updistance(nodes, [i.index]) + updistance(nodes, [i.index]) - sumJ
    if N_active_nodes == 2:
        return sum_du_ij
    return sum_du_ij / (N_active_nodes - 2)

def minimize_nj_criterion(nodes, index, T):
    """Returns i,j for which d(i, j) - r(i) -r(j) is minimal and corresponding new node of merging i,j
    """
    active_nodes = []
    for node in nodes:
        if node.active:
            active_nodes.append(node)

    min_dist = sys.float_info.max
    best_join = (0, 0)
    for i in active_nodes:
        for j in active_nodes:
            if i == j:
                continue
            temp_profile_new_node = averageProfile([i, j]) #calculates profile of potential merge
            if i.leaf == True and j.leaf == True:
                criterion = uncorDistance(nodes, [i, j]) - outDistanceNew(i, nodes, T) - outDistanceNew(j, nodes, T)
            elif i.leaf == True:
                criterion = uncorDistance(nodes, [j.leftchild, j.rightchild, i]) - outDistanceNew(i, nodes, T) - outDistanceNew(j, nodes, T)
            elif j.leaf == True:
                criterion = uncorDistance(nodes, [i.leftchild, i.rightchild, j]) - outDistanceNew(i, nodes, T) - outDistanceNew(j, nodes, T)
            else:
                criterion = uncorDistance(nodes, [i.leftchild, i.rightchild, j]) - outDistanceNew(i, nodes, T) - outDistanceNew(j, nodes, T)
           
            if criterion < min_dist: #if best join for now
                profile_new_node = temp_profile_new_node #sets profile of new node to profile of best join
                min_dist = criterion
                best_join = (i, j)  # saves best joining nodes
    
    # save just calculated profile of joining nodes to a Node with name containing both joined nodes and make this new node active
    # we should probably change the class Node as the sequence is not known for the merged nodes. I just made a beun oplossing. Don't know if it's good enough
    new_node = Node(str(best_join[0].name) + '&' + str(best_join[1].name), int(index), 'nosequence')
    new_node.profile = profile_new_node
    new_node.active = True
    # add indices of left child, right child
    new_node.leftchild = best_join[0].index
    new_node.rightchild = best_join[1].index
    nodes[best_join[0].index].parent = new_node.index
    nodes[best_join[1].index].parent = new_node.index
    return best_join, new_node


def CreateInitialTopology(nodes, T):
    numberLeaf = len(nodes)
    for i in range(numberLeaf - 1):
        minimized_join, new_node = minimize_nj_criterion(nodes, len(nodes), T)
        # make joined nodes inactive
        nodes[int(minimized_join[1].index)].active = False
        nodes[int(minimized_join[0].index)].active = False
        # append the newly joined node to list of nodes 
        BranchLength(minimized_join, numberLeaf, nodes, new_node)
        nodes.append(new_node)
   
def JC_distance(d_u: float) -> float:
    """Compute Jukes-Cantor distance of FastTree's uncorrected distance

    Defined on page 1643 as d = -(3/4)log(1 - (4/3)d_u).

    Important note: Page 1643-1644
    "For both nucleotide and protein sequences,
     FastTree truncates the corrected distances to a maximum of 3.0 substitutions per site,
     and for sequences that do not overlap because of gaps, FastTree uses this maximum distance."
    
    Args:
        d_u (float): FastTree's uncorrected distance, the fraction of positions that differ between sequences.

    Returns:
        (float): Jukes-Cantor distance.
    """
    # FastTree's max distance
    max_distance = 3.0

    # Distances are truncated to 3 substitutions per site
    # if d_u is 3/4 or larger then the log of a negative value will be calculated, instead return max_dist
    if d_u >= 0.75:
        return max_distance

    # Calculate Jukes-Cantor distance (d in paper)
    # Paper mentions log, literature uses ln.
    # Therefore we also use log base 10
    jd_d = -0.75 * math.log(1 - (4 / 3) * d_u)

    # For sequences that do not overlap, FastTree uses a max distance of 3.0
    if jd_d > max_distance:
        return max_distance

    return jd_d


def BranchLength(minimized_join, numberLeaf, nodes, new_node):
    n1 = minimized_join[0].index
    n2 = minimized_join[1].index
    if n1 < numberLeaf and n2 < numberLeaf:      #connect single leaf with other single leaf
        fraction = uncorDistance(nodes, [nodes[n1], nodes[n2]])
        nodes[n1].branchlength = JC_distance(fraction)
        nodes[n2].branchlength = JC_distance(fraction)
    elif n1 < numberLeaf and n2 >= numberLeaf:    #connect single leaf with other branch
        d12 = uncorDistance(nodes, [nodes[n1], nodes[nodes[n2].leftchild]])
        d13 = uncorDistance(nodes, [nodes[n1], nodes[nodes[n2].rightchild]])
        d23 = uncorDistance(nodes, [nodes[nodes[n2].leftchild], nodes[nodes[n2].rightchild]])
        nodes[n1].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23))/2
        nodes[n2].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23))/2
    elif n2 < numberLeaf and n1 >= numberLeaf:    #connect single leaf with other branch
        d12 = uncorDistance(nodes, [nodes[n2], nodes[nodes[n1].leftchild]])
        d13 = uncorDistance(nodes, [nodes[n2], nodes[nodes[n1].rightchild]])
        d23 = uncorDistance(nodes, [nodes[nodes[n2].leftchild], nodes[nodes[n2].rightchild]])
        nodes[n1].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23))/2
        nodes[n2].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23))/2
    else:                                          #connect two branches
        d13 = uncorDistance(nodes, [nodes[nodes[n1].leftchild],nodes[nodes[n2].leftchild]])
        d14 = uncorDistance(nodes, [nodes[nodes[n1].leftchild],nodes[nodes[n2].rightchild]])
        d23 = uncorDistance(nodes, [nodes[nodes[n1].rightchild],nodes[nodes[n2].leftchild]])
        d24 = uncorDistance(nodes, [nodes[nodes[n1].rightchild],nodes[nodes[n2].rightchild]])
        d12 = uncorDistance(nodes, [nodes[nodes[n1].leftchild],nodes[nodes[n1].rightchild]])
        d34 = uncorDistance(nodes, [nodes[nodes[n2].leftchild],nodes[nodes[n2].rightchild]])
        nodes[n1].branchlength = (JC_distance(d13) + JC_distance(d14) + JC_distance(d23) + JC_distance(d24))/4 - (JC_distance(d12) + JC_distance(d34))/2
        nodes[n2].branchlength = (JC_distance(d13) + JC_distance(d14) + JC_distance(d23) + JC_distance(d24))/4 - (JC_distance(d12) + JC_distance(d34))/2


def PrintNewick(nodes: list):
    """ Generates a Newick string based on the list of nodes.

        The Newick list will initially contain the index of each node.
        Cycle through all nodes to build the structure.
        The index of internal nodes (= has children) will be replaced by the index of their children.
        The index of leaf nodes (= has no children) will be replaced by the name of the leaf node.
        This results in a list with all the names of the nodes with the right hierarchy of brackets.

    Args:
        nodes (list): all nodes in the tree, including internal (merged) nodes

    Returns:
        Newick string (str): the desired format of the tree structure (newick string)
    """

    # Initiate newick string with root node
    for node in nodes:
        if node.parent is None:
            newick_list = ['(', node.leftchild, ',', node.rightchild, ')']
            queue = [node.leftchild, node.rightchild]   # indices of nodes that are waiting to be replaced
    # print('Newick string initialized')

    while len(queue) > 0:
        # print('New while iteration')
        add_to_queue = []
        for ii in queue:
            replace = newick_list.index(ii)             # Find where this node was located in the newick list, this entry should be replaced with the index of the children or name of the node
            if nodes[ii].leaf:                          # skip leaf nodes
                newick_list[replace] = nodes[ii].name   # Replace node index by it's name (this is a leaf node)
            else:   # If not a leaf node,
                newick_list[replace:replace+1] = ('(', nodes[ii].leftchild, ',', nodes[ii].rightchild, ')')   # Replace node index by the index of it's children
                add_to_queue.extend([nodes[ii].leftchild, nodes[ii].rightchild])
            queue.remove(ii)                            # Node is checked, remove from queue
        queue.extend(add_to_queue)  # Add all new nodes in the newick list to the queue
        # print('Newick list at end of iteration', newick_list)

    newick_str = "".join(newick_list)+';'
    print('\nNewick string:', newick_str)


def NNI(nodes):
    print('start NNI')
    nn = sum([node.leaf for node in nodes])  # number of leaf nodes = number of unique sequences

    # print('#nodes:', len(nodes))
    # print('#rounds:', round(math.log(nn) + 1))

    ########
    # Manual swap to see if NNI is swapping
    jj = 9
    kk = 4
    jj_parent = nodes[jj].parent
    # change indices of node jj to node from better topology
    nodes[jj].parent = nodes[kk].parent
    # find the node from the better topology and change indices to the ones from jj
    nodes[kk].parent = jj_parent

    # swap indices
    nodes[jj].index = nodes[kk].index
    nodes[kk].index = jj
    # swap positions in node list
    node_jj = nodes[jj]
    nodes[jj] = nodes[kk]
    nodes[kk] = node_jj

    # Repeat log2(N)+1 times
    for ii in range(round(math.log(nn,2) + 1)):
        # Loop over all nodes
        for node in nodes:
            # Find what other nodes it can be fixed with (and nodes that are attached to them so which ones to compare)
            # go from bottom to top
            if node.leaf:  # skip all leaf nodes as you cannot fix them to try a new topology
                continue
            if node.parent is None:  # skip if ii is the root node
                continue
            if nodes[node.parent].parent is None:  # if jj is root node: choose right child of the root as jj, and both its children as cc and dd
                jj = nodes[node.parent].rightchild
            else:
                # node ii is fixed together with its parent jj
                jj = node.parent

            # get the indices of the nodes that can be swapped
            aa = node.leftchild         # child of original node
            bb = node.rightchild        # child of original node
            cc = nodes[jj].rightchild   # child of second fixed node (jj)
            dd = nodes[jj].parent       # parent of the second fixed node (jj)

            print('NNI compares', nodes[aa].index, nodes[bb].index, nodes[cc].index, nodes[dd].index)

            # For each possible combination of fixed nodes, find the best topology
            best_top = MinimizedEvolution(nodes[aa], nodes[bb], nodes[cc], nodes[dd], nodes)
            print('NNI best topology', best_top[0][0].index, best_top[0][1].index, best_top[1][0].index,
                  best_top[1][1].index)

            # Do all switches

            # maximum of one switch can be made per round, stop checking if switch was found
            # ss is never switched, no need to check
            # if rr is switched, the switch is already taken care of when checking jj or kk, no need to check again

            # if node was switched, switch their parents
            if aa != best_top[0][0].index:
                # save indices of node jj
                aa_parent = nodes[aa].parent

                # change indices of node jj to node from better topology
                nodes[aa].parent = best_top[0][0].parent

                # find the node from the better topology and change indices to the ones from jj
                nodes[best_top[0][0].index].parent = aa_parent

                # swap indices
                nodes[aa].index = best_top[0][0].index
                nodes[best_top[0][0].index].index = aa

                # swap positions in node list
                node_aa = nodes[aa]
                nodes[aa] = nodes[best_top[0][0].index]
                nodes[best_top[0][0].index] = node_aa

                print("swapped nodes", nodes[aa].index, nodes[bb].index, nodes[cc].index, nodes[dd].index)

            elif bb != best_top[0][1].index:
                # save indices of node kk
                bb_parent = nodes[bb].parent

                # change indices of node kk to node from better topology
                nodes[bb].parent = best_top[0][1].parent

                # find the node from the better topology and change indices to the ones from kk
                nodes[best_top[0][1].index].parent = bb_parent

                # swap indices
                nodes[bb].index = best_top[0][1].index
                nodes[best_top[0][1].index].index = bb

                # swap positions in node list
                node_bb = nodes[bb]
                nodes[bb] = nodes[best_top[0][1].index]
                nodes[best_top[0][1].index] = node_bb

                print("swapped nodes", nodes[aa].index, nodes[bb].index, nodes[cc].index, nodes[dd].index)

        # Recompute profiles of internal nodes
        for node in nodes:
            if node.leaf:  # skip all leaf nodes
                continue
            node.profile = averageProfile([nodes[node.leftchild], nodes[node.rightchild]])
    # best_top = MinimizedEvolution(nodes[0], nodes[1], nodes[8], nodes[11])
    # print("Best topology", best_top[0][0].index, best_top[0][1].index, best_top[1][0].index, best_top[1][1].index)

    return nodes


def MinimizedEvolution(n1, n2, n3, n4, nodes):
    """ Evaluate all possible topologies with four nodes surrounding two fixed nodes
    and find the topology that minimizes the evolution citerion.

    Args:
        n1 (node): node

    Returns:
        best_topology (list): list of lists that contain the nodes that form the best topology

    """

    # All possible topologies involving these nodes
    option1 = [[n1, n2], [n3, n4]]
    option2 = [[n1, n3], [n2, n4]]
    option3 = [[n3, n2], [n1, n4]]
    options = [option1, option2, option3]

    # Calculate the evolution criterion for each possible topology
    dist_options = []
    for ii in range(len(options)):
        dist_a = JC_distance(uncorDistance(nodes, [options[ii][0][0], options[ii][0][1]]))
        dist_b = JC_distance(uncorDistance(nodes, [options[ii][1][0], options[ii][1][1]]))
        print(uncorDistance(nodes, [options[ii][0][0], options[ii][0][1]]))
        dist_options.append(dist_a + dist_b)
        # print(uncorDistance([options[ii][0][0].profile, options[ii][0][1].profile]))
        # print(uncorDistance([options[ii][1][0].profile, options[ii][1][1].profile]))

    # Choose the topology with the minimized criterion
    best_topology = options[dist_options.index(min(dist_options))]
    print("dist_option", dist_options)
    return best_topology