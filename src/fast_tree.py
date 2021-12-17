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
    for ii,seq in enumerate(sequences.keys()):
        node = Node(seq, ii, sequences[seq])
        node.leaf = True
        nodes.append(node)

    # Actual first step : Unique sequences ( page 1646, do later )

    # Step 1 of algorithm : Create total profile T
    # T = Profile(nodes)  # <----------------------------------------------------------------- why are we doing this??
    # print("Total profile T")
    # pp.pprint(T)

    # Step 2 : Top hits sequence
    # Skip for now

    # Step 3 : Create initial topology
    CreateInitialTopology(nodes)
    print("Initial topology is created:")
    PrintNewick(nodes)

    # Step 4 : Nearest Neighbor Interchanges
    NNI(nodes)
    print("NNI is finished:")

    # Final step: print tree topology as Newick string
    PrintNewick(nodes)


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


def uncorDistance(profiles: list) -> float:
    differ = 0
    length = len(profiles[1])
    for k in range(length):

        if profiles[0][k][:] != profiles[1][k][:]:
            differ += 1


    fraction = differ / length
    return fraction


def out_distance(i, nodes):
    """Calculates r distance of i : r(i)

    Args:
        i (Node) : 
    """
    N_active_nodes = 1  # i is always an active node; is 
    dist_to_others = 0
    dist_to_others1 = 0
    for j in nodes:
        if j.name == i:
            continue
        if not j.active:
            continue

        N_active_nodes += 1

        profile_i_j = averageProfile([i,j])
        # dist_to_others += uncorrectedDistance(profile_i_j)
        dist_to_others1 += uncorDistance([profile_i_j, i.profile])
    # print(dist_to_others-dist_to_others1)

    # Don't divide by 0
    if N_active_nodes == 2:
        return dist_to_others

    r = dist_to_others / (N_active_nodes - 2)
    # print("Out distance r({}) = ".format(i.name), r)
    return r


def minimize_nj_criterion(nodes, index):
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
            criterion = uncorDistance([temp_profile_new_node, i.profile]) - out_distance(i, nodes) - out_distance(j, nodes)

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

    # print("Minimized distance = ", min_dist, "of nodes ", best_join[0].name, best_join[1].name)
    return best_join, new_node


def CreateInitialTopology(nodes):
    numberLeaf = len(nodes)
    for i in range(numberLeaf - 1):
        minimized_join, new_node = minimize_nj_criterion(nodes, len(nodes))
        # make joined nodes inactive
        nodes[int(minimized_join[1].index)].active = False
        nodes[int(minimized_join[0].index)].active = False
        # append the newly joined node to list of nodes 
        BranchLength(minimized_join, numberLeaf, nodes, new_node)
        nodes.append(new_node)
        
        # print(new_node.name, nodes[minimized_join[0].index].branchlength)
        # print("Merged nodes to: " + new_node.name)
        # print("left child: " + str(nodes[len(nodes)-1].leftchild))
        # print("right child: " + str(nodes[len(nodes)-1].rightchild))
   
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
        fraction = uncorDistance([nodes[n1].profile, nodes[n2].profile])
        nodes[n1].branchlength = JC_distance(fraction)
        nodes[n2].branchlength = JC_distance(fraction)
    elif n1 < numberLeaf and n2 >= numberLeaf:    #connect single leaf with other branch
        d12 = uncorDistance([nodes[n1].profile, nodes[nodes[n2].leftchild].profile])
        d13 = uncorDistance([nodes[n1].profile, nodes[nodes[n2].rightchild].profile])
        d23 = uncorDistance([nodes[nodes[n2].leftchild].profile, nodes[nodes[n2].rightchild].profile])
        nodes[n1].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23))/2
        nodes[n2].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23))/2
    elif n2 < numberLeaf and n1 >= numberLeaf:    #connect single leaf with other branch
        d12 = uncorDistance([nodes[n2].profile, nodes[nodes[n1].leftchild].profile])
        d13 = uncorDistance([nodes[n2].profile, nodes[nodes[n1].rightchild].profile])
        d23 = uncorDistance([nodes[nodes[n2].leftchild].profile, nodes[nodes[n2].rightchild].profile])
        nodes[n1].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23))/2
        nodes[n2].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23))/2
    else:                                          #connect two branches
        d13 = uncorDistance([nodes[nodes[n1].leftchild].profile,nodes[nodes[n2].leftchild].profile])
        d14 = uncorDistance([nodes[nodes[n1].leftchild].profile,nodes[nodes[n2].rightchild].profile])
        d23 = uncorDistance([nodes[nodes[n1].rightchild].profile,nodes[nodes[n2].leftchild].profile])
        d24 = uncorDistance([nodes[nodes[n1].rightchild].profile,nodes[nodes[n2].rightchild].profile])
        d12 = uncorDistance([nodes[nodes[n1].leftchild].profile,nodes[nodes[n1].rightchild].profile])
        d34 = uncorDistance([nodes[nodes[n2].leftchild].profile,nodes[nodes[n2].rightchild].profile])
        nodes[n1].branchlength = (JC_distance(d13) + JC_distance(d14) + JC_distance(d23) + JC_distance(d24))/4 - (JC_distance(d12) + JC_distance(d34))/2
        nodes[n2].branchlength = (JC_distance(d13) + JC_distance(d14) + JC_distance(d23) + JC_distance(d24))/4 - (JC_distance(d12) + JC_distance(d34))/2


def PrintNewick(nodes: list):
    """ Generates a Newick string based on the list of nodes.
        Uses the left and right child of each node as structure.
        Assumes that the internal (merged) nodes are added to the nodes list in order of merging and that the children saved as indices.

        The Newick list will initially contain the index of each node.
        Cycle through all nodes from last added to the first to build the structure.
        The index of internal nodes (= has children) will be replaced by the index of their children.
        The index of leaf nodes (= has no children) will be replaces by the name of the leaf node.
        This results in a list with all the names of the nodes with the right hierarchy of brackets.

    Args:
        nodes (list): all nodes in the tree, including internal (merged) nodes

    Returns:
        Newick string (str): the desired format of the tree structure
    """

    # Initiate newick list for the root node (last added node to the list)
    newick_list = ['(', nodes[-1].leftchild, ',', nodes[-1].rightchild, ')']
    # print('Initial newick list', newick_list)
    # Update newick list for each node
    for node in reversed(nodes[:-1]):
        # print('search for', node.index)
        replace = newick_list.index(node.index) # Find where this node was located in the newick list, this entry should be replaced with the index of the children or name of the node

        # Check if node has children
        if node.leftchild is None:              # Right node not checked as it always has zero or two children
            newick_list[replace] = node.name     # Replace node index by it's name (this is a leaf node)
        # + ':' + str(round(node.branchlength,3))
        else:
            newick_list[replace:replace+1] = ('(', node.leftchild, ',', node.rightchild,')')   # Replace node index by the index of it's children
        # print('newick list at end of iteration', newick_list)  
    # print('Newick list', newick_list)

    newick_str = "".join(newick_list)+';'
    print('\nNewick string:', newick_str)


def NNI(nodes):
    print('start NNI')
    nn = sum([node.leaf for node in nodes])  # number of leaf nodes = number of unique sequences

    print('#nodes:', len(nodes))
    print('#rounds:', round(math.log(nn) + 1))

    # Manual swap to see if NNI is swapping
    jj = 8
    jj_parent = nodes[jj].parent
    print(jj_parent)
    # change indices of node jj to node from better topology
    nodes[jj].parent = nodes[9].parent
    print(nodes[9].parent)
    # find the node from the better topology and change indices to the ones from jj
    nodes[nodes[9].index].parent = jj_parent

    # swap indices
    nodes[jj].index = nodes[9].index
    nodes[nodes[9].index].index = jj

    # swap positions in node list
    node_jj = nodes[jj]
    nodes[jj] = nodes[nodes[9].index]
    nodes[nodes[9].index] = node_jj
    print(nodes[9].parent)
    print(nodes[8].parent)




    # Repeat log2(N)+1 times
    for ii in range(round(math.log(nn) + 1)):
        # Loop over all nodes
        for node in nodes:
            # Find what other nodes it can be fixed with (and nodes that are attached to them so which ones to compare)
            # go from bottom to top
            if node.leaf:  # skip all leaf nodes as you cannot fix them to try a new topology
                continue
            if node.parent is None:  # skip the root node
                continue
            if nodes[node.parent].parent is None:  # skip the nodes with root node as parent
                continue

            # node ii is fixed together with its parent qq
            qq = node.parent

            # get the indices of the nodes that can be swapped
            jj = node.leftchild     # child of original node
            kk = node.rightchild    # child of original node
            ss = nodes[qq].parent   # parent of the second fixed node
            # find the other child of node qq (original node is either left or right child)
            if nodes[qq].rightchild == node.index:
                rr = nodes[qq].leftchild    # child of second fixed node (qq)
            else:
                rr = nodes[qq].rightchild   # child of second fixed node (qq)
            print('NNI compares', jj, kk, ss, rr)

            # For each possible combination of fixed nodes, find the best topology
            best_top = MinimizedEvolution(nodes[jj], nodes[kk], nodes[ss], nodes[rr])
            print('NNI best topology', best_top[0][0].index, best_top[0][1].index, best_top[1][0].index,
                  best_top[1][1].index)

            # Do all switches

            # maximum of one switch can be made per round, stop checking if switch was found
            # ss is never switched, no need to check
            # if rr is switched, the switch is already taken care of when checking jj or kk, no need to check again

            # if node was switched, switch their parents
            if jj != best_top[0][0].index:
                # save indices of node jj
                jj_parent = nodes[jj].parent

                # change indices of node jj to node from better topology
                nodes[jj].parent = best_top[0][0].parent

                # find the node from the better topology and change indices to the ones from jj
                nodes[best_top[0][0].index].parent = jj_parent

                # swap indices
                nodes[jj].index = best_top[0][0].index
                nodes[best_top[0][0].index].index = jj

                # swap positions in node list
                node_jj = nodes[jj]
                nodes[jj] = nodes[best_top[0][0].index]
                nodes[best_top[0][0].index] = node_jj

                print("swapped nodes", nodes[jj].index, nodes[kk].index, nodes[rr].index, nodes[ss].index)

            elif kk != best_top[0][1].index:
                # save indices of node kk
                kk_parent = nodes[kk].parent
                # kk_leftchild = nodes[kk].leftchild
                # kk_rightchild = nodes[kk].rightchild

                # change indices of node kk to node from better topology
                nodes[kk].parent = best_top[0][1].parent
                # nodes[kk].leftchild = best_top[0][1].leftchild
                # nodes[kk].rightchild = best_top[0][1].rightchild

                # find the node from the better topology and change indices to the ones from kk
                nodes[best_top[0][1].index].parent = kk_parent
                # nodes[best_top[0][1].index].leftchild = kk_leftchild
                # nodes[best_top[0][1].index].rightchild = kk_rightchild

                # swap indices
                nodes[kk].index = best_top[0][1].index
                nodes[best_top[0][1].index].index = kk

                # swap positions in node list
                node_kk = nodes[kk]
                nodes[kk] = nodes[best_top[0][1].index]
                nodes[best_top[0][1].index] = node_kk

                print("swapped nodes", nodes[jj].index, nodes[kk].index, nodes[rr].index, nodes[ss].index)

        # Recompute profiles of internal nodes
        for node in nodes:
            if node.leaf:  # skip all leaf nodes
                continue
            node.profile = averageProfile([nodes[node.leftchild], nodes[node.rightchild]])
    # best_top = MinimizedEvolution(nodes[0], nodes[1], nodes[8], nodes[11])
    # print("Best topology", best_top[0][0].index, best_top[0][1].index, best_top[1][0].index, best_top[1][1].index)

    return nodes


def MinimizedEvolution(n1, n2, n3, n4):
    """ Evaluate all possible topologies with four nodes surrounding two fixed nodes
    and find the topology that minimizes the evolution citerion.

    Args:
        n1 (node): node

    Returns:
        best_topology (list): list of lists that contain the nodes that form the best topology

    """

    # All possible topologies involving these nodes
    option1 = [[n1, n2], [n3, n4]]
    option2 = [[n1, n4], [n3, n2]]
    option3 = [[n4, n2], [n3, n1]]
    options = [option1, option2, option3]

    # Calculate the evolution criterion for each possible topology
    dist_options = []
    for ii in range(len(options)):
        dist_a = JC_distance(uncorDistance([options[ii][0][0].profile, options[ii][0][1].profile]))
        dist_b = JC_distance(uncorDistance([options[ii][1][0].profile, options[ii][1][1].profile]))
        dist_options.append(dist_a + dist_b)
        # print(uncorDistance([options[ii][0][0].profile, options[ii][0][1].profile]))
        # print(uncorDistance([options[ii][1][0].profile, options[ii][1][1].profile]))

    # Choose the topology with the minimized criterion
    best_topology = options[dist_options.index(min(dist_options))]
    print("dist_option", dist_options)
    return best_topology