"""
Maybe some comments about implementation here.

References to page numbers in this code are referring to the paper:
[1] Price at al. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix.
    Molecular Biology and Evolution, vol 26 (7). 2009.
    Paper can be found on https://pubmed.ncbi.nlm.nih.gov/19377059/
"""
import math
import sys
from queue import PriorityQueue

import pprint as pp
from typing import Text  # For pretty printing (Replace with own code before submission)
import numpy as np
import argparse

from src.node import Node
from src.tree import Tree
import src.tophits as tophits
# from src.NNI import NNI

def fast_tree(args: argparse.Namespace, sequences: dict) -> str:
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
    # Find a way to track T
    T = update_T(nodes)

    # print("Total profile T")
    # pp.pprint(T)

    #### Heuristics for neighbor joining (with or without top-hits)

    # Top hits sequence
    if not args.no_top_hits:
        nodes = tophits.top_hits_init(nodes, verbose=args.verbose)

    # FastNJ best joins
    fast_NJ_best_joins(nodes, args.no_top_hits, verbose=args.verbose)
  
    #### End of neighbor joining heuristic preparation

    # # # Step 3 : Create initial topology
    CreateInitialTopology(nodes, args.no_top_hits, verbose=args.verbose)
    # print("Initial topology is created:")
    # PrintNewick(nodes)

    # # Step 4 : Nearest Neighbor Interchanges
    # NNI(nodes)
    # print("NNI is finished:")

    # # Final step: print tree topology as Newick string
    # PrintNewick(nodes)

def update_T(nodes: list) -> list:
    """
    We gaan ff breaken
    """
    all_profiles = []
    for node in nodes:
        if node.active:
            all_profiles.append(node.profile)

    print(all_profiles)
    print("")
    print([sum(y) for x in zip(*all_profiles) for y in x])
    print([sum(y) for x in zip(*all_profiles) for y in x] / len(all_profiles))

    return [sum(y) for x in zip(*all_profiles) for y in x] / len(all_profiles)




def fast_NJ_best_joins(nodes: list, no_top_hits: bool, verbose: int=0) -> None:
    """
    The key idea in FastNJ is to store the best join for each node.
    
    The best join for each leaf is determined before the joins begin, and the best join for
    each new interior node is determined when that node is created. When searching for the best join overall,
    FastNJ considers only best join for each node, or n candidates. Thus, FastNJ requires a total of O(N 2)
    time.

    Note: when using top-hits, then at the beginning the best join for each node is found in their top-hits lists.
    """
    if no_top_hits:
        for i in nodes:
            best_join = (sys.maxsize / 2, 0)  # (distance, node index)
            for j in nodes:
                if i.index == j.index:
                    continue
                temp_profile_new_node = averageProfile([i, j])  # calculates profile of potential merge
                criterion = uncorDistance([temp_profile_new_node, i.profile]) - out_distance(i, nodes) - out_distance(j, nodes)

                if criterion < best_join[0]: #if best join for now
                    best_join = (criterion, j.index)

            # Found the best join for i, now update it
            i.best_join = best_join

    else:
        for node in nodes:
            # Nice, the best join for each node is found in their top-hits list already!
            # But wait, while computing the top hits of A, we may discover that A,B is a better join than B,best(B).

            # So if the best_join was already set when A,B is a better join than B,best(B) was true, move on to the next
            if node.best_join is not None:
                continue

            # Okay after that check, we can use tophits
            best_join_dist, best_join = node.tophits.list[0]
            node.best_join = (best_join_dist, best_join)

            best_B_dist, best_B = nodes[best_join].tophits.list[0]

            # A, B is a better join than B, Best(B)
            if best_B_dist > best_join_dist:
                nodes[best_join].best_join = (best_join_dist, node.index)

    if verbose == 1:
        for node in nodes:
            print('FastNJ best join for Node ', node.index, 'is', node.best_join)

def tophits_new_node(new_node: Node) -> PriorityQueue():
    """
    After a join, FastTree computes the top-hits list for the new node in O(mLa) time
    by comparing the node to all entries in the top-hits lists of its children.

    Args:
        new_node (Node) : The newly created inner node after a join
    
    Returns:
        PriorityQueue : the tophits for the new node.
    """

    top_hits = PriorityQueue()

    for th in new_node.leftchild.tophits.queue:
        tophits.put(th)

    for th in new_node.rightchild.tophits.queue:
        tophits.put(th)

    return tophits



def averageProfile(nodes: list) ->  list:
    '''Calculates the average of profiles of internal nodes

        profiles of internalNodes (list): The profiles of internal nodes
    Returns:
        Average profile (list): the profile matrix containing average of two profiles
    '''
    p1 = nodes[0].profile
    p2 = nodes[1].profile
    ptotal = []
    for i in range(len(p1)):
        pbase = []
        for j  in range((4)):
            pbase.append((p1[i][j] + p2[i][j]) / 2)
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


"""Following 25 rules do the same as uncorrectedDistance but uses sequences of nodes as input and returns list of distances

    """


def makeCombisofChildren(children: list) -> list:
    """Make combinations of all children
    
    put combinations of 2 sequences of all nodes in list. 

    Args:
        children (list): sequences of all nodes

    Returns:
        (list): combinations of all sequences
    """
    combi = []
    for i, v1 in enumerate(children):
        for j in range(i + 1, len(children)):
            combi.append([v1, children[j]])
    return combi


def HammingDistance(combi: list) -> list:
    """Calculate hamming distance between 2 nodes
    
    Args:
        combi (list): list with combinations of all nodes 

    Returns:
        hamming distance (list): hamming distances of all input nodes
    """
    distance = []
    for j in range(len(combi)):
        hammingdistance = 0
        for i in range(len(combi[0][0])):
            if combi[j][0][i] != combi[j][1][i]:
                hammingdistance += 1
        distance.append(hammingdistance)
    return distance


def SequenceDistance(combi: list, k: int) -> list:
    """calculate the sequence distance between combinations of all sequences
    
    Args:
        combi (list): list with combinations of all nodes 
        k (int): length of 1 sequence

    Returns:
        SequenceDistance (list): ratio of hamming distance/sequence length for each combination
    """
    ham = HammingDistance(combi)
    seqDis = []
    for i in range(len(ham)):
        seqDis.append(ham[i] / int(k))
    return seqDis


def out_distance(i, nodes):
    """Calculates r distance of i : r(i)

    Args:
        i (Node) : 
    """
    active_nodes = 1  # i is always an active node
    dist_to_others = 0
    dist_to_others1 = 0
    for j in nodes:
        if j.name == i:
            continue
        if not j.active:
            continue

        active_nodes += 1

        profile_i_j = averageProfile([i,j])
        # dist_to_others += uncorrectedDistance(profile_i_j)
        dist_to_others1 += uncorDistance([profile_i_j, i.profile])
    # print(dist_to_others-dist_to_others1)

    # Don't divide by 0
    if active_nodes == 2:
        return dist_to_others

    r = dist_to_others / (active_nodes - 2)
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

    # # Update top-hits for new node
    # new_node.tophits = tophits_new_node(new_node)


    # print("Minimized distance = ", min_dist, "of nodes ", best_join[0].name, best_join[1].name)
    return best_join, new_node


def create_join(nodes, minimized_join, new_node, verbose: int=0):
    """
    Join two nodes.


    """
    # make joined nodes inactive
    nodes[int(minimized_join[1].index)].active = False
    nodes[int(minimized_join[0].index)].active = False

    # append the newly joined node to list of nodes 
    nodes.append(new_node)

    if verbose == 1:
        print(new_node.name, nodes[minimized_join[0].index].branchlength)
        print("Merged nodes to: " + new_node.name)
        print("left child: " + str(nodes[len(nodes)-1].leftchild))
        print("right child: " + str(nodes[len(nodes)-1].rightchild))


    # Recalculate total profile T
    # When we do a join, we also need to update the total profile (the average over all active nodes). To
    # compute this profile takes O(nLa) time. However, we can subtract the joined nodes and add the new
    # node to the total profile in O(La) time. (FastTree) recomputes the total profile from scratch every 200
    # iterations to avoid roundoff errors from accumulating, where the choice of 200 is arbitrary. This adds 
    # another O(N 2La/200) work.)
    update_T(nodes)

    # And update top-hits of new join

    return nodes


def CreateInitialTopology(nodes: list, no_top_hits: bool, verbose: int=0) -> list:

    if no_top_hits:
        numberLeaf = len(nodes)
        for i in range(numberLeaf - 1):
            minimized_join, new_node = minimize_nj_criterion(nodes, len(nodes))
            nodes = create_join(nodes, minimized_join, new_node, verbose)
        
    else:
        # With top hits, we can use We use the FastNJ and local hill-climbing heuristics,
        # with the further restriction that we consider only the top m candidates at each step

        active_nodes = []
        for node in nodes:
            if node.active:
                active_nodes.append(node)

        # First find the best m joins among the best-hit entries for the n active nodes
        # FastTree simply sorts all n entries, but the paper suggests a speed up using a PriorityQueue!
        best_m_joins = PriorityQueue()
        for node in active_nodes:

            # Put in the best join pair (node, best_join)
            best_m_joins.put((node.index, node.best_join[1]))
        
        # Then, we compute the current value of the neighbor-joining criterion for those m candidates,
        # which takes O(mLa) time
        m = active_nodes[0].tophits.m
        best_candidate = (0, 0)  # (distance, node index)]
        min_dist = sys.maxsize / 2
        for _ in range(m):
            candidate = best_m_joins.get()

            i = nodes[candidate[0]]
            j = nodes[candidate[1]]

            temp_profile_new_node = averageProfile([i, j])  # calculates profile of potential merge
            criterion = uncorDistance([temp_profile_new_node, i.profile]) - out_distance(i, nodes) - out_distance(j, nodes)

            if criterion < min_dist: #if best join for now
                best_candidate = candidate
                min_dist = criterion

        print("The best candidate = ", best_candidate)

        # Given this candidate join, we do a local hill-climbing search for a better join
        best_join = local_hill_climbing_tophits(nodes, best_candidate, min_dist, verbose=verbose)

        # Make the join
        create_join(nodes,)



def local_hill_climbing_tophits(nodes: list, best_candidate: tuple, best_dist: float, verbose: int=0) -> tuple:
    """
    Perform Local Hill Climbing with the top-hits heuristic

    Using the top-hits heuristic, we only search within the top-hit lists rather than
    comparing the two nodes to all other active nodes.
    """
    for hit in nodes[best_candidate[0]].tophits.list:
        if hit[0] < best_dist:
            best_candidate = (best_candidate[0], hit[1])
            best_dist = hit[0]

    for hit in nodes[best_candidate[1]].tophits.list:
        if hit[0] < best_dist:
            best_candidate = (hit[1], best_candidate[1])
            best_dist = hit[0]

    return (nodes[best_candidate[0]], nodes[best_candidate[1]])
    
   
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

            # if node was switched, switch parent and children
            if jj != best_top[0][0].index:
                # save indices of node jj
                jj_parent = nodes[jj].parent
                jj_leftchild = nodes[jj].leftchild
                jj_rightchild = nodes[jj].rightchild

                # change indices of node jj to node from better topology
                nodes[jj].parent = best_top[0][0].parent
                nodes[jj].leftchild = best_top[0][0].leftchild
                nodes[jj].rightchild = best_top[0][0].rightchild

                # find the node from the better topology and change indices to the ones from jj
                nodes[best_top[0][0].index].parent = jj_parent
                nodes[best_top[0][0].index].leftchild = jj_leftchild
                nodes[best_top[0][0].index].rightchild = jj_rightchild

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
                kk_leftchild = nodes[kk].leftchild
                kk_rightchild = nodes[kk].rightchild

                # change indices of node kk to node from better topology
                nodes[kk].parent = best_top[0][1].parent
                nodes[kk].leftchild = best_top[0][1].leftchild
                nodes[kk].rightchild = best_top[0][1].rightchild

                # find the node from the better topology and change indices to the ones from kk
                nodes[best_top[0][1].index].parent = kk_parent
                nodes[best_top[0][1].index].leftchild = kk_leftchild
                nodes[best_top[0][1].index].rightchild = kk_rightchild

                # swap indices
                nodes[kk].index = best_top[0][1].index
                nodes[best_top[0][1].index].index = kk

                # swap positions in node list
                node_kk = nodes[kk]
                nodes[kk] = nodes[best_top[0][1].index]
                nodes[best_top[0][1].index] = node_kk

                print("swapped nodes", nodes[jj].index, nodes[kk].index, nodes[rr].index, nodes[ss].index)

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