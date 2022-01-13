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
from typing import Tuple

import pprint as pp
import argparse
from operator import itemgetter

from src.node import Node
from src.tree import Tree
import src.heuristics as heuristics
import src.neighbor_joining as neighbor_joining
import src.util as util


def fast_tree(args: argparse.Namespace, sequences: dict) -> str:
    """FastTree Algorithm.

    This main function will follow the steps as discussed in the paper Figure 1:

    1) Create alignment Total Profile T ( Done inside Tree Object )
    2) Initialize Top-Hits (and FastNJ heuristic)
    3) Initial Topology
    4) NNI, containing Log2(N) + 1 rounds of NNIs
    5) Computing the final branch lengths
    6) Returning the Newick String representation of the phylogenetic tree.

    Args:
        args: Namespace of the argparse library, containing all user inputs
        sequences (dict): Mapping of sequences to their names, as was provided in the program's input

    Returns:
        (str): A phylogenetic tree in Newick format.
    """
    # Create list of Nodes representing the sequences
    nodes = []
    for ii, seq in enumerate(sequences.keys()):
        node = Node(seq, ii, sequences[seq])
        node.leaf = True
        nodes.append(node)

    # Create Tree object
    ft = Tree(nodes, args)

    if ft.verbose == 1:
        print("The sequences entered into the program : ")
        for key, value in sequences.items():
            print(key, ':', value)

    # Actual first step : Unique sequences ( page 1646, do later )

    # Heuristics for neighbor joining (with or without top-hits)

    # Top hits sequence
    if not ft.no_top_hits:
        heuristics.TopHits.top_hits_init(ft)

    # Initialize FastNJ Best-hit heuristic
    heuristics.fastNJ_init(ft)

    # End of neighbor joining heuristic preparation

    # Create initial topology
    CreateInitialTopology(ft)

    if ft.verbose == 1:
        print("Initial topology is created containing", len(ft.nodes), "nodes : ")
        ft.newick_str()
        print()

    # Nearest Neighbor Interchanges
    NNI(ft)
    if ft.verbose == 1:
        print("NNI is finished")

    # Final step: print tree topology as Newick string
    return ft.newick_str()

def average_profile(nodes: list, lambda1: float) -> list:
    """Calculates the average of profiles of internal nodes

    Args:
        nodes: List of two Nodes whose average profile is calculated
        lambda1 (float): lambda value BIONJ
    Returns:
        Average profile (list): the profile matrix containing average of two profiles
    """
    p1 = nodes[0].profile
    p2 = nodes[1].profile

    ptotal = []
    for i in range(len(p1)):
        pbase = []
        for j in range(4):
            pbase.append((p1[i][j] * lambda1 + p2[i][j] * (1 - lambda1)))
        ptotal.append(pbase)

    return ptotal


def findAllKids(nodes: list, nodeIndex: int) -> list:
    """Finds all children of input nodeIndex
        
        Input:
            All nodes (list): list of all nodes
            Index of node (int): index of node of which you want to know the children
        
        Returns:
            All children (list): list with indexes of all children of specific node
    """
    kids = []  # list of all kids of node with index nodeIndex
    tempKids = [nodeIndex]  # queue of list with all kids
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
        if len(tempKids) == 0:  # if no queue left
            return kids


def create_join(ft: Tree, best_join) -> None:
    """Create a new Node and join the best two nodes under it.
    Join two nodes.

    Args:
        ft (Tree): A Tree object
        best_join (tuple(Node, Node)): The two Nodes to be joined

    Returns:
        None
    """
    # Save just calculated profile of joining nodes to a Node with name containing both joined nodes and make this new node active
    # we should probably change the class Node as the sequence is not known for the merged nodes. I just made a beun oplossing. Don't know if it's good enough
    new_node = Node(str(best_join[0].name) + '&' + str(best_join[1].name), len(ft.nodes), 'nosequence')
    new_node.profile = average_profile([best_join[0], best_join[1]], ft.lambda1)

    # add indices of left child, right child
    new_node.leftchild = best_join[0].index
    new_node.rightchild = best_join[1].index
    ft.nodes[best_join[0].index].parent = new_node.index
    ft.nodes[best_join[1].index].parent = new_node.index

    # make joined nodes inactive
    ft.nodes[int(best_join[0].index)].active = False
    ft.nodes[int(best_join[1].index)].active = False

    BranchLength(ft, best_join)

    # append the newly joined node to list of nodes
    ft.nodes.append(new_node)

    if ft.verbose == 1:
        print(new_node.name, ft.nodes[best_join[0].index].branchlength)
        print("Merged nodes to: " + new_node.name)
        print("left child: " + str(new_node.leftchild))
        print("right child: " + str(new_node.rightchild))

    # Recalculate total profile T
    # When we do a join, we also need to update the total profile (the average over all active nodes). To
    # compute this profile takes O(nLa) time. However, we can subtract the joined nodes and add the new
    # node to the total profile in O(La) time. (FastTree) recomputes the total profile from scratch every 200
    # iterations to avoid roundoff errors from accumulating, where the choice of 200 is arbitrary. This adds 
    # another O(N 2La/200) work.)
    ft.update_T()

    # Update Top-Hits heuristic if applicable
    if not ft.no_top_hits:
        new_node.tophits = heuristics.TopHits(ft.m)
        new_node.tophits.tophits_new_node(ft, new_node)

        # No more top-hits means all nodes have been joined!
        if len(new_node.tophits.list) == 0:
            if ft.verbose == 1:
                print("Newly created node ", new_node.index, " has no top-hits. This means this was the last join!")
                print()
            return

    # Update the best-hit according to FastNJ
    heuristics.fastNJ_update(ft, new_node)


def CreateInitialTopology(ft: Tree) -> None:
    """
    Create the initial topology given a list with all input nodes. 

    Args:
        ft (Tree): Tree object
    """

    # Original number of nodes (nr of sequences)
    nr_leafs = len(ft.nodes)
    for i in range(nr_leafs - 1):
        if ft.verbose == 1:

            active_nodes = []
            for node in ft.nodes:
                if node.active:
                    active_nodes.append(node.index)

            print("Active nodes remaining in initial topology creation : ", active_nodes)

        if ft.no_top_hits:
            # TODO: Implement FastNJ and local hill-climbing for non-top-hits usage
            minimized_join = neighbor_joining.minimize_nj_criterion(ft)
            create_join(ft, minimized_join)

        else:
            # With top hits, we can use We use the FastNJ and local hill-climbing heuristics,
            # with the further restriction that we consider only the top m candidates at each step

            # First find the best m joins among the best-hit entries for the n active nodes
            # FastTree simply sorts all n entries, but the paper suggests a speed-up using a PriorityQueue!
            best_m_joins = PriorityQueue()
            for node in ft.nodes:
                if not node.active:
                    continue

                # Fix FastNJ references to inactive notes the "lazy" way
                if not ft.nodes[node.best_join[1]].active:
                    heuristics.fastNJ_update(ft, node)

                # Put in the best join pair (node, best_join)
                best_m_joins.put((node.best_join[0], (node.index, node.best_join[1])))

            # Then, we compute the current value of the neighbor-joining criterion for those m candidates,
            # which takes O(mLa) time
            best_candidate = (0, 0)  # (distance, node index)]
            min_dist = sys.maxsize / 2
            for _ in range(ft.m):
                candidate = best_m_joins.get()[1]

                i = ft.nodes[candidate[0]]
                j = ft.nodes[candidate[1]]

                criterion = neighbor_joining.nj_criterion(ft, i, j)

                if criterion < min_dist:  # if best join for now
                    best_candidate = candidate
                    min_dist = criterion

            # Given this candidate join, we do a local hill-climbing search for a better join
            best_join = heuristics.local_hill_climb(ft, best_candidate, min_dist)

            # Make the join
            create_join(ft, best_join)





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


def BranchLength(ft: Tree, minimized_join: list):
    """Compute Branch lengths for each node 
    
    Args:
        minimized_join (list): containing the just joined nodes
        numberLeaf (int): total number of leafs (e.g. number of input nodes)
        nodes (list): all nodes including joined ones
        lambda1 (float): lambda value BIONJ 
    """
    nr_leafs = len(ft.nodes)

    n1 = minimized_join[0].index
    n2 = minimized_join[1].index
    # connect single leaf with other single leaf
    if n1 < nr_leafs and n2 < nr_leafs:
        fraction = util.uncorrected_distance(ft, [ft.nodes[n1], ft.nodes[n2]])
        ft.nodes[n1].branchlength = JC_distance(fraction)
        ft.nodes[n2].branchlength = JC_distance(fraction)
    # connect single leaf with other branch
    elif n1 < nr_leafs <= n2:
        d12 = util.uncorrected_distance(ft, [ft.nodes[n1], ft.nodes[ft.nodes[n2].leftchild]])
        d13 = util.uncorrected_distance(ft, [ft.nodes[n1], ft.nodes[ft.nodes[n2].rightchild]])
        d23 = util.uncorrected_distance(ft, [ft.nodes[ft.nodes[n2].leftchild], ft.nodes[ft.nodes[n2].rightchild]])
        ft.nodes[n1].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23)) / 2
        ft.nodes[n2].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23)) / 2
    # connect single leaf with other branch
    elif n2 < nr_leafs <= n1:
        d12 = util.uncorrected_distance(ft, [ft.nodes[n2], ft.nodes[ft.nodes[n1].leftchild]])
        d13 = util.uncorrected_distance(ft, [ft.nodes[n2], ft.nodes[ft.nodes[n1].rightchild]])
        d23 = util.uncorrected_distance(ft, [ft.nodes[ft.nodes[n2].leftchild], ft.nodes[ft.nodes[n2].rightchild]])
        ft.nodes[n1].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23)) / 2
        ft.nodes[n2].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23)) / 2
    # connect two branches
    else:
        d13 = util.uncorrected_distance(ft, [ft.nodes[ft.nodes[n1].leftchild], ft.nodes[ft.nodes[n2].leftchild]])
        d14 = util.uncorrected_distance(ft, [ft.nodes[ft.nodes[n1].leftchild], ft.nodes[ft.nodes[n2].rightchild]])
        d23 = util.uncorrected_distance(ft, [ft.nodes[ft.nodes[n1].rightchild], ft.nodes[ft.nodes[n2].leftchild]])
        d24 = util.uncorrected_distance(ft, [ft.nodes[ft.nodes[n1].rightchild], ft.nodes[ft.nodes[n2].rightchild]])
        d12 = util.uncorrected_distance(ft, [ft.nodes[ft.nodes[n1].leftchild], ft.nodes[ft.nodes[n1].rightchild]])
        d34 = util.uncorrected_distance(ft, [ft.nodes[ft.nodes[n2].leftchild], ft.nodes[ft.nodes[n2].rightchild]])
        ft.nodes[n1].branchlength = (JC_distance(d13) + JC_distance(d14) + JC_distance(d23) + JC_distance(d24)) / 4 - (
                JC_distance(d12) + JC_distance(d34)) / 2
        ft.nodes[n2].branchlength = (JC_distance(d13) + JC_distance(d14) + JC_distance(d23) + JC_distance(d24)) / 4 - (
                JC_distance(d12) + JC_distance(d34)) / 2


def NNI(ft: Tree):
    nn = sum([node.leaf for node in ft.nodes])  # number of leaf nodes = number of unique sequences

    # print('#nodes:', len(nodes))
    # print('#rounds:', round(math.log(nn) + 1))

    ########
    # Manual swap to see if NNI is swapping
    jj = 9
    kk = 4
    jj_parent = ft.nodes[jj].parent
    # change indices of node jj to node from better topology
    ft.nodes[jj].parent = ft.nodes[kk].parent
    # find the node from the better topology and change indices to the ones from jj
    ft.nodes[kk].parent = jj_parent

    # swap indices
    ft.nodes[jj].index = ft.nodes[kk].index
    ft.nodes[kk].index = jj
    # swap positions in node list
    node_jj = ft.nodes[jj]
    ft.nodes[jj] = ft.nodes[kk]
    ft.nodes[kk] = node_jj

    # Repeat log2(N)+1 times
    for ii in range(round(math.log(nn, 2) + 1)):
        # Loop over all nodes
        for node in ft.nodes:
            # Find what other nodes it can be fixed with (and nodes that are attached to them so which ones to compare)
            # go from bottom to top
            if node.leaf:  # skip all leaf nodes as you cannot fix them to try a new topology
                continue
            if node.parent is None:  # skip if ii is the root node
                continue
            if ft.nodes[
                node.parent].parent is None:  # if jj is root node: choose right child of the root as jj, and both its children as cc and dd
                jj = ft.nodes[node.parent].rightchild
            else:
                # node ii is fixed together with its parent jj
                jj = node.parent

            # get the indices of the nodes that can be swapped
            aa = node.leftchild  # child of original node
            bb = node.rightchild  # child of original node
            cc = ft.nodes[jj].rightchild  # child of second fixed node (jj)
            dd = ft.nodes[jj].parent  # parent of the second fixed node (jj)
            # print(aa, bb, cc, dd)
            # print(jj)
            # print(node.index)
            # print(nodes[jj].leftchild)
            # error because index of node is equal to index of parent -> it seems like something goes wrong with indices after swapping
            try:
                if ft.verbose == 1:
                    print('NNI compares', ft.nodes[aa].index, ft.nodes[bb].index, ft.nodes[cc].index, ft.nodes[dd].index)
                    print()

                # For each possible combination of fixed nodes, find the best topology
                best_top = MinimizedEvolution(ft, ft.nodes[aa], ft.nodes[bb], ft.nodes[cc], ft.nodes[dd])
            except:
                break

            if ft.verbose == 1:
                print('NNI best topology', best_top[0][0].index, best_top[0][1].index, best_top[1][0].index,
                      best_top[1][1].index)
                print()

            # Do all switches

            # maximum of one switch can be made per round, stop checking if switch was found
            # ss is never switched, no need to check
            # if rr is switched, the switch is already taken care of when checking jj or kk, no need to check again

            # if node was switched, switch their parents
            if aa != best_top[0][0].index:
                # save indices of node jj
                aa_parent = ft.nodes[aa].parent

                # change indices of node jj to node from better topology
                ft.nodes[aa].parent = best_top[0][0].parent

                # find the node from the better topology and change indices to the ones from jj
                ft.nodes[best_top[0][0].index].parent = aa_parent

                # swap indices
                ft.nodes[aa].index = best_top[0][0].index
                ft.nodes[best_top[0][0].index].index = aa

                # swap positions in node list
                node_aa = ft.nodes[aa]
                ft.nodes[aa] = ft.nodes[best_top[0][0].index]
                ft.nodes[best_top[0][0].index] = node_aa

                if ft.verbose == 1:
                    print("swapped nodes1", ft.nodes[aa].index, ft.nodes[bb].index, ft.nodes[cc].index, ft.nodes[dd].index)
                    print()

            elif bb != best_top[0][1].index:
                # save indices of node kk
                bb_parent = ft.nodes[bb].parent

                # change indices of node kk to node from better topology
                ft.nodes[bb].parent = best_top[0][1].parent

                # find the node from the better topology and change indices to the ones from kk
                ft.nodes[best_top[0][1].index].parent = bb_parent

                # swap indices
                ft.nodes[bb].index = best_top[0][1].index
                ft.nodes[best_top[0][1].index].index = bb

                # swap positions in node list
                node_bb = ft.nodes[bb]
                ft.nodes[bb] = ft.nodes[best_top[0][1].index]
                ft.nodes[best_top[0][1].index] = node_bb

                if ft.verbose == 1:
                    print("swapped nodes2", ft.nodes[aa].index, ft.nodes[bb].index, ft.nodes[cc].index, ft.nodes[dd].index)
                    print()

        # Recompute profiles of internal nodes
        for node in ft.nodes:
            if node.leaf:  # skip all leaf nodes
                continue
            node.profile = average_profile([ft.nodes[node.leftchild], ft.nodes[node.rightchild]], ft.lambda1)
    # best_top = MinimizedEvolution(ft, nodes[0], nodes[1], nodes[8], nodes[11])
    # print("Best topology", best_top[0][0].index, best_top[0][1].index, best_top[1][0].index, best_top[1][1].index)

    return ft.nodes


def MinimizedEvolution(ft: Tree, n1, n2, n3, n4):
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
        dist_a = JC_distance(util.uncorrected_distance(ft, [options[ii][0][0], options[ii][0][1]]))
        dist_b = JC_distance(util.uncorrected_distance(ft, [options[ii][1][0], options[ii][1][1]]))
        # print(uncorrected_distance(nodes, [options[ii][0][0], options[ii][0][1]], lambda1))
        # if ii == 0:
        #     join = [options[ii][0][0], options[ii][0][1]]
        #     indices = [join[0].index, join[1].index]
        #     print(nodes[indices[0]].profile)
        #     print( nodes[indices[1]].profile)
        #     profiles = [nodes[indices[0]].profile, nodes[indices[1]].profile]
        #     value = 0
        #     for L in range(len(profiles)):
        #         for a in range(4):
        #             for b in range(4):
        #                 if profiles[0][L][a] > 0 and profiles[1][L][a] > 0:
        #                     P = 0
        #                 else:
        #                     P = 1

        #                 print(profiles[0][L])
        #                 value += profiles[0][L][a] * profiles[1][L][b] * P   
        #             print(value)        

        #     print(value / len(profiles))

        #     del_ij = profileDistanceNew([nodes[indices[0]].profile, nodes[indices[1]].profile])
        #     print(del_ij)

        # print(dist_a)
        dist_options.append(dist_a + dist_b)
        # print(uncorrected_distance([options[ii][0][0].profile, options[ii][0][1].profile]))
        # print(uncorrected_distance([options[ii][1][0].profile, options[ii][1][1].profile]))

    # Choose the topology with the minimized criterion
    best_topology = options[dist_options.index(min(dist_options))]

    return best_topology
