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

from src.node import Node
from src.tree import Tree
import src.tophits as tophits


def fast_tree(args: argparse.Namespace, sequences: dict) -> str:
    """FastTree Algorithm.

    Args:
        args: Namespace of the argparse library, containing all user inputs
        sequences (dict): Mapping of sequences to their names, as was provided in the program's input

    Returns:
        (str): A phylogenetic tree in Newick format.
    """

    print("The sequences entered into the program : ")
    pp.pprint(sequences)

    # Create list of Nodes representing the sequences
    nodes = []
    for ii, seq in enumerate(sequences.keys()):
        node = Node(seq, ii, sequences[seq])
        node.leaf = True
        nodes.append(node)

    # Create Tree object
    ft = Tree(nodes, args)

    # Actual first step : Unique sequences ( page 1646, do later )

    # Heuristics for neighbor joining (with or without top-hits)

    # Top hits sequence
    if not ft.no_top_hits:
        tophits.top_hits_init(ft)

    # FastNJ best joins
    fast_nj_best_joins(ft)

    # End of neighbor joining heuristic preparation

    # Create initial topology
    CreateInitialTopology(ft)
    print("Initial topology is created:")
    print(len(ft.nodes))
    ft.newick_str()

    # Nearest Neighbor Interchanges
    NNI(ft)
    print("NNI is finished:")

    # Final step: print tree topology as Newick string
    return ft.newick_str()


def fast_nj_best_joins(ft: Tree) -> None:
    """
    The key idea in FastNJ is to store the best join for each node.

    The best join for each leaf is determined before the joins begin, and the best join for
    each new interior node is determined when that node is created. When searching for the best join overall,
    FastNJ considers only best join for each node, or n candidates. Thus, FastNJ requires a total of O(N^2)
    time.

    Note: when using top-hits, then at the beginning the best join for each node is found in their top-hits lists.
    """
    if ft.no_top_hits:
        for i in ft.nodes:
            best_join = (sys.maxsize / 2, 0)  # (distance, node index)
            for j in ft.nodes:
                if i.index == j.index:
                    continue

                criterion = uncorrected_distance(ft, [i, j]) - out_distance_new(ft, i) - out_distance_new(ft, j)

                if criterion < best_join[0]:  # if best join for now
                    best_join = (criterion, j.index)

            # Found the best join for i, now update it
            i.best_join = best_join

    else:
        for node in ft.nodes:
            # Nice, the best join for each node is found in their top-hits list already!
            # But wait, while computing the top hits of A, we may discover that A,B is a better join than B,best(B).

            # So if the best_join was already set when A,B is a better join than B,best(B) was true, move on to the next
            if node.best_join is not None:
                continue

            # Okay after that check, we can use tophits
            best_join_dist, best_join = node.tophits.list[0]
            node.best_join = (best_join_dist, best_join)

            best_B_dist, best_B = ft.nodes[best_join].tophits.list[0]

            # A, B is a better join than B, Best(B)
            if best_B_dist > best_join_dist:
                ft.nodes[best_join].best_join = (best_join_dist, node.index)

    if ft.verbose == 1:
        for node in ft.nodes:
            print('FastNJ best join for Node ', node.index, 'is', node.best_join)


def tophits_new_node(new_node: Node) -> PriorityQueue:
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
        top_hits.put(th)

    for th in new_node.rightchild.tophits.queue:
        top_hits.put(th)

    return top_hits


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


def profile_distance_new(profiles: list) -> float:
    """Calculates “profile distance” that is the average distance between profile characters over all positions
    
    input: 
        list of profiles of which you want to calculate profile distance (i,j)
    
    returns:
        profile distance ∆(i, j)
    """
    value = 0
    for L in range(len(profiles[0])):
        for a in range(4):
            for b in range(4):
                if profiles[0][L][a] > 0 and profiles[1][L][a] > 0:
                    P = 0
                else:
                    P = 1
                value += profiles[0][L][a] * profiles[1][L][b] * P

    return value / len(profiles)


def findAllKids(nodes: list, nodeIndex: int) -> list:
    '''Finds all children of input nodeIndex
        
        Input:
            All nodes (list): list of all nodes
            Index of node (int): index of node of which you want to know the children
        
        Returns:
            All children (list): list with indexes of all children of specific node
    '''
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


def uncorrected_distance(ft: Tree, join: list):
    """Uncorrected distance between joined nodes (node_ij) and other nodes (node_k)
    du(ij, k) = ∆(ij, k) − u(ij) − u(k)
    args;
        list with input nodes
        indices of nodes which are joined
    returns:
        uncorrected Distance 
    """
    if len(join) == 2:
        indices = [join[0].index, join[1].index]
        # du(i, j) = ∆(i, j)−u(i)−u(j), where u(i) = 0 and u(j) = 0 as i and j are leaves
        del_ij = profile_distance_new([ft.nodes[indices[0]].profile, ft.nodes[indices[1]].profile])
        return del_ij
    if len(join) == 3:
        indices = [join[0], join[1], join[2].index]
        del_ijk = profile_distance_nodes(ft, indices)
        u_ij = updistance(ft.nodes, [indices[0], indices[1]])
        u_k = updistance(ft.nodes, [indices[2]])
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
        return profile_distance_new([nodes[ijk[0]].profile, nodes[ijk[1]].profile]) / 2
    elif nodes[ijk[0]].leaf == True:
        return 0
    else:
        return profile_distance_new(
            [nodes[nodes[ijk[0]].rightchild].profile, nodes[nodes[ijk[0]].leftchild].profile]) / 2


def out_distance_new(ft: Tree, i: Node) -> float:
    """The average profile distance between a node and all other
       nodes can be inferred from the total profile T: r(i) = (n∆(i, T) − ∆(i, i) − (n − 1)u(i) + u(i) − sum u(j))/(n-2)
        
    Args:
        active nodes; list of nodes; T total profile of current topology
    returns:
        out distance of one node 
    """

    N_active_nodes = 1  # i is always an active node; is 
    sumJ = 0
    for j in ft.nodes:
        if j.name == i:
            continue
        if not j.active:
            continue
        N_active_nodes += 1
        sumJ += updistance(ft.nodes, [j.index])

    # ∆(i, i) is the average distance between children of i, including self-comparisons.
    if i.leaf == True:
        d_ii = profile_distance_new([i.profile, i.profile]) # is always 0 I think but ok
        
    else:
        d_ii = profile_distance_new([ft.nodes[i.rightchild].profile, ft.nodes[i.leftchild].profile])
        # print('out_dist', d_ii)
    # d_ii = profile_distance_new([i.profile, i.profile])

    n_iT = N_active_nodes * profile_distance_new([i.profile, ft.T])
    n_u_i = (N_active_nodes - 1) * updistance(ft.nodes, [i.index])
    u_i = updistance(ft.nodes, [i.index])

    #(n∆(i, T) − ∆(i, i) − (n − 1)u(i) + u(i) − sum u(j))
    sum_du_ij = n_iT - d_ii - n_u_i + u_i - sumJ
    if N_active_nodes == 2:
        return sum_du_ij
    return sum_du_ij / (N_active_nodes - 2)


def profile_distance_nodes(ft: Tree, indices: list) -> float:
    """∆(ij, k)
        indices [i, j, k]
        ∆(ij, k) = λ∆(i, k) + (1 − λ)∆(j, k)
        profile of parent after joining i and j

        Args:
            ft: Tree object
            indices: index of nodes [i, j, k]
        Returns:
            float: Profile distance between joined nodes and other nodes
        """

    profile_dist = ft.lambda1 * \
                   (profile_distance_new([ft.nodes[indices[0]].profile, ft.nodes[indices[2]].profile])) + \
                   (1 - ft.lambda1) * \
                   (profile_distance_new([ft.nodes[indices[1]].profile, ft.nodes[indices[2]].profile]))
    return profile_dist


def minimize_nj_criterion(ft: Tree, index: int) -> Tuple[tuple, Node]:
    """Returns i,j for which d(i, j) - r(i) -r(j) is minimal and corresponding new node of merging i,j
    
    Args:
        ft (Tree): Tree object
        index (int): index of new node
    Returns:
        best_joins (list): list containing the joined node objects
        new_node (Node): Node object of the new node
    """
    active_nodes = []
    for node in ft.nodes:
        if node.active:
            active_nodes.append(node)

    min_dist = sys.float_info.max
    best_join = (0, 0)
    for i in active_nodes:
        for j in active_nodes:
            if i == j:
                continue
            temp_profile_new_node = average_profile([i, j], ft.lambda1)  # calculates profile of potential merge
            ft.update_T()
            if i.leaf and j.leaf:
                criterion = uncorrected_distance(ft, [i, j]) - out_distance_new(ft, i) - out_distance_new(ft, j)
            elif i.leaf:
                criterion = uncorrected_distance(ft, [j.leftchild, j.rightchild, i]) - out_distance_new(ft, i) - out_distance_new(ft, j)
            elif j.leaf:
                criterion = uncorrected_distance(ft, [i.leftchild, i.rightchild, j]) - out_distance_new(ft, i) - out_distance_new(ft, j)
            else:
                criterion = uncorrected_distance(ft, [i.leftchild, i.rightchild, j]) - out_distance_new(ft, i) - out_distance_new(ft, j)

            if criterion < min_dist:  # if best join for now
                profile_new_node = temp_profile_new_node  # sets profile of new node to profile of best join
                min_dist = criterion
                best_join = (i, j)  # saves best joining nodes

    # Save just calculated profile of joining nodes to a Node with name containing both joined nodes and make this new node active
    # we should probably change the class Node as the sequence is not known for the merged nodes. I just made a beun oplossing. Don't know if it's good enough
    new_node = Node(str(best_join[0].name) + '&' + str(best_join[1].name), int(index), 'nosequence')
    new_node.profile = profile_new_node

    # add indices of left child, right child
    new_node.leftchild = best_join[0].index
    new_node.rightchild = best_join[1].index
    ft.nodes[best_join[0].index].parent = new_node.index
    ft.nodes[best_join[1].index].parent = new_node.index

    # # Update top-hits for new node
    # new_node.tophits = tophits_new_node(new_node)

    # print("Minimized distance = ", min_dist, "of nodes ", best_join[0].name, best_join[1].name)
    return best_join, new_node


def create_join(ft: Tree, minimized_join, new_node):
    """
    Join two nodes.


    """
    # make joined nodes inactive
    ft.nodes[int(minimized_join[1].index)].active = False
    ft.nodes[int(minimized_join[0].index)].active = False

    BranchLength(ft, minimized_join)

    # append the newly joined node to list of nodes
    ft.nodes.append(new_node)

    if ft.verbose == 1:
        print(new_node.name, ft.nodes[minimized_join[0].index].branchlength)
        print("Merged nodes to: " + new_node.name)
        print("left child: " + str(ft.nodes[len(ft.nodes) - 1].leftchild))
        print("right child: " + str(ft.nodes[len(ft.nodes) - 1].rightchild))

    # Recalculate total profile T
    # When we do a join, we also need to update the total profile (the average over all active nodes). To
    # compute this profile takes O(nLa) time. However, we can subtract the joined nodes and add the new
    # node to the total profile in O(La) time. (FastTree) recomputes the total profile from scratch every 200
    # iterations to avoid roundoff errors from accumulating, where the choice of 200 is arbitrary. This adds 
    # another O(N 2La/200) work.)
    ft.update_T()

    # TODO And update top-hits of new join

    return ft.nodes


def CreateInitialTopology(ft: Tree) -> list:
    """
    Create the initial topology given a list with all input nodes. 

    Args:
        ft (Tree): Tree object
    """

    if ft.no_top_hits:
        # Original number of nodes (nr of sequences)
        nr_leafs = len(ft.nodes)
        for i in range(nr_leafs - 1):
            minimized_join, new_node = minimize_nj_criterion(ft, len(ft.nodes))
            create_join(ft, minimized_join, new_node)

    else:
        # With top hits, we can use We use the FastNJ and local hill-climbing heuristics,
        # with the further restriction that we consider only the top m candidates at each step

        active_nodes = []
        for node in ft.nodes:
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

            i = ft.nodes[candidate[0]]
            j = ft.nodes[candidate[1]]

            criterion = uncorrected_distance(ft, [i, j]) - out_distance_new(ft, i) - out_distance_new(ft, j)

            if criterion < min_dist:  # if best join for now
                best_candidate = candidate
                min_dist = criterion

        print("The best candidate = ", best_candidate)

        # Given this candidate join, we do a local hill-climbing search for a better join
        best_join = local_hill_climbing_tophits(ft.nodes, best_candidate, min_dist)

        # Make the join
        # @TODO create_join(nodes,)


def local_hill_climbing_tophits(nodes: list, best_candidate: tuple, best_dist: float) -> tuple:
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
        fraction = uncorrected_distance(ft, [ft.nodes[n1], ft.nodes[n2]])
        ft.nodes[n1].branchlength = JC_distance(fraction)
        ft.nodes[n2].branchlength = JC_distance(fraction)
    # connect single leaf with other branch
    elif n1 < nr_leafs <= n2:
        d12 = uncorrected_distance(ft, [ft.nodes[n1], ft.nodes[ft.nodes[n2].leftchild]])
        d13 = uncorrected_distance(ft, [ft.nodes[n1], ft.nodes[ft.nodes[n2].rightchild]])
        d23 = uncorrected_distance(ft, [ft.nodes[ft.nodes[n2].leftchild], ft.nodes[ft.nodes[n2].rightchild]])
        ft.nodes[n1].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23)) / 2
        ft.nodes[n2].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23)) / 2
    # connect single leaf with other branch
    elif n2 < nr_leafs <= n1:
        d12 = uncorrected_distance(ft, [ft.nodes[n2], ft.nodes[ft.nodes[n1].leftchild]])
        d13 = uncorrected_distance(ft, [ft.nodes[n2], ft.nodes[ft.nodes[n1].rightchild]])
        d23 = uncorrected_distance(ft, [ft.nodes[ft.nodes[n2].leftchild], ft.nodes[ft.nodes[n2].rightchild]])
        ft.nodes[n1].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23)) / 2
        ft.nodes[n2].branchlength = (JC_distance(d12) + JC_distance(d13) - JC_distance(d23)) / 2
    # connect two branches
    else:
        d13 = uncorrected_distance(ft, [ft.nodes[ft.nodes[n1].leftchild], ft.nodes[ft.nodes[n2].leftchild]])
        d14 = uncorrected_distance(ft, [ft.nodes[ft.nodes[n1].leftchild], ft.nodes[ft.nodes[n2].rightchild]])
        d23 = uncorrected_distance(ft, [ft.nodes[ft.nodes[n1].rightchild], ft.nodes[ft.nodes[n2].leftchild]])
        d24 = uncorrected_distance(ft, [ft.nodes[ft.nodes[n1].rightchild], ft.nodes[ft.nodes[n2].rightchild]])
        d12 = uncorrected_distance(ft, [ft.nodes[ft.nodes[n1].leftchild], ft.nodes[ft.nodes[n1].rightchild]])
        d34 = uncorrected_distance(ft, [ft.nodes[ft.nodes[n2].leftchild], ft.nodes[ft.nodes[n2].rightchild]])
        ft.nodes[n1].branchlength = (JC_distance(d13) + JC_distance(d14) + JC_distance(d23) + JC_distance(d24)) / 4 - (
                JC_distance(d12) + JC_distance(d34)) / 2
        ft.nodes[n2].branchlength = (JC_distance(d13) + JC_distance(d14) + JC_distance(d23) + JC_distance(d24)) / 4 - (
                JC_distance(d12) + JC_distance(d34)) / 2

def NNI(ft: Tree):
    print('start NNI')
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

            # error because index of node is equal to index of parent -> it seems like something goes wrong with indices after swapping
            print('NNI compares', ft.nodes[aa].index, ft.nodes[bb].index, ft.nodes[cc].index, ft.nodes[dd].index)

            # For each possible combination of fixed nodes, find the best topology
            best_top = MinimizedEvolution(ft, ft.nodes[aa], ft.nodes[bb], ft.nodes[cc], ft.nodes[dd])
            print('NNI best topology', best_top[0][0].index, best_top[0][1].index, best_top[1][0].index,
                  best_top[1][1].index)

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

                print("swapped nodes1", ft.nodes[aa].index, ft.nodes[bb].index, ft.nodes[cc].index, ft.nodes[dd].index)

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

                print("swapped nodes2", ft.nodes[aa].index, ft.nodes[bb].index, ft.nodes[cc].index, ft.nodes[dd].index)

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
        dist_a = JC_distance(uncorrected_distance(ft, [options[ii][0][0], options[ii][0][1]]))
        dist_b = JC_distance(uncorrected_distance(ft, [options[ii][1][0], options[ii][1][1]]))
        dist_options.append(dist_a + dist_b)

    # Choose the topology with the minimized criterion
    best_topology = options[dist_options.index(min(dist_options))]
    print("dist_option", dist_options)
    return best_topology
