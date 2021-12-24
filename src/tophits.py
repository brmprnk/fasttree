"""


References to page numbers in this code are referring to the paper or its supplementary material:
[1] Price at al. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix.
    Molecular Biology and Evolution, vol 26 (7). 2009.
    Paper can be found on https://pubmed.ncbi.nlm.nih.gov/19377059/
"""

import math
from operator import itemgetter


class TopHits:
    """
    Class attached to each node to implement top-hits heuristic.

    Note on implementation:
        This heuristic seems to lend itself well to using a PriorityQueue.
        However, it is often required to get the first m elements from the Queue.
        In Python, a priority queue is heap-based,
        and the only way to efficiently get the top m items is to pop m times, time complexity O(m log n).
        A PriorityQueue object in Python can also not efficiently pop without removing from the list,
        so inserting back would also take O(m log m) time.
        Therefore, a simple list is used, and just like the FastTree algorithm as described in the paper,
        this list is periodically sorted in O(N log N) time.
    """

    def __init__(self, m):
        self.list = []

        self.m = m
        self.age = 0

        # A top-hit list is too short if it has less than 0.8m entries,
        # where 0.8 is a tuning parameter and m is the size of the top-hit lists (see Supplement.)
        self.refreshFactor = 0.8


def top_hits_init(nodes: list, verbose: int=0) -> list:
    """
    Create top-hits list for all N nodes before joining.

    Paper excerpt:
    Before doing any joins, FastTree estimates these lists for all N sequences by assuming that,
    if A and B have similar sequences, then the top-hits lists of A and B will largely overlap.
    More precisely, FastTree computes the 2m top hits of A, where the factor of two is a safety factor.
    Then, for each node B within the top m hits of A that does not already have a top-hits list,
    FastTree estimates the top hits of B by comparing B to the top 2m hits of A.
    """
    N = len(nodes)
    m = int(math.sqrt(N))

    # For each node, FastTree records a top-hits node.
    # Should this be randomized ? random.shuffle(nodes)
    for A in nodes:

        # Since we infer top-hits list of B through A, a node B might already have a tophits list, so skip.
        if A.tophits is not None:
            continue

        # Compute the 2m tophits of a node A (2 is a safety factor)
        A.tophits = TopHits(m)

        # Get top-hits sorted
        for node in nodes:
            if A.index == node.index:
                continue

            temp_profile_new_node = averageProfile([A, node]) #calculates profile of potential merge

            # closest is defined according to the Neighbor-Joining criterion
            criterion = uncorDistance(
                [temp_profile_new_node, A.profile]) - out_distance(A, nodes) - out_distance(node, nodes)
            
            A.tophits.list.append((criterion, node.index))

        # For a Node A, only the top 2m tophits are required
        A.tophits.list = sorted(A.tophits.list, key=itemgetter(0))
        A.tophits.list = A.tophits.list[:2 * m]

        # Then, for each node B within the top m hits of A that does not already have a top-hits list,
        # FastTree estimates the top hits of B by comparing B to the top 2m hits of A.

        # For top m hits of A
        for ii in range(m):  
            # Make sure A has at least m hits
            if ii >= len(A.tophits.list) - 1: 
                break

            # top-hits are stored as tuple, (distance, node_index)
            B_index = A.tophits.list[ii][1]
            B = nodes[B_index]

            # That does not already have a top-hits list
            if B.tophits is not None:
                continue

            # Before FastTree estimates the top-hits of B from the top-hits of A, 
            # FastTree requires that du(A,B) ≤ 0.75·du(A,H2m), where H2m is A’s 2m-th best hit. (See supplement)
            close_enough_factor = 0.75
            du_A_B = uncorDistance([A.profile, B.profile])
            H2m = nodes[A.tophits.list[2 * m - 1][1]]
            du_A_H2m = uncorDistance([A.profile, H2m.profile])


            if du_A_B > close_enough_factor * du_A_H2m:
                # du(AB) was not smaller than or equal to 0.75·du(A,H2m), so B was not close enough to perform this heuristic.
                break
            
            # Top hits of B are found in the top 2m hits of A
            B.tophits = TopHits(m)
            for jj in range(2 * m):

                # Make sure A has a hit
                if jj > len(A.tophits.list) - 1:
                    break
                node_index = A.tophits.list[jj][1]
                node = nodes[node_index]

                # Don't add yourself
                if B.index == node.index:
                    continue

                temp_profile_new_node = averageProfile([A, node]) #calculates profile of potential merge

                # closest is defined according to the Neighbor-Joining criterion
                criterion = uncorDistance(
                    [temp_profile_new_node, A.profile]) - out_distance(A, nodes) - out_distance(node, nodes)
                
                B.tophits.list.append((criterion, node.index))
        
    # Finally, some nodes will have been considered as A, having a tophits list of length 2m,
    # And some Nodes B, inferred from A, will have smaller ones.
    # But, "For each node, FastTree records a top-hits list: the nodes that are the closest m neighbors of that node"
    for node in nodes:
        node.tophits.list = node.tophits.list[:m]

        if verbose == 1:
            print("Tophits of node", node.index)
            for th in node.tophits.list:
                print(th)

            print()

    return nodes

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