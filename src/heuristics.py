"""
File representing the FastNJ, Local Hill-Climbing and Top-hits heuristics.
All relevant classes or functions will be defined here, with possible references to utility functions i.e. distances.

References to page numbers in this code are referring to the paper or its supplementary material:
[1] Price at al. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix.
    Molecular Biology and Evolution, vol 26 (7). 2009.
    Paper can be found on https://pubmed.ncbi.nlm.nih.gov/19377059/
"""
import math
import sys
from operator import itemgetter

from src.tree import Tree
from src.node import Node
import src.neighbor_joining as neighbor_joining
import src.util as util


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

    @classmethod
    def top_hits_init(cls, ft: Tree) -> list:
        """
        Create top-hits list for all N nodes before joining.

        Paper excerpt:
        Before doing any joins, FastTree estimates these lists for all N sequences by assuming that,
        if A and B have similar sequences, then the top-hits lists of A and B will largely overlap.
        More precisely, FastTree computes the 2m top hits of A, where the factor of two is a safety factor.
        Then, for each node B within the top m hits of A that does not already have a top-hits list,
        FastTree estimates the top hits of B by comparing B to the top 2m hits of A.
        """
        # For each node, FastTree records a top-hits node.
        # Should this be randomized ? random.shuffle(nodes)
        for A in ft.nodes:

            # Since we infer top-hits list of B through A, a node B might already have a tophits list, so skip.
            if A.tophits is not None:
                continue

            # Compute the 2m tophits of a node A (2 is a safety factor)
            A.tophits = TopHits(ft.m)

            # Get top-hits sorted
            for node in ft.nodes:
                if A.index == node.index:
                    continue

                # closest is defined according to the Neighbor-Joining criterion
                criterion = neighbor_joining.nj_criterion(ft, A, node)

                A.tophits.list.append((criterion, node.index))

            # For a Node A, only the top 2m tophits are required
            A.tophits.list = sorted(A.tophits.list, key=itemgetter(0))
            A.tophits.list = A.tophits.list[:2 * ft.m]

            # Then, for each node B within the top m hits of A that does not already have a top-hits list,
            # FastTree estimates the top hits of B by comparing B to the top 2m hits of A.

            # For top m hits of A
            for ii in range(ft.m):
                # Make sure A has at least m hits
                if ii >= len(A.tophits.list) - 1:
                    break

                # top-hits are stored as tuple, (distance, node_index)
                B_index = A.tophits.list[ii][1]
                B = ft.nodes[B_index]

                # That does not already have a top-hits list
                if B.tophits is not None:
                    continue

                # Before FastTree estimates the top-hits of B from the top-hits of A,
                # FastTree requires that du(A,B) ≤ 0.75·du(A,H2m), where H2m is A’s 2m-th best hit. (See supplement)
                close_enough_factor = 0.75
                du_A_B = util.uncorrected_distance(ft, [A, B])
                H2m = ft.nodes[A.tophits.list[2 * ft.m - 1][1]]
                du_A_H2m = util.uncorrected_distance(ft, [A, H2m])

                if du_A_B > close_enough_factor * du_A_H2m:
                    # du(AB) wasn't smaller than or equal to 0.75·du(A,H2m) -> B wasn't close enough for this heuristic.
                    break

                # Top hits of B are found in the top 2m hits of A
                B.tophits = TopHits(ft.m)
                for jj in range(2 * ft.m):
                    # Make sure A has a hit
                    if jj > len(A.tophits.list) - 1:
                        break
                    node_index = A.tophits.list[jj][1]
                    node = ft.nodes[node_index]

                    # Don't add yourself
                    if B.index == node.index:
                        node = A
                        # continue

                    # closest is defined according to the Neighbor-Joining criterion
                    criterion = neighbor_joining.nj_criterion(ft, B, node)

                    B.tophits.list.append((criterion, node.index))

        # Finally, some nodes will have been considered as A, having a tophits list of length 2m,
        # And some Nodes B, inferred from A, will have smaller ones.
        # "For each node, FastTree records a top-hits list: the nodes that are the closest m neighbors of that node"
        for node in ft.nodes:
            node.tophits.list = sorted(node.tophits.list, key=itemgetter(0))
            node.tophits.list = node.tophits.list[:ft.m]

            if ft.verbose == 1:
                print("Tophits of node", node.index)
                for th in node.tophits.list:
                    print(th)

                print()

        return ft.nodes

    @classmethod
    def tophits_refresh(cls, ft: Tree, outdated_node: Node) -> list:
        """Refresh Top-Hits List of a Node.

        "We recompute the top-hit list for the new joined node and we update the top-hit lists of the
        new node’s top hits."

        Refreshing a top-hit list takes O(nLa + m2La) = O(N La) time and ensures that the top-hit lists
        of O(m = √N ) other nodes reach size m. Thus, the refreshes take a total of O(N √N La) time.

        Args:
            ft (Tree) : The Tree Object
            outdated_node (Node) : Node with outdated top-hits list

        Returns:
            list : The Tophits list of the node to be updated
        """
        # To refresh, we compare the new node to all n − 1 other active nodes
        new_top_hits = []
        for node in ft.nodes:
            if outdated_node.index == node.index:
                continue

            if node.active:
                criterion = neighbor_joining.nj_criterion(ft, outdated_node, node)
                new_top_hits.append((criterion, node.index))

        new_top_hits = sorted(new_top_hits, key=itemgetter(0))

        # Then, we compare the close neighbors of the outdated node (the top m hits) to the top 2m hits of the
        # outdated node, and we update the close neighbors’ top-hit lists by merging.
        for i in range(ft.m):
            # New top hits has less than m entries, don't forget to refresh if the size is below 0.8m
            if i >= len(new_top_hits):
                break

            close_neighbor = ft.nodes[new_top_hits[i][1]]

            comparison_top_hits = []
            for j in range(0, 2 * ft.m):
                if j >= len(new_top_hits):  # No more than m + j top-hits found
                    break

                # Don't add yourself
                if i == j:
                    continue

                criterion = neighbor_joining.nj_criterion(ft, outdated_node, ft.nodes[new_top_hits[j][1]])
                comparison_top_hits.append((criterion, ft.nodes[new_top_hits[j][1]].index))

            # Update by merging best m candidates
            close_neighbor.tophits.list = sorted(comparison_top_hits, key=itemgetter(0))[:ft.m]

            if ft.verbose == 1:
                print("Refreshed the tophits of Node ", close_neighbor.index, 'to be', close_neighbor.tophits.list)
                print()

            # Set age of updated top-hits list to 0
            close_neighbor.tophits.age = 0

        # We save the top m hits of the new node's tophits
        new_top_hits = new_top_hits[:ft.m]

        if ft.verbose == 1:
            print("Refreshed the tophits of Node ", outdated_node.index, 'to be', new_top_hits)
            print()

        return new_top_hits

    def tophits_new_node(self, ft: Tree, new_node: Node) -> None:
        """
        After a join, FastTree computes the top-hits list for the new node in O(mLa) time
        by comparing the node to all entries in the top-hits lists of its children.

        Args:
            ft (Tree) : The Tree Object
            new_node (Node) : The newly created inner node after a join

        Returns:
            None : the new Node gets updated
        """
        new_node.tophits = TopHits(ft.m)

        # The top hits list of a newly joined node AB is the top m hits from the tophits of A and B
        # Merge the tophits lists of A and B to get up to 2(m - 1) candidates, and remove duplicates
        tophits_A = ft.nodes[new_node.leftchild].tophits.list
        tophits_B = ft.nodes[new_node.rightchild].tophits.list

        tophits_AB = tophits_A + list(set(tophits_B) - set(tophits_A))

        # Note that we must remove A and B from this new list (they are joined and no longer active)
        tophits_AB_cleaned = []
        for hit in tophits_AB.copy():
            if hit[1] != new_node.leftchild and hit[1] != new_node.rightchild:  # No removal required
                tophits_AB_cleaned.append(hit)

        # Doing it this way instead of .remove in the loop reduces time complexity from O(N^2) to O(N),
        # but takes up a O(2m) more memory
        tophits_AB = tophits_AB_cleaned

        # And store the top m hits
        tophits_AB = sorted(tophits_AB, key=itemgetter(0))
        tophits_AB = tophits_AB[:ft.m]

        # Set the age of a new node to one plus the maximum of its children's ages
        new_node.tophits.list = \
            1 + max(ft.nodes[new_node.leftchild].tophits.age, ft.nodes[new_node.rightchild].tophits.age)

        # If after a join, either
        # i) the top-hit list has shrunk too much (below 0.8m, where 0.8 is an arbitrary parameter),
        # or ii) the age is above 1 + log2(m)
        # then it does a refresh.
        if len(tophits_AB) < (self.refreshFactor * ft.m) or new_node.tophits.age > (1 + math.log2(ft.m)):
            tophits_AB = self.tophits_refresh(ft, new_node)

        new_node.tophits.list = tophits_AB

        if ft.verbose == 1:
            print("Tophits of the new node ", new_node.index, '=', new_node.tophits.list)


def fastNJ_init(ft: Tree) -> None:
    """
    The key idea in FastNJ is to store the best join for each node.

    The best join for each leaf is determined before the joins begin, and the best join for
    each new interior node is determined when that node is created. When searching for the best join overall,
    FastNJ considers only best join for each node, or n candidates. Thus, FastNJ requires a total of O(N^2)
    time.

    Note: when using top-hits, then at the beginning the best join for each node is found in their top-hits lists.

    Args:
        ft (Tree): Tree object
    """
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
            if node.active:
                print('FastNJ best join for Node ', node.index, 'is', node.best_join)
        print()


def fastNJ_update(ft: Tree, node: Node):
    """Calculate the FastNJ Best-hit for a Node

    Args:
        ft (Tree): Tree Object
        node (Node): Newly created join, or Node with outdated Best-hit

    Returns:
        None -> the new_node object gets updated
    """
    if ft.verbose == 1:
        print('Old FastNJ best join for Node ', node.index, 'is', node.best_join)
        print()

    # Nice, the best join for each node is found in their top-hits list already!
    best_join_dist, best_join = node.tophits.list[0]

    # Refresh top-hits lists the "lazy" way, after coming across an inactive node
    if not ft.nodes[best_join].active:
        # Update with reference to active parent
        best_join = ft.nodes[best_join].parent
        best_join_dist = neighbor_joining.nj_criterion(ft, node, ft.nodes[best_join])

    node.best_join = (best_join_dist, best_join)

    # But if a node has no top-hits, this means we have reached the end of the program!
    if len(ft.nodes[best_join].tophits.list) == 0:
        if ft.verbose == 1:
            print("Newly created node ", best_join, " has no top-hits. This means this was the last join!")
            print()

        return
    best_B_dist, best_B = ft.nodes[best_join].tophits.list[0]

    # A, B is a better join than B, Best(B)
    if best_B_dist > best_join_dist:
        ft.nodes[best_join].best_join = (best_join_dist, node.index)

    if ft.verbose == 1:
        print('FastNJ best join for Node ', node.index, 'is', node.best_join)
        print()


def local_hill_climb(ft: Tree, best_candidate: tuple, best_dist: float) -> tuple:
    """
    Perform Local Hill Climbing with or without the top-hits heuristic.

    Without top-hits:
    Given an (arbitrary) node A, it will find the best join partner B for A, and then the best join partner C for B.
    If A=C, then A and B are each other’s best hit and it has reached a local optimum; otherwise it continues searching
    from (B,C). To avoid very poor local optima, local hill climb also adds a check
    to ensure that is not lengthening the tree. (If it is, it starts over with another node.)
    Local hill climb takes O(N 2 log N ) time.

    Using the top-hits heuristic, we only search within the top-hit lists rather than
    comparing the two nodes to all other active nodes.
    """

    for hit_idx, hit in enumerate(ft.nodes[best_candidate[0]].tophits.list):

        # Lazy: when it encounters a hit to a joined node, it replaces that with a hit to the active ancestor
        if not ft.nodes[hit[1]].active:
            # Update with reference to active parent
            best_join = ft.nodes[hit[1]].parent
            best_join_dist = neighbor_joining.nj_criterion(ft, ft.nodes[best_candidate[0]], ft.nodes[best_join])
            ft.nodes[best_candidate[0]].tophits.list[hit_idx] = (best_join_dist, best_join)

            # Update the hit we encountered with the updated version
            hit = ft.nodes[best_candidate[0]].tophits.list[hit_idx]

        if hit[0] < best_dist:
            best_candidate = (best_candidate[0], hit[1])
            best_dist = hit[0]

    for hit_idx, hit in enumerate(ft.nodes[best_candidate[1]].tophits.list):

        # Lazy: when it encounters a hit to a joined node, it replaces that with a hit to the active ancestor
        if not ft.nodes[hit[1]].active:
            # Update with reference to active parent
            best_join = ft.nodes[hit[1]].parent
            best_join_dist = neighbor_joining.nj_criterion(ft, ft.nodes[best_candidate[1]], ft.nodes[best_join])
            ft.nodes[best_candidate[1]].tophits.list[hit_idx] = (best_join_dist, best_join)

            # Update the hit we encountered with the updated version
            hit = ft.nodes[best_candidate[1]].tophits.list[hit_idx]

        if hit[0] < best_dist:
            best_candidate = (hit[1], best_candidate[1])
            best_dist = hit[0]

    return ft.nodes[best_candidate[0]], ft.nodes[best_candidate[1]]
