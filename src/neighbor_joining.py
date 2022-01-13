"""

"""
import sys

from typing import Tuple

from src.tree import Tree
from src.node import Node
import src.util as util

def minimize_nj_criterion(ft: Tree) -> Tuple[tuple, Node]:
    """Returns i,j for which d(i, j) - r(i) -r(j) is minimal and corresponding new node of merging i,j

    Args:
        ft (Tree): Tree object
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

            ft.update_T()
            criterion = nj_criterion(ft, i, j)
            print("criterion ", i.index, j.index, criterion)

            if criterion < min_dist:  # if best join for now
                min_dist = criterion
                best_join = (i, j)  # saves best joining nodes

    return best_join


def nj_criterion(ft: Tree, i: Node, j: Node) -> float:
    """
    Calculates the Neighbour Joining criterion between two nodes, based on Supplement 1.

    Args:
        ft: Tree
        i: Node
        j: Node

    Returns:
        criterion(float)
    """
    if i.leaf and j.leaf:
        criterion = util.uncorrected_distance(ft, [i, j]) - util.out_distance_new(ft, i) - util.out_distance_new(ft, j)
    elif i.leaf:
        criterion = util.uncorrected_distance(ft, [j.leftchild, j.rightchild, i]) - util.out_distance_new(ft,
                                                                                                i) - util.out_distance_new(
            ft, j)
    elif j.leaf:
        criterion = util.uncorrected_distance(ft, [i.leftchild, i.rightchild, j]) - util.out_distance_new(ft,
                                                                                                i) - util.out_distance_new(
            ft, j)
    else:
        criterion = util.uncorrected_distance(ft, [i.leftchild, i.rightchild, j]) - util.out_distance_new(ft,
                                                                                                i) - util.out_distance_new(
            ft, j)

    return criterion
