"""
Maybe some comments about implementation here.

References to page numbers in this code are referring to the paper:
[1] Price at al. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix.
    Molecular Biology and Evolution, vol 26 (7). 2009.
    Paper can be found on https://pubmed.ncbi.nlm.nih.gov/19377059/
"""
import math
import pprint as pp # For pretty printing (Replace with own code before submission)
import sys

from src.node import Node
from src.node import Node1

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
    i = 0
    for seq in sequences.keys():
        node = Node(seq, i, sequences[seq])
        i +=1
        nodes.append(node)

    # Actual first step : Unique sequences ( page 1646, do later )

    # Step 1 of algorithm : Create total profile T
    T = Profile(nodes)
    # print("Total profile T")
    # pp.pprint(T)

    # Step 2 : Top hits sequence
    # Skip for now

    # Step 3 : Create initial topology
    CreateInitialTopology(nodes)     
  
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

def uncorrectedDistance(profile: list) -> float: 
    """Calculate uncorrected distance
    
    The fraction of position that differ by using profiles. It's also #differences/sequence length

    Args:
        profile (list): The profile matrix of sequences of internal nodes containing ratios

    Returns:
        (float): the distance between internal nodes
    """
    k = len(profile)
    differ = 0
    for i in range(k):
        if 1.0 not in profile[i]:
            differ += 1
    fraction = differ / k
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
        for j in range(i+1, len(children)):
            combi.append( [v1, children[j]])
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
        seqDis.append(ham[i]/int(k))
    return seqDis

def out_distance(i, nodes):
    """Calculates r distance of i : r(i)

    Args:
        i (Node) : 
    """
    active_nodes = 1  # i is always an active node
    dist_to_others = 0
    for j in nodes:
        if j.name == i:
            continue
        if not j.active:
            continue

        active_nodes += 1

        profile_i_j = Profile([i, j])

        dist_to_others += uncorrectedDistance(profile_i_j)

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
    for i  in active_nodes:
        for j in active_nodes:
            if i == j:
                continue
            temp_profile_new_node = Profile([i, j]) #calculates profile of potential merge
            criterion = uncorrectedDistance(temp_profile_new_node) - out_distance(i, nodes) - out_distance(j, nodes)
            if criterion < min_dist: #if best join for now
                profile_new_node = temp_profile_new_node #sets profile of new node to profile of best join
                min_dist = criterion
                best_join = (i, j) #saves best joining nodes


    #save just calculated profile of joining nodes to a Node with name containing both joined nodes and make this new node active
    #we should probably change the class Node as the sequence is not known for the merged nodes. I just made a beun oplossing. Don't know if it's good enough
    new_node = Node(str(best_join[0].name) + '&' + str(best_join[1].name), int(index), 'nosequence') 
    new_node.profile = profile_new_node
    new_node.active = True
    #add indices of left child, right child 
    new_node.leftchild = best_join[0].index
    new_node.rightchild = best_join[1].index
    print("Minimized distance = ", min_dist, "of nodes ", best_join[0].name, best_join[1].name)
    return best_join, new_node

def CreateInitialTopology(nodes):
    for i in range(len(nodes)-1):
        minimized_join, new_node = minimize_nj_criterion(nodes, len(nodes))
        # make joined nodes inactive
        nodes[int(minimized_join[1].index)].active = False
        nodes[int(minimized_join[0].index)].active = False
        # append the newly joined node to list of nodes
        nodes.append(new_node)
        
        
        
        print("Merged nodes to: " + new_node.name)
        print("left child: " + str(nodes[len(nodes)-1].leftchild))
        print("right child: " + str(nodes[len(nodes)-1].rightchild))
        # I don't know what to return and save so heelpppp

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
    jd_d = -0.75 * math.log(1 - (4/3) * d_u)

    # For sequences that do not overlap, FastTree uses a max distance of 3.0
    if jd_d > max_distance:
        return max_distance
    
    return jd_d
