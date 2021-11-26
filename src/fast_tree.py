"""
Maybe some comments about implementation here.

References to page numbers in this code are referring to the paper:
[1] [1] Price at al. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix.
    Molecular Biology and Evolution, vol 26 (7). 2009.
    Paper can be found on https://pubmed.ncbi.nlm.nih.gov/19377059/
"""
import math
import pprint  # For pretty printing (Replace with own code before submission)

def fast_tree(sequences) -> str:
    """FastTree Algorithm.

    Args:
        sequences (dict): Mapping of sequences to their names, as was provided in the program's input

    Returns:
        (str): A phylogenetic tree in Newick format.
    """

    print("The sequences entered into the program : ")
    pprint.pprint(sequences)

    # Actual first step : Unique sequences ( page 1646, do later )

    # Step 1 of algorithm : Create total profile T
    all_sequences = list(sequences.values())
    T = Profile(all_sequences)

    # Step 2 : Top hits sequence
    # Skip for now

    # Step 3 : Create initial topology

    # calculate Profile for internal nodes
    Hypothetical_internal_node = [sequences['2'] , sequences['4'], sequences['5']]
    Node_profile = Profile(Hypothetical_internal_node)
    Sequence_distance = uncorrectedDistance(Node_profile)

    combi = makeCombisofChildren(Hypothetical_internal_node)
    Sequence_distance2 = SequenceDistance(combi, len(sequences['2']))  #calculates list of all combi's of children
    # print('Profile between 2 and 4 is: ' + str(Node_profile))
    # print('Their sequence distance (uncorrected distance) is: ' + str(Sequence_distance))
    # print('Their sequence distance (uncorrected distance) is: ' + str(Sequence_distance2)) 


# Calculate Profile of internal nodes
def Profile(internalNodes):
    columns = [''.join(seq) for seq in zip(*internalNodes)]
    return [[float(col.count(base)) / float(len(col)) for base in 'ACGT'] for col in columns]

def uncorrectedDistance(profile): #i.e. the fraction of positions that differ aka #difference/sequence length by using profiles
    k = len(profile)
    differ = 0
    for i in range(k):
        if 1.0 not in profile[i]:
            differ += 1 
    fraction = differ / k                 
    return fraction

def makeCombisofChildren(children): #make combinations to calculate the hamming distance between children
    combi = []
    for i, v1 in enumerate(children):
        for j in range(i+1, len(children)):
            combi.append( [v1, children[j]])
    return combi


def HammingDistance(combi): #calculate hamming distance between combinations of children
    distance = []
    for j in range(len(combi)):    
        hammingdistance = 0
        for i in range(len(combi[0][0])):
            if combi[j][0][i] != combi[j][1][i]:
                hammingdistance += 1
        distance.append(hammingdistance)
    return distance

def SequenceDistance(combi, k): #ratio of #difference/sequence length by using hamming distance
    ham = HammingDistance(combi)
    seqDis = []
    for i in range(len(ham)):
        seqDis.append(ham[i]/int(k))
    return seqDis

def JC_distance(d_u: float) -> float:
    """Compute Jukes-Cantor distance of FastTree's uncorrected distance

    Defined on page 1643 as d = -(3/4)log(1 - (4/3)d_u).

    Important note: Page 1643-1644
    "For both nucleotide and protein sequences, 
     FastTree truncates the corrected distances to a maximum of 3.0 substitutions per site, 
     and for sequences that do not overlap because of gaps, FastTree uses this maximum distance."
    
    Args:
        d_u (float): FastTree's uncorrected distancgie, the fraction of positions that differ between sequences.

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
