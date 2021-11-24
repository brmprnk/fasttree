"""
Maybe some comments about implementation here.
"""
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
    # calculate Profile for internal nodes
    Hypothetical_internal_node = [sequences['2'] , sequences['4'], sequences['5']]
    Node_profile = Profile(Hypothetical_internal_node)
    Sequence_distance = uncorrectedDistance(Node_profile)

    combi = makeCombisofChildren(Hypothetical_internal_node)
    Sequence_distance2 = SequenceDistance(combi, len(sequences['2']))  #calculates list of all combi's of children
    print('Profile between 2 and 4 is: ' + str(Node_profile))
    print('Their sequence distance (uncorrected distance) is: ' + str(Sequence_distance))
    print('Their sequence distance (uncorrected distance) is: ' + str(Sequence_distance2)) 


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