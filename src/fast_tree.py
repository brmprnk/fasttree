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

    # Calculate Jukes-Cantor distance (d in paper)
    # What log should we use?
    jd_d = -0.75 * math.log(1 - (4/3) * d_u, math.e)

    # Distances are truncated to 3 substitutions per site
    # Does that mean if d_u >= 0.75 then jd_d = 3.0?

    # For sequences that do not overlap, FastTree uses a max distance of 3.0
    if jd_d > max_distance:
        return max_distance
    
    return jd_d
