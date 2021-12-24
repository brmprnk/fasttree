"""


References to page numbers in this code are referring to the paper or its supplementary material:
[1] Price at al. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix.
    Molecular Biology and Evolution, vol 26 (7). 2009.
    Paper can be found on https://pubmed.ncbi.nlm.nih.gov/19377059/
"""

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
        self.tophits = []

        self.m = m
        self.age = 0

        # A top-hit list is too short if it has less than 0.8m entries,
        # where 0.8 is a tuning parameter and m is the size of the top-hit lists (see Supplement.)
        self.refreshFactor = 0.8