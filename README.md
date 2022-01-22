# FASTTREE implementation           <img src="https://www.abuse.nl/assets/logos/tudelft.png" width="100" height="50">

[![Python Version](https://img.shields.io/static/v1.svg?label=Python%20Version&message=3.8+&color=blue)](https://www.python.org/downloads)

Implementation of the FastTree algoritm by Price et. al. (Price et al. 2009). Our FastTree algorithm computes a tree given aligned sequences without gaps. Neighbor Joining, FastNJ, Top hits heuristics, Local hill climbing, and Nearest Neighbor Interchange are used. BIONJ profile weights (Gascuel et. al. 1997) and handling of duplicate sequences is implemented. 


## Installation
<!---
This section should contain installation, testing, and running instructions for people who want to get started with the project. 
- These instructions should work on a clean system.
- These instructions should work without having to install an IDE.
- You can specify that the user should have a certain operating system.
--->
This project is built solely on built-in libraries. It has been tested on Python versions 3.8 and above.
Since this project relies on the argparse library, which has been included since Python 3.2, that should be the absolute lower limit on version.

Recommended is running from a conda environment, ```conda create --name myenv python=3.8```.

Then to run the project, cd into the project's root dir, then simply call the ```main.py``` wrapper with the desired input file like so:
```bash
python main.py -f data/fasttree-input.aln
```
Simply calling 
```bash
python main.py
```
will run the program on the test-small.aln dataset, if the test-small.aln is found in the /data folder in the project's root.

To see additional arguments, use argparse's built in -h command:
```bash
python main.py -h
```

Running the program with the argument ```--verbose=1``` includes additional print statements during runtime.

Output is printed to the terminal.

<!-- TODO: run verbose -->


## General Status

This algorithm will follow the steps as discussed in the original FastTree paper (Price et al. 2009):
- Make sure all sequences are unique
- Initialize Top-Hits and FastNJ heuristics
- Generate an initial topology using Neighbor Joining and the Top-Hits heuristic (which uses FastNJ and Local Hill Climbing)
- Revise topology with Nearest Neighbor Interchange 
- Return Newick String representation of the phylogenetic tree including final branch lengths

## Implementation details

#### Neighbor Joining
Neighbor joining uses aligned sequences to make an initial tree with minimized distance. The minimization criterion is given by the profile distance - outdistance for both joined nodes. 

#### FastNJ
Each node keeps track of it's best join, which makes looking for the best join in NJ a task ran in O(m) time, since we look for m candidates.

#### Top hits heuristics
FastTree reduces the number of considered joins during Neighbor Joining. TODO

#### Local hill-climbing
Make sure the considered join is a local optimum. In other words, if there is a join with a lower NJ criterion found in the node's top-hits lists, prefer that one over the current considered join.

#### Nearest Neighbor Interchange
Refine the initial topology, by considering if alternate topologies of nodes ((A, B), (C, D)) reduce the minimum evolution criterion.

#### Lambda
Lambda is calculated using BIONJ (Gascuel et. al. 1997) and is updated every time two nodes are joined. It determines the weight of the new node by estimating the variance correction for each newly formed node.

## File Structure
![UML Diagram.png](UML%20Diagram.png)

```
.
├── LICENSE
├── UML Diagram.png
├── README.md
├── .gitignore
├── data
│       └── fasttree-input.aln
│       └── test-duplicates.aln
│       └── test-small.aln
├── main.py
└── src
    ├── fast_tree.py
    ├── heuristics.py
    ├── neighbor_joining.py
    ├── node.py
    ├── tree.py
    └── util.py
```
## Authors
    - Roos Bressers           4570405        R.Bressers@student.tudelft.nl
    - Bram Pronk              4613066        I.B.Pronk@student.tudelft.nl
    - Klarinda de Zwaan       4657136        B.K.dezwaan@student.tudelft.nl

## Citations
- Price et al. 2009. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix. Molecular Biology and Evolution, 14:2009.
- Gascuel et. al. 1997. BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. Molecular Biology and Evolution, 14:685–695.
