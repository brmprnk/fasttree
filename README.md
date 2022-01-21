# FASTTREE implementation           <img src="https://www.abuse.nl/assets/logos/tudelft.png" width="100" height="50">

[![Python Version](https://img.shields.io/static/v1.svg?label=Python%20Version&message=3.10&color=blue)](https://www.python.org/downloads)

Implementation of the FastTree algoritm by Price et. al. (Price et al. 2009). Our FastTree algorithm computes a tree given aligned sequences without gaps. Neighbor Joining, FastNJ, Top hits heuristics, Local hill climbing, and Nearest Neighbor Interchange are used. BIONJ profile weights (Gascuel et. al. 1997) and handling of duplicate sequences is implemented. 


## Installation
<!---
This section should contain installation, testing, and running instructions for people who want to get started with the project. 
- These instructions should work on a clean system.
- These instructions should work without having to install an IDE.
- You can specify that the user should have a certain operating system.
--->
TODO help ik snap het niet
Recommended installation uses a new Anaconda environment. To ease the process, this project includes an environment file.
This can be plugged into Anaconda [following this short tutorial](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).

Then to run the project, simply call the ```run.py``` wrapper with the desired config file like so:
```bash
python run.py -c configs/geme.yaml
```
This default call to run.py will run all implemented models. For specific models, use one or combine the 
```-poe -moe -mofa -mvib -cgae``` flags.
To see additional arguments, use argparse's built in -h command:
```bash
python run.py -h
```


- The experiment name is set in the config file, or in the command line using ```-e experiment_name```
- All results and model outputs will be stored in the project's ```results``` directory
- All metrics are written to a TensorBoard file, also found in the results folder
- All informational print statements are saved to a log.txt file, using a logger found in the util folder.



<!-- TODO: run verbose -->
## To-do's
TODO: moeten we dit doen?
To see some open issues, have a look at this Drive doc
## General Status

This algorithm will follow the steps as discussed in the original FastTree paper (Price et al. 2009):
- Make sure all sequences are unique
- Initialize Top-Hits and FastNJ heuristics
- Generate an initial topology using Neighbor Joining
- Revise topology with Nearest Neighbor Interchange 
- Return Newick String representation of the phylogenetic tree including final branch lengths

## Implementation details

#### Neighbor Joining
Neighbor joining uses aligned sequences to make an initial tree with minimized distance. The minimization criterion is given by the profile distance - outdistance for both joined nodes. 

#### FastNJ
TODO
#### Top hits heuristics
FastTree reduces the number of considered joins during Neighbor Joining. TODO

#### Local hill-climbing
TODO

#### Nearest Neighbor Interchange
TODO

#### Lambda
Lambda is calculated using BIONJ (Gascuel et. al. 1997) and is updated every time two nodes are joined. It determines the weight of the new node by estimating the variance correction for each newly formed node.

## File Structure
TODO UML in document
```
.
├── LICENSE
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
    - Roos Bressers                   R.Bressers@student.tudelft.nl
    - Bram Pronk                      I.B.Pronk@student.tudelft.nl
    - Klarinda de Zwaan               B.K.dezwaan@student.tudelft.nl

## Citations
- Price et al. 2009. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix. Molecular Biology and Evolution, 14:2009.
- Gascuel et. al. 1997. BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. Molecular Biology and Evolution, 14:685–695.
