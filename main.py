"""

Based on the FastTree algorithm as introduced by Price et al. [1]
[1] Price at al. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix. Molecular Biology and Evolution, vol 26 (7). 2009.


Authors:
    Roos Bressers
    Klarinda de Zwaan
    Bram Pronk
"""
import sys
import os
import argparse

# This is the Project Root, important that run.py is inside this folder for accurate ROOT_DIR
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

# Argument Parser
PARSER = argparse.ArgumentParser(prog='main.py', 
                                 description="Runner of Group 4's implementation of the FastTree algorithm")

PARSER.add_argument('--file', '-f',
                    dest='aln_file',
                    metavar='FILE_PATH',
                    help="path to the input .aln file",
                    default='data/test-small.aln')
PARSER.add_argument('--experiment', '-e',
                    help="Name of experiment",
                    default="experiment")

def main() -> None:
    """
    The start of the program. This function processes user input and calls the FastTree algorithm.

    @return: None
    """
    # Print all the user entered arguments in a neatly organized fashion
    args = PARSER.parse_args()

    if args.aln_file is None:
        print("No input .aln file specified. Use -f '/path/to/sequences.aln' to insert one. Exiting program.")
        sys.exit()

    with open(args.aln_file, 'r') as file:
        try:
            data = file.read().split('\n')
        except yaml.YAMLError as exc:
            print("Incorrect .aln file detected!")
            print(exc)
    
    # Process the .aln file into a dictionary mapping {name -> sequence}
    sequences = {}
    

if __name__ == "__main__":
   main()
