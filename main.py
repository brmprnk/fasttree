"""

Based on the FastTree algorithm as introduced by Price et al. [1]
[1] Price at al. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix.
    Molecular Biology and Evolution, vol 26 (7). 2009.


Authors:
    Roos Bressers
    Klarinda de Zwaan
    Bram Pronk
"""
import sys
import os
import argparse

from src.fast_tree import fast_tree


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
    """Program entry function.

    This function processes user input and calls the FastTree algorithm.

    Returns:
        (None)
    """
    # Get arguments from parser
    args = PARSER.parse_args()

    # Check for input file in args
    if args.aln_file is None:
        print("No input .aln file specified. Use -f '/path/to/sequences.aln' to insert one. Exiting program.")
        sys.exit()

    # Fetch sequences from input file
    sequences = readinput(args.aln_file)

    # Run FastTree
    fast_tree(sequences)


def readinput(filepath):
    """Processes .aln file into a mapping.

    Reads in a file whose path was provided in the program's arguments.
    Assumes the input is .aln formatted, then loops through the file
    and creates mappings of the sequences to their names.

    Args:
        filepath (str): String containg the path of the input file

    Returns:
        (dict): A mapping of sequences to their names, as found in the input file
    """
    # Check if .aln format
    if not filepath.endswith('.aln'):
        print("Input file should be .aln formatted!")
        sys.exit()

    # Try reading in the file
    try:
        file = open(filepath, 'r')
    except IOError as exception:
        print("Input file could not be opened : ")
        print(exception)
        sys.exit()

    # File has successfully been loaded in. Now process data from it
    with open(filepath, 'r') as file:
        data = file.read().splitlines()  # Read in file, get rid of '\n' characters

    # Process the .aln file into a dictionary mapping {name -> sequence}
    sequences = {}
    for i in range(0, len(data) - 1, 2):  # Since file is in format (name, sequence pairs), move in steps of 2
        name = data[i][1:]  # [1:] to get rid of '>' character
        sequence = data[i + 1]

        sequences[name] = sequence

    return sequences

if __name__ == "__main__":
    main()
