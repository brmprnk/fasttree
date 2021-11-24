"""
File that holds all unit tests for the project.
If the size becomes too large, it is recommended to split into a tests folder.

This file supports both 'nose' and running the tests via unit test.

To run:
With nose:
    $project-root: pip install nose
    $project-root: nosetests

Normally:
    $project-root: python src/test_unit.py
"""

import unittest
# import src.fast_tree.JC_distance
from src.fast_tree import JC_distance


class TestJukesCantor(unittest.TestCase):

    def test_jc1(self):
        """Simple test of positive d_u, answer found by online calculator"""
        self.assertEqual(round(JC_distance(5/31), 4), 0.1816)


if __name__ == '__main__':
    unittest.main()