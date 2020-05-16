"""
    Tests for the plugin Tree submodule

    Isaac Stead, May 2020
"""
import unittest
from os import path
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

TESTFILES = {
    'indentfile': path.join(path.dirname(__file__), 'test_tree.txt'),
    'newickfile': path.join(path.dirname(__file__), 'test_tree.nwk')
}

class TreeTest(unittest.TestCase):

    def setUp(self):
        from phylo_tree.trees import newick, indent, drawtree
        tree_indent = indent.read(TESTFILES['indentfile'])
        tree_newick = newick.read(TESTFILES['newickfile'])

    def test_indent(self):
        pass

    def test_newick(self):
        pass

    def test_drawtree(self):
        drawtree.layout(tree_indent)
        self.assertFalse(intersection((-1, 0.5), (1, 0.5), (0, 1), (0, 2)))
        segments = [((node.x, node.y), (node.parent.x, node.parent.y))
                    for node in tree_indent.walk()]
        lines = LineCollection(segments, linewidth=1, linestyle='solid')
        figure, axis = plt.subplots()
        axis.add_collection(lines)
        plt.show()

