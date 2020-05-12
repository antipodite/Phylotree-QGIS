"""
    Tests for the plugin Tree submodule

    Isaac Stead, May 2020
"""
from os import path
from unittest import TestCase
from phylo_tree.trees import newick, indent, drawtree

TESTFILES = {
    'indentfile': path.join(path.dirname(__file__), 'test_tree.txt'),
    'newickfile': path.join(path.dirname(__file__), 'test_tree.nwk')
}

class TreeTest(TestCase):

    def setUp(self):
        tree_indent = indent.read(TESTFILES['indentfile'])
        tree_newick = newick.read(TESTFILES['newickfile'])

    def test_indent(self):
        pass

    def test_newick(self):
        pass

    def test_drawtree(self):
        drawtree.layout(tree_indent)
        self.assertFalse(intersection((-1, 0.5), (1, 0.5), (0, 1), (0, 2)))
        node_coords = [(node.x, node.y) for node in drawtree.walk()]