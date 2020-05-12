"""
    Tests for the plugin Tree submodule

    Isaac Stead, May 2020
"""
import unittest
from os import path
from phylo_tree.trees import newick, indent, drawtree

TESTFILES = {
    'indentfile': path.join(path.dirname(__file__), 'test_tree.txt'),
    'newickfile': path.join(path.dirname(__file__), 'test_tree.nwk')
}

class TreeTest(unittest.TestCase):

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
        try:
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
            segments = [((node.x, node.y), (node.parent.x, node.parent.y))
                        for node in tree_indent.walk()]
            lines = LineCollection(segments, linewidth=1, linestyle='solid')
            figure, axis = plt.subplots()
            axis.add_collection(lines)
            plt.show()
        except ModuleNotFoundError:
            print('matplotlib not installed, skipping...')

if __name__ == '__main__':
    unittest.main()
