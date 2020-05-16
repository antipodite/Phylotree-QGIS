"""
    Contains Bill Mill's implementation of the Rheingold-Tilford algorithm
    for drawing balanced trees, from https://github.com/llimllib/pymag-trees

    Will work on any tree representation that has a list of descendants
    for each node as a property of that node.

    Modifications (c) Isaac Stead 2020
"""
import os
from phylo_tree.trees.indent import read as read_indent
from phylo_tree.trees.newick import read as read_newick

class DrawTree(object):
    def __init__(self, tree, parent=None, depth=0, number=1):
        self.x = -1.
        self.y = depth
        self.tree = tree
        self.children = [DrawTree(c, self, depth+1, i+1) 
                         for i, c
                         in enumerate(tree.descendants)]
        self.parent = parent
        self.thread = None
        self.mod = 0
        self.ancestor = self
        self.change = self.shift = 0
        self._lmost_sibling = None
        #this is the number of the node in its group of siblings 1..n
        self.number = number
        # The data from self.tree
        self.name = None
        self.length = None

    def left(self): 
        return self.thread or len(self.children) and self.children[0]

    def right(self):
        return self.thread or len(self.children) and self.children[-1]

    def lbrother(self):
        n = None
        if self.parent:
            for node in self.parent.children:
                if node == self: return n
                else:            n = node
        return n

    def get_lmost_sibling(self):
        if not self._lmost_sibling and self.parent and self != \
        self.parent.children[0]:
            self._lmost_sibling = self.parent.children[0]
        return self._lmost_sibling
    lmost_sibling = property(get_lmost_sibling)

    # Isaac's mods start here

    def walk(self):
        """
        Traverse the whole tree and yield each visited node in preorder
        """
        yield self
        for node in self.children:
            for n in node.walk():
                yield n

    def all_positions(self):
        """
        Return a flat dict of self.tree: self.coords. Used to match
        calculated positions with nodes in the tree this DrawTree was 
        created from.
        """
        return {node.tree: node.coords for node in self.walk()}

    @property
    def boundingbox(self):
        """
        The central point of a box defined by the max and min coordinates
        of the tree.
        """
        # Find the corners of the tree
        all_coords = [(node.x, node.y) for node in self.walk()]
        xs, ys = zip(*all_coords)
        a = Line( (min(xs), min(ys)), (max(xs), max(ys)) )
        b = Line( (min(xs), max(ys)), (max(xs), min(ys)) )
        
        return a.intersection(b)

    def translate(self, center):
        """
        Move all subtree's coordinates such that the central point
        of the tree's bounding box is now `center`
        """
        old_center = self.boundingbox
        dx, dy = Line(old_center, center).slope_xy()
        for node in self.walk():
            node.x = node.x - dx
            node.y = node.y - dy

    @property
    def coords(self):
        return (self.x, self.y)

    def __str__(self): return "%s: x=%s mod=%s" % (self.tree, self.x, self.mod)
    def __repr__(self): return self.__str__()

def buchheim(tree):
    dt = firstwalk(DrawTree(tree))
    min = second_walk(dt)
    if min < 0:
        third_walk(dt, -min)
    return dt

def third_walk(tree, n):
    tree.x += n
    for c in tree.children:
        third_walk(c, n)

def firstwalk(v, distance=1.):
    if len(v.children) == 0:
        if v.lmost_sibling:
            v.x = v.lbrother().x + distance
        else:
            v.x = 0.
    else:
        default_ancestor = v.children[0]
        for w in v.children:
            firstwalk(w)
            default_ancestor = apportion(w, default_ancestor, distance)
        execute_shifts(v)

        midpoint = (v.children[0].x + v.children[-1].x) / 2

        ell = v.children[0]
        arr = v.children[-1]
        w = v.lbrother()
        if w:
            v.x = w.x + distance
            v.mod = v.x - midpoint
        else:
            v.x = midpoint
    return v

def apportion(v, default_ancestor, distance):
    w = v.lbrother()
    if w is not None:
        #in buchheim notation:
        #i == inner; o == outer; r == right; l == left; r = +; l = -
        vir = vor = v
        vil = w
        vol = v.lmost_sibling
        sir = sor = v.mod
        sil = vil.mod
        sol = vol.mod
        while vil.right() and vir.left():
            vil = vil.right()
            vir = vir.left()
            vol = vol.left()
            vor = vor.right()
            vor.ancestor = v
            shift = (vil.x + sil) - (vir.x + sir) + distance
            if shift > 0:
                move_subtree(ancestor(vil, v, default_ancestor), v, shift)
                sir = sir + shift
                sor = sor + shift
            sil += vil.mod
            sir += vir.mod
            sol += vol.mod
            sor += vor.mod
        if vil.right() and not vor.right():
            vor.thread = vil.right()
            vor.mod += sil - sor
        else:
            if vir.left() and not vol.left():
                vol.thread = vir.left()
                vol.mod += sir - sol
            default_ancestor = v
    return default_ancestor

def move_subtree(wl, wr, shift):
    subtrees = wr.number - wl.number
    wr.change -= shift / subtrees
    wr.shift += shift
    wl.change += shift / subtrees
    wr.x += shift
    wr.mod += shift

def execute_shifts(v):
    shift = change = 0
    for w in v.children[::-1]:
        w.x += shift
        w.mod += shift
        change += w.change
        shift += w.shift + change

def ancestor(vil, v, default_ancestor):
    #the relevant text is at the bottom of page 7 of
    #"Improving Walker's Algorithm to Run in Linear Time" by Buchheim et al, (2002)
    #http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.16.8757&rep=rep1&type=pdf
    if vil.ancestor in v.parent.children:
        return vil.ancestor
    else:
        return default_ancestor

def second_walk(v, m=0, depth=0, min=None):
    v.x += m
    v.y = depth

    if min is None or v.x < min:
        min = v.x

    for w in v.children:
        min = second_walk(w, m + v.mod, depth+1, min)

    return min

def layout(nodetree):
    """
    Calculate and add positions for python-newick Node tree
    """
    drawtree = buchheim(nodetree)
    positions = drawtree.all_positions()
    for node in nodetree.walk():
        node.coord = positions[node]

def buildtree(path):
    """The entry point into this module.

    Takes a path to a tree in a supported file format and returns
    a DrawTree object with coordinates and labels set up, ready
    for use in QGIS API.
    """
    ext = os.path.splitext(path)
    if ext == '.nwk':
        nodetree = read_newick(path)
    elif ext == '.txt':
        nodetree = read_indent(path)
    else:
        raise ValueError('Unsupported file type {}'.format(ext))

    drawtree = buchheim(nodetree)
    for node in drawtree.walk():
        node.name = node.tree.name
        node.length = node.tree.length

    return drawtree

class Line(object):
    """Quick class to represent a line for geometry operations"""
    
    def __init__(self, point_a, point_b):
        self.a = point_a
        self.b = point_b

    @property
    def slope(self):
        """The line's slope"""
        xa, ya = self.a
        xb, yb = self.b
        return (ya - yb) / (xa - xb)

    def slope_xy(self):
        """Return the rise and run of the line segment a->b"""
        xa, ya = self.a
        xb, yb = self.b
        return (xa - xb, ya - yb)

    @property
    def intercept(self):
        """The line's y-intercept"""
        x, y = self.a
        return y - (self.slope * x)

    def intersection(self, line):
        """Return the point of intersection this and `line`"""
        try:
            x = (self.intercept - line.intercept) / (line.slope - self.slope)
            y = self.slope * (x + self.intercept)
        except ZeroDivisionError:
            return False
        return (x, y)

def bestfit(points):
    """Compute best fit line for points using least square method"""
    pass
