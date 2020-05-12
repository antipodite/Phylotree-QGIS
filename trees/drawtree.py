"""
    Contains Bill Mill's implementation of the Rheingold-Tilford algorithm
    for drawing balanced trees, from https://github.com/llimllib/pymag-trees

    Will work on any tree representation that has a list of descendants
    for each node as a property of that node.

    Modifications (c) Isaac Stead 2020
"""

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
    def max_xy(self):
        all_xy = [(node.x, node.y) for node in self.walk()]
        return map(max(zip(*all_xy)))

    @property
    def min_xy(self):
        all_xy = [(node.x, node.y) for node in self.walk()]
        return map(min(zip(*all_xy)))

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

def translate(point, x, y):
    px, py = point
    return (px + x, py + y)

def slope(a, b):
    xa, ya = a
    xb, yb = b
    return (ya - yb) / (xa - xb)

def intercept(point, slope):
    x, y = point
    return y - (slope * x)

def line(a, b):
    s = slope(a, b)
    i = intercept(a, s)
    return s, i

def intersection(a, b, c, d):
    try:
        ab_s, ab_i = line(a, b)
        cd_s, cd_i = line(c, d)
        x = (cd_i - ab_i) / (ab_s - cd_s)
        y = ab_s * ((cd_i - ab_i) / (ab_s - cd_s)) + ab_i
    except ZeroDivisionError:
        return False
    return x, y

# These should really all be methods of DrawTree its getting kinda
# confusing 

def treecenter(tree):
    """
    Calculate the central point of the box formed by the maximum
    extent of the tree
    """
    coords = [node.coord for node in tree.walk()]
    xs = [coord[0] for coord in coords]
    ys = [coord[1] for coord in coords]
    max_x, min_x = max(xs), min(xs)
    max_y, min_y = max(ys), min(ys)
    return intersection(
        (min_x, min_y), (max_x, max_y),
        (min_x, max_y), (max_x, min_y)
    )

def shiftxy(tree, point):
    """
    Calculate the x and y values for translation of an entire
    tree such that the center of its bounding box now lies
    over `point`.
    """
    new_x, new_y = point
    old_x, old_y = treecenter(tree)
    return new_x - old_x, new_y - old_y

def translate_tree(tree):
    shift_x, shift_y = shiftxy(tree)
    for node in tree.walk():
        old_x, old_y = node.coord
        node.coord = (old_x - shift_x, old_y - shift_y)
        