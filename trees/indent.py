"""
    Parse language trees in indented format.
    (c) Isaac Stead 2020

    Expects a plain text file with the tree specified like this:

    Austronesian
        Rukai
        Atayalic
        Malayo-Polynesian
            Philippines
                Tagalog
            Palauan

    At this point, indents should be 4 spaces.
"""
from collections import namedtuple
from phylo_tree.trees.newick import Node

INDENT_CHARS = ('\t', ' ')
INDENT_LEVEL = 4

def count_initial_whitespace(line):
    count = 0
    for i, char in enumerate(line):
        if char not in INDENT_CHARS:
            count = i
            break
    return count

def process_line(line, indent_size):
    ws = count_initial_whitespace(line)
    level = ws // indent_size
    name = line[ws:].strip('\n')
    return (level, name)

def load_treefile(path, indent_size=INDENT_LEVEL):    
    with open(path) as f:
        lines = f.readlines()   
    return [process_line(l, indent_size) for l in lines]

def build_tree(lines):

    # Make a copy of items so this function has no
    # bug inducing side effects from .pop()
    items = lines.copy()

    def inner(items, depth=0, level=[]):
        while items:
            d, name = items[0]
            if d == depth:
                new = Node(name)
                if not level:
                    level.append(new)
                else:
                    level.add_descendant(new)
                items.pop(0)
            elif d > depth:
                inner(items, d, new)
            elif d < depth:
                return
        return level

    return inner(items)

def read(path):
    return build_tree(load_treefile(path))
