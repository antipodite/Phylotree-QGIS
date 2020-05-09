INDENT_CHARS = ('\t', ' ')
TEST_PATH = 'test/test_subgrouping.txt'

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

def load_treefile(path, indent_size=4):    
    with open(path) as f:
        lines = f.readlines()   
    return [process_line(l, indent_size) for l in lines]

def build_tree(items, depth=0, level=None):
    if level is None:
        level = []
    while items:
        d, name = items[0]
        if d == depth:
            children = []
            level.append((name, children))
            items.pop(0)
        elif d > depth:
            build_tree(items, d, children)
        elif d < depth:
            return
    return level[0] if level else None

def build_tree_from_file(path):
    return build_tree(load_treefile(path))