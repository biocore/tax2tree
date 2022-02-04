from functools import reduce
from operator import or_


def index_backbone(tree):
    """Construct covered tip attributes for named internal nodes

    WARNING: operates inplace

    Parameters
    ----------
    tree : skbio.TreeNode
        The tree to operate on, assumes some nodes are named
    """
    # this feels more awkward to do in a single tree traversal,
    # so let's just do two and call it good.

    # pass 1: post order reduction of tipnames
    for node in tree.postorder(include_self=True):
        if node.is_tip():
            node.covered_tips = frozenset([node.name, ])
        else:
            node.covered_tips = reduce(or_, [c.covered_tips
                                             for c in node.children])

    # pass 2: remove cover from unnamed and tips
    for node in tree.traverse(include_self=True):
        if node.is_tip():
            node.covered_tips = frozenset()
        else:
            if node.name is None:
                node.covered_tips = frozenset()


def indexed_to_name(tree):
    """Pull the .covered_tips attribute mapping to node name

    Parameters
    ----------
    tree : skbio.TreeNode
        The tree to operate on, assumes some nodes have the .covered_tips
        attribute

    Raises
    ------
    KeyError
        If a nonunique mapping will be created
    ValueError
        If a node with covered_tips lacks a name

    Returns
    -------
    dict
        A mapping of frozenset : node name.
    """
    mapping = {}
    for node in tree.non_tips(include_self=True):
        if node.covered_tips:
            if node.name is None:
                raise ValueError("node is lacking a name")
            if node.covered_tips in mapping:
                raise KeyError("cannot construct unique mapping")
            mapping[node.covered_tips] = node.name
    return mapping
