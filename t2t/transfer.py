from functools import reduce
from operator import or_
import numpy as np


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


def transfer(backbone, with_placement):
    """Transfer names from a backbone to another tree

    WARNING: operates inplace

    Parameters
    ----------
    backbone : skbio.TreeNode
        The backbone tree which we assume has useful names on it
    with_placement : bp.BP
        The tree derived from the backbone but which is much larger

    Raises
    ------
    KeyError
        If the backbone and placement tree appear inconsistent
    """
    index_backbone(backbone)
    mapping = indexed_to_name(backbone)
    backbone_names = frozenset({n.name for n in backbone.tips()})

    tip_with_ranknames = {}  # TODO: pull rank data from tips if/where present
    placement_backbone_positions = {}

    placement_names = np.array([with_placement.name(i)
                                for i in range(with_placement.B.size)])

    # for each node in postorder
    for i in range(with_placement.B.sum()):
        node_idx = with_placement.postorderselect(i + 1)

        # if we're a leaf, check if we have a backbone name, and index if so
        if with_placement.isleaf(node_idx):
            name = with_placement.name(node_idx)
            if name in backbone_names:
                placement_backbone_positions[node_idx] = [name, ]

                # TODO: is our tip special?
                if name in tip_with_ranknames:
                    # deepest ancestor where this and only this tip is
                    # represented, set rank names on that node
                    raise ValueError()

        # else we are not a leaf, check if we have backbone names,
        # and see if we've found the correct node to set a name
        else:
            # gather names from immediate children if any
            names = []
            child = with_placement.fchild(node_idx)
            while child != 0:
                if child in placement_backbone_positions:
                    names.extend(placement_backbone_positions[child]) # pop(child))
                child = with_placement.nsibling(child)

            if names:
                placement_backbone_positions[node_idx] = names

                # this is conservative: set the name as low as we can
                names_as_set = frozenset(names)
                if names_as_set in mapping:
                    placement_names[node_idx] = mapping.pop(names_as_set)

    if len(mapping) > 0:
        problems = '\n'.join([str(k) + ':' + str(v)
                              for k, v in mapping.items()])
        raise KeyError(f"Inconsistency detected, not transfered: "
                       f"{problems}")

    with_placement.set_names(placement_names)
    return with_placement

    # separate:
    # we need to account for species / genus names AT THE TIPS of the backbone
    # which will occur for single entry names
    # do a separate method, check for missing lower level names, place at tip


