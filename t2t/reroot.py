#!/usr/bin/env python

"""Reroot a tree based on tips"""

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def reroot(tree, tipnames, tmp_nodename="TEMPORARY_ROOT_NODE_NAME"):
    """Returns a tree rerooted based on tipnames"""
    node = tree.lca(tipnames)

    # make a new node that sits inbetween LCA and parent
    parent = node.parent
    parent.remove(node)
    new_node = parent.__class__(name=tmp_nodename, length=node.length / 2.0)
    node.length = node.length / 2.0

    # add node back to tree and reconnect LCA
    parent.append(new_node)
    new_node.append(node)

    # root at the new node, unset its temporary name
    new_tree = tree.root_at(tmp_nodename)
    new_root = new_tree.find(tmp_nodename)
    new_root.name = None

    # collapse single descendents if they exist
    new_tree.prune()

    return new_tree
