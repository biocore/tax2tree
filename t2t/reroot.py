#!/usr/bin/env python

"""Reroot a tree based on tips"""

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def reroot(tree, tipnames, tmp_nodename="TEMPORARY_ROOT_NODE_NAME"):
    """Returns a tree rerooted based on tipnames"""
    node = tree.lowestCommonAncestor(tipnames)

    # make a new node that sits inbetween LCA and parent
    parent = node.Parent
    parent.removeNode(node)
    node.Parent = None
    new_node = parent.__class__()
    new_node.Name = tmp_nodename
    new_node.Length = node.Length / 2.0
    node.Length = node.Length / 2.0

    # add node back to tree and reconnect LCA
    parent.append(new_node)
    new_node.append(node)

    # root at the new node, unset its temporary name
    new_tree = tree.rootedAt(tmp_nodename)
    new_root = new_tree.getNodeMatchingName(tmp_nodename)
    new_root.Name = None
    
    # collapse single descendents if they exist
    new_tree.prune()

    return new_tree
    

