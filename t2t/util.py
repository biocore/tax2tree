#!/usr/bin/env python

from skbio.parse.sequences import parse_fasta

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def combine_alignments(fp1, fp2):
    """take two filepointers, combine the files"""
    seqs1 = dict(parse_fasta(fp1))
    seqs2 = dict(parse_fasta(fp2))

    if set(seqs1.keys()).intersection(set(seqs2.keys())):
        raise ValueError, "Conflicting sequence ids in fp1 and fp2"

    combined = seqs1
    combined.update(seqs2)

    return combined

def reroot(tree, tipnames, tmp_nodename="TEMPORARY_ROOT_NODE_NAME"):
    """Returns a tree rerooted based on tipnames"""
    node = tree.lowest_common_ancestor(tipnames)

    # make a new node that sits inbetween LCA and parent
    parent = node.parent
    parent.remove(node)
    node.parent = None
    new_node = parent.__class__()
    new_node.name = tmp_nodename

    if hasattr(new_node, 'length') and new_node.length:
        new_node.length = node.length / 2.0
        node.length = node.length / 2.0

    # add node back to tree and reconnect LCA
    parent.append(new_node)
    new_node.append(node)

    # root at the new node, unset its temporary name
    new_tree = tree.root_at(tmp_nodename)
    new_root = new_tree.find(tmp_nodename)
    new_tree.name = None
    new_root.name = None

    # collapse single descendents if they exist
    new_tree.prune()

    return new_tree

def unzip(items):
    """The inverse of zip

    Parameters
    ----------
    items : a nested iteratable

    Returns
    -------
    list
        The unzipped items

    Examples
    --------
    >>> from skbio.util.misc import unzip
    >>> unzip([[1, 2], ['a', 'b']])
    [[1, 'a'], [2, 'b']]

    """
    if items:
        return [list(i) for i in zip(*items)]
    else:
        return []
