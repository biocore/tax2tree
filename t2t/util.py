#!/usr/bin/env python
import skbio

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"


# this unrooted_copy is an EXACT copy of the code from scikit-bio
# (BSD license) with the exception that we carry edge_num if set
# we are setup to monkey patch as we cannot easily control who
# is constructing the object.
class _TreeNode(skbio.TreeNode):
    def unrooted_copy(self, parent=None):
        r"""Walks the tree unrooted-style and returns a copy

        Perform a copy of self and return a new copy of the tree as an
        unrooted copy. This is useful for defining new roots of the tree as
        the `TreeNode`.

        This method is recursive.

        Warning, this is _NOT_ a deepcopy

        Parameters
        ----------
        parent : TreeNode or None
            Used to avoid infinite loops when performing the unrooted traverse

        Returns
        -------
        TreeNode
            A new copy of the tree

        See Also
        --------
        copy
        unrooted_deepcopy
        root_at

        Examples
        --------
        >>> from skbio import TreeNode
        >>> tree = TreeNode.read(["((a,(b,c)d)e,(f,g)h)i;"])
        >>> new_tree = tree.find('d').unrooted_copy()
        >>> print(new_tree)
        (b,c,(a,((f,g)h)e)d)root;
        <BLANKLINE>

        """
        neighbors = self.neighbors(ignore=parent)
        children = [c.unrooted_copy(parent=self) for c in neighbors]

        # we might be walking UP the tree, so:
        if parent is None:
            # base edge
            edgename = None
            length = None
            edge_num = None
        elif parent.parent is self:
            # self's parent is becoming self's child
            edgename = parent.name
            length = parent.length
            edge_num = parent.edge_num
        else:
            assert parent is self.parent
            edgename = self.name
            length = self.length
            edge_num = self.edge_num

        result = self.__class__(name=edgename, children=children,
                                length=length)
        result.edge_num = edge_num

        if parent is None:
            result.name = "root"

        return result


def _convert_to_local_treenode(tree):
    for n in tree.traverse(include_self=True):
        # this is dirty
        n.__class__ = _TreeNode


def _convert_to_skbio_treenode(tree):
    for n in tree.traverse(include_self=True):
        n.__class__ = skbio.TreeNode


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
    else:
        new_node.length = 0.0

    # add node back to tree and reconnect LCA
    parent.append(new_node)
    new_node.append(node)

    # monkey patch? before we root
    _convert_to_local_treenode(tree)
    _edge_label(tree)

    # root at the new node, unset its temporary name
    new_tree = tree.root_at(tmp_nodename)
    new_root = new_tree.find(tmp_nodename)
    new_tree.name = None
    new_root.name = None

    # collapse single descendents if they exist
    new_tree.prune()

    # unmonkey patch
    _convert_to_skbio_treenode(new_tree)
    _edge_label(new_tree)

    return new_tree


def _edge_label(tree):
    # put edge labels where they don't exist use unique ones while we're
    # at it
    no_edge_label = []
    max_edge_label = -1
    for n in tree.traverse(include_self=True):
        if hasattr(n, 'edge_num') and n.edge_num is not None:
            max_edge_label = max(max_edge_label, n.edge_num)
        else:
            no_edge_label.append(n)
    for n in no_edge_label:
        max_edge_label += 1
        n.edge_num = max_edge_label


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
    >>> from t2t.util import unzip
    >>> unzip([[1, 2], ['a', 'b']])
    [[1, 'a'], [2, 'b']]

    """
    if items:
        return [list(i) for i in zip(*items)]
    else:
        return []
