#!/usr/bin/env python

from collections import defaultdict
from operator import itemgetter
from numpy import argmin, array, where
from skbio import TreeNode
from t2t.util import unzip
import re

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"


RANK_ORDER = ['d', 'p', 'c', 'o', 'f', 'g', 's']
BAD_NAMES = [
    'environmental sample', 'uncultured', 'UNNAMEABLE', 'unclassified',
    'unidentified', 'cluster', 'isolate', 'environmental samples']

BAD_NAMES_REGEX = re.compile("(%s)" % ')|('.join(map(str.lower, BAD_NAMES)))


def set_rank_order(order):
    """Reset the global RANK_ORDER"""
    global RANK_ORDER
    RANK_ORDER = order


def determine_rank_order(con):
    """Determines dynamically rank order based on first input con string"""
    order = [s.strip()[0] for s in con.split(';')]
    global RANK_ORDER
    RANK_ORDER = order

    return order


def has_badname(name):
    """Boolean, if name contains a badname"""
    return len(BAD_NAMES_REGEX.findall(name)) > 0


def load_consensus_map(lines, append_rank, check_bad=True,
                       check_min_inform=True, assert_nranks=True,
                       verbose=False, check_euk_unc=False):
    """Input is tab delimited mapping from tipname to a consensus string

    tipname is the tipnames in the loaded tree
    consensus string must be len(RANK_ORDER), and ';' delimited

    check_bad : check for bad names
    check_min_inform: check for informative information below domain

    If append_rank is True, rank information will be appended on. For instance,
    the name at the 0-index position of the consensus will be joined with
    RANK_ORDER[0]

    check_euk_unc : check for eukarayota or unclassified, set to none if found
    and true

    Output is a dictionary mapping tipname to consensus strings split into
    a list.
    """
    if verbose:
        print("loading consensus map...")

    mapping = {}
    n_ranks = len(RANK_ORDER)
    for line in lines:
        id_, consensus = line.strip().split('\t')
        id_ = id_.strip()

        names = [n.strip() for n in consensus.split(';')]

        if check_euk_unc and 'Eukaryota' in names[0] or \
                'Unclassified' in names[0]:
            names = [None] * n_ranks

        if assert_nranks:
            if len(names) != n_ranks:
                raise ValueError

        # clean up missing names
        for idx in range(n_ranks):
            if names[idx] == '' or names[idx] == 'None':
                names[idx] = None
            if names[idx] is not None and \
                    '__' in names[idx] and \
                    names[idx].split('__')[1] == '':
                names[idx] = None

        if check_min_inform:
            if names[1] is None:
                names = [None] * n_ranks

        # clean up bad names
        if check_bad:
            for idx, name in enumerate(names):
                if name is not None and has_badname(name.lower()):
                    names[idx] = None

        # append rank if needed
        if append_rank:
            for idx in range(n_ranks):
                if names[idx] is not None:
                    names[idx] = '__'.join([RANK_ORDER[idx], names[idx]])
                else:
                    names[idx] = "%s__" % RANK_ORDER[idx]
        mapping[id_] = names

    return mapping


def load_tree(tree, tipname_map):
    """Returns a PhyloNode tree decorated with helper attrs

    The following attributes and descriptions are decorated onto the tree:

        Consensus
            If the node is a tip, the corresponding taxonomy string is placed
            here otherwise [None] * number_of_ranks is stored

            If the node is internal, then [None] * number_of_ranks is stored

        TipStart
            The left most tip

        TipStop
            The right most tip

    Parameters
    ----------
    tree : str or TreeNode
        A newick string or a TreeNode
    tipname_map : dict
        {id_: [tax, string]}

    Returns
    -------
    TreeNode

    """
    if not isinstance(tree, TreeNode):
        tree = TreeNode.read(tree, convert_underscores=False)

    n_ranks = len(RANK_ORDER)

    missing_tax = [None] * n_ranks

    for idx, tip in enumerate(tree.tips()):
        if tip.name:
            tip.name = tip.name.replace("'", "")

        tip.TipStart = idx
        tip.TipStop = idx
        tip.Consensus = tipname_map.get(tip.name, missing_tax)

    for node in tree.postorder(include_self=True):
        if node.is_tip():
            continue

        node.TipStart = node.children[0].TipStart
        node.TipStop = node.children[-1].TipStop
        node.Consensus = missing_tax

        if node.name is None:
            node.Bootstrap = None
        else:
            try:
                node.Bootstrap = float(node.name)
                node.name = None
            except ValueError:
                node.Bootstrap = None

    return tree


def collect_names_at_ranks_counts(tree):
    """Returns total name counts for a given name at a given rank

    Assumes the Consensus attribute is present on the tips

    Parameters
    ----------
    tree : TreeNode

    Returns
    -------
    dict of dict
        Returns a 2d dict, [RANK][name] -> count

    """
    total_counts = {i: defaultdict(int) for i in range(len(RANK_ORDER))}

    for consensus in (tip.Consensus for tip in tree.tips()):
        for rank, name in enumerate(consensus):
            if name is None:
                continue
            total_counts[rank][name] += 1
    return total_counts


def decorate_name_relative_freqs(tree, total_counts, min_count):
    """Decorates relative frequency information for names on the tree

    Adds on the attribute ConsensusRelFreq which is a 2d dict containing
    the relative frequency of each name at each rank for the subtree that
    descends from a given node.

    Adds on the attribute ValidRelFreq which is a 2d dict containing
    the valid tip frequency of each name at each rank for the subtree that
    descends from a given node.

    Both these attributes will be None on tips.

    Parameters
    ----------
    tree : TreeNode
    total_counts : dict of dict
        The return data from collect_names_at_ranks_counts
    min_count : int
        is the minimum number of tips that must represent a name for that
        frequency to be retained

    """
    tips = list(tree.tips())
    for tip in tips:
        tip.ConsensusRelFreq = None

    n_ranks = len(RANK_ORDER)
    n_ranks_it = range(n_ranks)

    for n in tree.traverse(include_self=True):
        if n.is_tip():
            n.ConsensusRelFreq = None
            n.ValidRelFreq = None
            continue

        counts = {i: defaultdict(int) for i in n_ranks_it}

        # build of counts of the names at the tips per rank
        cons_at_tips = (t.Consensus for t in tips[n.TipStart:n.TipStop + 1])
        for con in cons_at_tips:
            for cur_rank, cur_name in enumerate(con):
                if cur_name is None:
                    continue
                counts[cur_rank][cur_name] += 1

        res_freq = {i: {} for i in n_ranks_it}
        res_valid = {i: {} for i in n_ranks_it}

        # collect frequency information of the names per rank
        for rank, names in counts.items():
            for name, name_counts in counts[rank].items():
                if name_counts < min_count:
                    continue

                relfreq = float(name_counts) / total_counts[rank][name]
                validfreq = float(name_counts) / n.NumTips
                res_freq[rank][name] = relfreq
                res_valid[rank][name] = validfreq

        n.ConsensusRelFreq = res_freq
        n.ValidRelFreq = res_valid


def decorate_name_counts(tree):
    """Decorates count information for names on the tree

    Adds on the attribute TaxaCount which is a 2d dict containing
    the count of each name at each rank for the subtree that
    descends from a given node.

    Parameters
    ----------
    tree : TreeNode
    total_counts : dict of dict
        The return data from collect_names_at_ranks_counts
    """

    tips = list(tree.tips())
    n_ranks = len(RANK_ORDER)
    n_ranks_it = range(n_ranks)

    for n in tree.traverse(include_self=True):
        counts = {i: defaultdict(int) for i in n_ranks_it}

        # build of counts of the names at the tips per rank
        cons_at_tips = (t.Consensus for t in tips[n.TipStart:n.TipStop + 1])
        for con in cons_at_tips:
            for cur_rank, cur_name in enumerate(con):
                if cur_name is None:
                    continue
                counts[cur_rank][cur_name] += 1

        n.TaxaCount = counts


def set_ranksafe(tree):
    """Determines what ranks are safe for a given node

    RankSafe is a len(RANK_ORDER) boolean list. True means at that rank, there
    is only a single name with > 50% relative abundance

    Parameters
    ----------
    tree : TreeNode

    """
    ranksafe = [False] * len(RANK_ORDER)
    for node in tree.traverse(include_self=True):
        node.RankSafe = ranksafe[:]

        if node.is_tip():
            continue

        for rank, names in node.ConsensusRelFreq.items():
            # this is strict
            if sum(x >= 0.5 for x in names.values()) == 1:
                node.RankSafe[rank] = True


def decorate_ntips(tree):
    """Cache the number of informative tips on the tree.

    This method will set NumTips as the number of informative tips that descend
    from a given node. If the node is a tip, and it is informative, it will
    have a count of 1. Informative is based on the presence of taxonomy
    information at a tip.

    Parameters
    ----------
    tree : TreeNode

    """
    n_ranks = len(RANK_ORDER)
    missing = [None] * n_ranks

    for node in tree.postorder(include_self=True):
        if node.is_tip():
            node.NumTips = node.Consensus != missing
        else:
            node.NumTips = sum(c.NumTips for c in node.children)


def decorate_ntips_rank(tree):
    """Cache the number of informative tips for each rank for each node.

    This method will set NumTipsRank as the number of informative tips that
    descend from a given node for each rank. Informative is based on the
    presence of taxonomy information at a tip for a give rank.

    Parameters
    ----------
    tree : TreeNode

    """
    n_ranks = len(RANK_ORDER)

    for node in tree.postorder(include_self=True):
        counts = defaultdict(int)
        for r in range(n_ranks):
            if node.is_tip():
                counts[r] = node.Consensus[r] is not None
            else:
                counts[r] = sum(c.NumTipsRank[r] for c in node.children)

        node.NumTipsRank = counts


def pick_names(tree):
    """Does an initial decoration of names on the tree

    The best name by relative frequency is placed, from kingdom -> species,
    up until the first uninformative name following an informative name is
    placed

    """
    names_prealloc = [None] * len(RANK_ORDER)

    for node in tree.non_tips(include_self=True):
        names = names_prealloc[:]
        count = 0

        # set names at ranksafe nodes, stop if we've set a name and descendent
        # rank names are not safe.
        for rank, is_safe in enumerate(node.RankSafe):
            if is_safe:
                # place best name
                count += 1
                relfreq = node.ConsensusRelFreq[rank].items()
                names[rank] = sorted(relfreq, key=itemgetter(1))[-1][0]
            else:
                # if we've had one or more useless rank, set remaining to None
                if count >= 1:
                    break

        node.RankNames = names


def fmeasure(precision, recall):
    """Returns the fmeasure (or F1-score)

    http://en.wikipedia.org/wiki/F1_score
    """
    return 2.0 * ((precision * recall) / (precision + recall))


def fpoint5measure(precision, recall):
    """Returns the f0.5measure (or F0.5-score)"""
    beta = 0.5
    betasqrd = beta ** 2

    tmp = (precision * recall) / ((betasqrd * precision) + recall)
    return (1 + betasqrd) * tmp


def f2measure(precision, recall):
    """Returns the f2measure (or F2-score)"""
    return (1.0 + (2 ** 2)) * ((precision * recall) /
                               ((2 ** 2 * precision) + recall))


def min_tips(nodes):
    """For a list of nodes, return the node with the fewest tips

    implemented naively...
    """
    scores = []
    for n in nodes:
        if n is None:
            scores.append(99999999999)
        else:
            scores.append(len(list(n.tips())))
    return nodes[argmin(scores)]


def name_node_score_fold(tree, score_f=fmeasure, tiebreak_f=min_tips,
                         verbose=False):
    """Compute name scores for internal nodes, pick the 'best'

    For this method, we traverse the tree once building up a dict of scores
    for names and nodes, we can then pick the 'best' node out of the dict
    to avoid horrible lookups in the tree
    """

    if verbose:
        print("Starting name_node_score_fold...")

    name_node_score = {i: {} for i in range(len(RANK_ORDER))}
    n_ranks = len(RANK_ORDER)

    for node in tree.non_tips(include_self=True):
        node.RankNameScores = [None] * n_ranks

        for rank, name in enumerate(node.RankNames):
            if name is None:
                continue

            # precision in this case is the percent of informative tips that
            # descend that are of the name relative to the number of
            # informative tips that descend
            precision = node.ValidRelFreq[rank][name]

            # recall in this case is the percent of informative tips that
            # descent that are of the name relative to the total number of
            # tips in the tree with name
            recall = node.ConsensusRelFreq[rank][name]

            # calculate score and save it for the corrisponding rank position
            # so that these values can be examined later in other contexts
            score = score_f(precision, recall)
            node.RankNameScores[rank] = score

            if name not in name_node_score[rank]:
                name_node_score[rank][name] = []
            name_node_score[rank][name].append((node, score))

    # run through the built up dict and pick the best node for a name
    used_scores = {}
    for rank, names in name_node_score.items():
        used_scores[rank] = []

        for name, node_scores in names.items():
            node_scores_sorted = sorted(node_scores, key=itemgetter(1))[::-1]
            nodes, scores = unzip(node_scores_sorted)
            scores = array(scores)

            used_scores[rank].append((name, scores[0]))

            # if there is a tie in scores...
            if sum(scores == scores[0]) > 1:
                # ugly hack to get around weird shape mismatch
                indices = where(scores == scores[0], range(len(nodes)), None)
                tie_nodes = []
                for i in indices:
                    if i is not None:
                        tie_nodes.append(nodes[i])
                    else:
                        tie_nodes.append(None)
                node_to_keep = tiebreak_f(tie_nodes)
                for node, score in node_scores_sorted:
                    if node == node_to_keep:
                        continue
                    else:
                        node.RankNames[rank] = None
            else:
                for node, score in node_scores_sorted[1:]:
                    node.RankNames[rank] = None

    return used_scores


def score_tree(tree, verbose=False):
    """Scores the tree based on RankNameScores and tip coverage

    if the node has a name in RankNames, multiply score in RankNameScores @ idx
    by the number of tips that are "valid".. Sum these scores, divide sum by
    total number of tips covered (a tip can be counted multiple times)
    """
    total_score = 0.0
    tip_count = 0

    if verbose:
        print("Scoring tree...")

    for n in tree.non_tips(include_self=True):
        for idx, name in enumerate(n.RankNames):
            if name is None:
                continue
            score = n.RankNameScores[idx]
            total_score += score * n.NumTips
            tip_count += n.NumTips
    return total_score / tip_count


def set_preliminary_name_and_rank(tree):
    """Sets names and rank at a node

    This method is destructive: will destroy the Name attribute on tree
    """

    n_ranks = len(RANK_ORDER)
    empty_ranknames = [None] * n_ranks

    for node in tree.non_tips(include_self=True):
        node.name = None
        node.Rank = None

        if empty_ranknames == node.RankNames:
            continue

        # take the "deepest" name
        for idx, name in enumerate(node.RankNames[::-1]):
            if name is None:
                continue
            node.name = name
            node.Rank = n_ranks - (idx + 1)  # adjust for 0-based index
            break


def make_consensus_tree(cons_split, check_for_rank=True, tips=None):
    """Returns a mapping by rank for names to their parent names and counts"""

    god_node = TreeNode(name=None)
    god_node.Rank = None

    base = list(cons_split)[0]
    #  SMJ base = next(iter(cons_split))
    cur_node = god_node

    # create a base path in the tree
    for rank, name in enumerate(base):
        new_node = TreeNode(name=name)
        new_node.Rank = rank
        cur_node.append(new_node)
        cur_node = new_node

    # setup the initial childlookup structure so taht we don't have to
    # always iterate over .children
    for n in god_node.traverse(include_self=True):
        if n.is_tip():
            n.ChildLookup = {}
            continue
        n.ChildLookup = {n.children[0].name: n.children[0]}

    # for every consensus string, start at the "god" node
    for idx, con in enumerate(cons_split):
        cur_node = god_node

        # for each name, see if we've seen it, if not, add that puppy on
        for rank, name in enumerate(con):
            if name in cur_node.ChildLookup:
                cur_node = cur_node.ChildLookup[name]
            else:
                new_node = TreeNode(name=name)
                new_node.Rank = rank
                new_node.ChildLookup = {}
                cur_node.append(new_node)
                cur_node.ChildLookup[name] = new_node
                cur_node = new_node

        if tips is not None:
            cur_node.append(TreeNode(name=tips[idx]))

    # build an assist lookup dict
    lookup = {}
    for node in god_node.traverse():
        if node.name is None:
            continue
        if check_for_rank and '__' in node.name and \
                node.name.split('__')[1] == '':
            continue
        lookup[node.name] = node

    return god_node, lookup


def get_nearest_named_ancestor(node):
    """Returns node of nearest .name'd ancestor

    returns None if there does not exist a named ancestor
    """
    curr = node.parent
    while curr is not None and curr.name is None:
        curr = curr.parent
    return curr


def backfill_names_gap(tree, consensus_lookup, verbose=False):
    """Fill in missing names

    use the consensus tree mapped by nodes_lookup to attempt to fill in missing
    consensus names. For instance, if we know, by the the consensus tree,
    that foo has the parent bar, then if we see the name foo in the real tree
    and it does not have a parent named bar, we can add that name on.

    We set the attribute BackFillNames here as we want to attempt to collapse
    names later if we sanely and easily can
    """
    for node in tree.non_tips(include_self=True):

        if node.name is not None:
            node.BackFillNames = [node.name]
        else:
            node.BackFillNames = []
            continue

        if node.parent is None:
            continue
        if node.Rank is None:
            continue
        if node.is_tip():
            continue
        # if node is kingdom, more on
        if node.Rank == 0:
            continue

        # find nearest ranked parent
        named_ancestor = get_nearest_named_ancestor(node)

        if named_ancestor is None:
            if verbose:
                print("Unable to find a named parent for %s" % (node.name))
            continue

        levels = node.Rank - named_ancestor.Rank

        # this is indicative of a problem with the consensus strings
        if levels < 0:
            print(node.name, node.Rank, named_ancestor.Rank)
            print('\t', node.RankSafe)
            print('\t', node.RankNames)
            print('\t', named_ancestor.RankSafe)
            print('\t', named_ancestor.RankNames)
            node.BackFillNames = []
            continue
        elif levels == 1:
            continue

        # walk consensus tree for missing names
        names = walk_consensus_tree(consensus_lookup, node.name, levels)
        node.BackFillNames = names


def walk_consensus_tree(lookup, name, levels, reverse=True, verbose=False):
    """Walk up the consensus tree for n levels, return names

    if reverse is True, names are [::-1]
    """
    node = lookup[name]
    names = [name]
    curr = node.parent

    for _ in range(1, levels):
        if curr.Rank is None:
            # at root...
            break

        if curr is None:
            if verbose:
                print("Possible malformed consensus tree! See node %s" % name)
            continue
        if curr.name is None:
            if verbose:
                print("Gap in consensus tree! See node %s" % name)
            names.append('%s__' % RANK_ORDER[curr.Rank])
        else:
            names.append(curr.name)
        curr = curr.parent

    if reverse:
        names = names[::-1]

    return names


def commonname_promotion(tree):
    """Promote names if possible from BackFillNames"""
    for node in tree.preorder(include_self=True):
        queue = [c for c in node.children[:] if not c.is_tip()]
        backfill_nodes = []
        push_name = True

        # collect nearest descendants backfilled names
        while queue:
            curr = queue.pop()
            if len(curr.BackFillNames) == 0:
                queue.extend([c for c in curr.children[:] if not c.is_tip()])
            elif len(curr.BackFillNames) == 1:
                push_name = False
                break
            else:
                backfill_nodes.append(curr)

        # see if there is 100% overlap in a name at a given rank, if so
        # we can and should put it on the deeper node
        while push_name:
            cur_name = [n.BackFillNames[0] for n in backfill_nodes]
            if len(set(cur_name)) == 1:
                node.BackFillNames.append(cur_name[0])

                for n in backfill_nodes:
                    n.BackFillNames.pop(0)
                    if len(n.BackFillNames) <= 1:
                        push_name = False
            else:
                push_name = False

    # set the .name attribute on the tree based on .BackFillNames
    for node in tree.preorder(include_self=True):
        if node.is_tip():
            continue

        if node.BackFillNames:
            unique = []
            for name in node.BackFillNames:
                unique.append(TaxaName.getTaxaName(name))
            node.BackFillNames = unique
        if len(node.BackFillNames) > 1:
            node.name = '; '.join(node.BackFillNames)
        elif len(node.BackFillNames) == 1:
            node.name = TaxaName.getTaxaName(node.BackFillNames[0])
        else:
            node.name = None


class TaxaName(object):

    """Support object to provide unique taxa names"""
    _names = {}

    def __init__(self):
        raise NotImplementedError("Object not meant for instantiation")

    @classmethod
    def getTaxaName(cls, request):
        return request


def make_names_unique(tree, append_suffix=True, suffix_glue_char='_',
                      verbose=False):
    """Appends on a unique number if multiple of the same names exist

    ordered by number of tips, ie, _1 has more tips that _2

    expects .BackFillNames to be set
    """
    if verbose:
        print("Assigning unique tags to duplicate names...")

    # build up a dict of the names and how many tips descend
    name_lookup = {}
    for node in tree.non_tips(include_self=True):
        if node.name is None:
            continue
        else:
            for idx, name in enumerate(node.BackFillNames):
                if name not in name_lookup:
                    name_lookup[name] = []
                name_info = ((node.TipStop - node.TipStart), idx, node)
                name_lookup[name].append(name_info)

    # assign unique numbers based on the number of tips that descend
    for name, scores_and_nodes in name_lookup.items():
        sorted_scores = sorted(scores_and_nodes, key=lambda x: x[0])[::-1]
        for count, (score, idx, node) in enumerate(sorted_scores):
            # only assign a number of we have more than 1
            if count > 0:
                if node.BackFillNames[idx].split('__')[1] != '':
                    if append_suffix:
                        unique_name = suffix_glue_char.join(
                            [node.BackFillNames[idx], str(count)])
                        node.BackFillNames[idx] = unique_name

    # should probably be refactored, but assign .name based on .BackFillNames
    for node in tree.non_tips(include_self=True):
        if len(node.BackFillNames) == 0:
            node.name = None
        elif len(node.BackFillNames) == 1:
            node.name = node.BackFillNames[0]
        else:
            node.name = '; '.join(node.BackFillNames)


def pull_consensus_strings(tree, verbose=False, append_prefix=True):
    """Pulls consensus strings off of tree

    assumes .name is set
    """
    if verbose:
        print("Pulling consensus strings...")

    constrings = []
    rank_order_rev = {r: i for i, r in enumerate(RANK_ORDER)}
    # start at the tip and travel up
    for tip in tree.tips():
        if append_prefix:
            consensus_string = ['%s__' % r for r in RANK_ORDER]
        else:
            consensus_string = ['' for r in RANK_ORDER]

        tipid = tip.name
        n = tip.parent

        # walk up the tree filling in the consensus string
        while n.parent:
            if n.name:
                if ';' in n.name:
                    names = [r.strip() for r in n.name.split(';')]
                    for name in names:
                        rank_idx = rank_order_rev[name[0]]
                        consensus_string[rank_idx] = name
                else:
                    rank_idx = rank_order_rev[n.name[0]]
                    consensus_string[rank_idx] = n.name
            n = n.parent

        # if there is a name at the root we need to make sure we grab it
        if n.name:
            if ';' in n.name:
                names = [r.strip() for r in n.name.split(';')]
                for name in names:
                    rank_idx = rank_order_rev[name[0]]
                    consensus_string[rank_idx] = name
            else:
                rank_idx = rank_order_rev[n.name[0]]
                consensus_string[rank_idx] = n.name

        # join strings with tip id
        constrings.append('\t'.join([tipid, '; '.join(consensus_string)]))
    return constrings


def save_bootstraps(tree, verbose=False):
    """Retains .Bootstrap if set in .name"""
    if verbose:
        print("Attempting to retain bootstrap values")
    for n in tree.non_tips(include_self=True):
        if n.Bootstrap is not None:
            if n.name is None:
                n.name = str(n.Bootstrap)
            else:
                if ';' in n.name:
                    n.name = ':'.join([str(n.Bootstrap), "%s" % n.name])
                else:
                    n.name = ':'.join([str(n.Bootstrap), n.name])


def getpath(foo):
    while foo.parent:
        print(foo.name)
        if hasattr(foo, 'RankNames'):
            print(foo.RankNames)
        foo = foo.parent
    print(foo.name)
    if hasattr(foo, 'RankNames'):
        print(foo.RankNames)


def is_float(s):
    """Returns True if the value in string s is a float"""
    if s is None:
        return False

    if '.' not in s:
        return False

    try:
        float(s)
        return True
    except ValueError:
        return False


def validate_all_paths(tree):
    """Walk each path in the tree and make sure there aren't any conflicts"""
    # helper getpath method
    def getpath_f(n):
        path = []
        n = n.parent
        while n.parent:
            if n.name is not None:
                if ':' in n.name:
                    names = n.name.split(':', 1)
                    if not is_float(names[1]):
                        path.append(names[1])
                else:
                    if not is_float(n.name):
                        path.append(n.name)
            n = n.parent
        if n.name is not None:
            path.append(n.name)

        clean_path = []
        for p in path:
            if ';' in p:
                clean_path.extend(p.split(';')[::-1])
            else:
                clean_path.append(p)
        return clean_path

    rank_order_rev = {r: i for i, r in enumerate(RANK_ORDER)}

    bad_tips = []

    for tip in tree.tips():
        tip_path = getpath_f(tip)

        # look up the ranks on the names
        vals = [rank_order_rev[p[0]] for p in tip_path]

        # sorting shouldn't change order, if it does, ranks are out of order
        if sorted(vals)[::-1] != vals:
            bad_tips.append(tip)

        # casting to a set shouldn't change the number of ranks present,
        # if it does we have a duplicate rank
        elif len(vals) != len(set(vals)):
            bad_tips.append(tip)

    return bad_tips
