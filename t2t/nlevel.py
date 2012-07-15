#!/usr/bin/env python

from optparse import OptionParser, make_option
from cogent.parse.tree import DndParser
from string import lower
from operator import itemgetter,add
from cogent.core.tree import TreeNode
from numpy import argmin, array, where
from cogent.util.misc import unzip
import re

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

options = [make_option('--consensus-map', dest='consensus_map', \
                       type='string'),
           make_option('--append-rank', dest='append_rank', \
                       action='store_true', default=False),
           make_option('--tree', dest='tree', type='string'),
           make_option('--output', dest='output', type='string'),
           make_option('--verbose', dest='verbose', action='store_true',\
                       default=False)]

RANK_ORDER = ['k', 'p', 'c', 'o', 'f', 'g', 's']
BAD_NAMES = ['environmental sample', 'uncultured', 'UNNAMEABLE', 'unclassified',
             'unidentified','cluster','isolate', 'environmental samples']

BAD_NAMES_REGEX = re.compile("(%s)" % ')|('.join(map(lower, BAD_NAMES)))

def set_rank_order(order):
    """Reset the global RANK_ORDER"""
    global RANK_ORDER
    RANK_ORDER = order

def determine_rank_order(con):
    """Determines dynamically rank order based on first input con string"""
    order = [s[0] for s in con.split('; ')]
    global RANK_ORDER
    RANK_ORDER = order

def has_badname(name):
    """Boolean, if name contains a badname"""
    return len(BAD_NAMES_REGEX.findall(name)) > 0

def load_consensus_map(lines, append_rank, check_bad=True, check_min_inform=True, assert_nranks=True, verbose=False,
        check_euk_unc=False):
    """Input is tab delimited mapping from tipname to a consensus string

    tipname is the tipnames in the loaded tree
    consensus string must be len(RANK_ORDER), and '; ' delimited
    
    check_bad : check for bad names
    check_min_inform: check for informative information below domain

    If append_rank is True, rank information will be appended on. For instance,
    the name at the 0-index position of the consensus will be joined with 
    RANK_ORDER[0] 

    check_euk_unc : check for eukarayota or unclassified, set to none if found and true

    Output is a dictionary mapping tipname to consensus strings split into
    a list.
    """
    if verbose:
        print "loading consensus map..."

    mapping = {}
    n_ranks = len(RANK_ORDER)
    for line in lines:
        id_, consensus = line.strip().split('\t')

        names = consensus.split('; ')

        if check_euk_unc and 'Eukaryota' in names[0] or 'Unclassified' in names[0]:
            names = [None] * n_ranks

        if assert_nranks:
            assert len(names) == n_ranks

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

def load_tree(input, tipname_map, verbose=False):
    """Returns a PhyloNode tree decorated with helper attrs
    
    Helper attrs include Consensus, TipStart and TipStop. Nontips and tips that
    do not have consensus information will have [None] * len(RANK_ORDER) set 
    as Consensus
    """
    if verbose:
        print "loading tree..."
    if isinstance(input, TreeNode):
        tree = input
    else:
        tree = DndParser(input)

    tips = tree.tips()
    n_ranks = len(RANK_ORDER)

    for idx, tip in enumerate(tips):
        tip.TipStart = idx
        tip.TipStop = idx
        tip.Consensus = tipname_map.get(tip.Name, [None] * 7)

        if verbose and tip.Consensus is None:
            print "No consensus for %s" % tip.Name

    for node in tree.postorder(include_self=True):
        if node.istip():
            continue
        node.TipStart = node.Children[0].TipStart
        node.TipStop = node.Children[-1].TipStop
        node.Consensus = [None] * n_ranks

        if node.Name is None:
            node.Bootstrap = None
        else:
            try:
                node.Bootstrap = float(node.Name)
                node.Name = None
            except:
                if verbose:
                    print "Could not save bootstrap %s, node is root: %s" % \
                                       (node.Name, str(node.Parent == None))
                node.Bootstrap = None

    for tip in tree.tips():
        if tip.Name:
            tip.Name = tip.Name.replace("'","")
    return tree

def collect_names_at_ranks_counts(tree, verbose=False):
    """Returns total name counts for a given name at a given rank
    
    Assumes the Consensus attribute is present on the tips
    
    Returns a 2d dict, [RANK][name] -> count
    """
    if verbose:
        print "collecting total counts..."
    list_of_con = [tip.Consensus for tip in tree.tips()]
    total_counts = dict([(i,{}) for i in range(len(RANK_ORDER))])

    for consensus in list_of_con:
        for rank, name in enumerate(consensus):
            if not name:
                continue
            if name not in total_counts[rank]:
                total_counts[rank][name] = 0
            total_counts[rank][name] += 1
    return total_counts

def decorate_name_relative_freqs(tree, total_counts, min_count, verbose=False):
    """Decorates ConsensusRelFreq and ValidRelFreq on tree

    Adds on the attribute ConsensusRelFreq which is a 2d dict containing
    the relative frequency of each name at each rank

    Adds on the attribute ValidRelFreq which is a 2d dict containing
    the valid tip frequency of each name at each rank

    min_count is the minimum number of tips that must represent a name for that
    frequency to be retained

    Tips will have attr as None
    """
    if verbose:
        print "HYBRID! decorating relative frequencies..."
    tips = tree.tips()
    for tip in tips:
        tip.ConsensusRelFreq = None

    n_ranks = len(RANK_ORDER)

    for n in tree.nontips(include_self=True):
        counts = dict([(i, {}) for i in range(n_ranks)])
       
        # build of counts of the names at the tips per rank
        cons_at_tips = [tip.Consensus for tip in tips[n.TipStart:n.TipStop+1]]
        for con in cons_at_tips:
            for cur_rank, cur_name in enumerate(con):
                if cur_name is None:
                    continue
                if cur_name not in counts[cur_rank]:
                    counts[cur_rank][cur_name] = 0
                counts[cur_rank][cur_name] += 1
        
        res_freq = dict([(i, {}) for i in range(n_ranks)])
        res_valid = dict([(i, {}) for i in range(n_ranks)])

        # collect frequency information of the names per rank
        for rank,names in counts.items():
            for name, name_counts in counts[rank].items():
                if name_counts < min_count:
                    continue
                relfreq = float(name_counts) / total_counts[rank][name]
                validfreq = float(name_counts) / n.NumTips
                res_freq[rank][name] = relfreq
                res_valid[rank][name] = validfreq

        n.ConsensusRelFreq = res_freq
        n.ValidRelFreq = res_valid
def set_ranksafe(tree, verbose=False):
    """Decorates RankSafe on tree

    RankSafe is a len(RANK_ORDER) boolean list. True means at that rank, there 
    is only a single name with > 50% relative abundance

    NOTE: tree is modified in place
    """
    if verbose:
        print "setting ranksafe nodes..."
    for node in tree.traverse(include_self=True):
        node.RankSafe = [False] * 7
        if node.istip():
            continue

        for rank, names in node.ConsensusRelFreq.items():
            # this is strict  
            if sum(map(lambda x: x >= 0.5, names.values())) == 1:
                node.RankSafe[rank] = True

def decorate_ntips(tree):
    """intelligently set the NumTips attribute on the tree"""
    n_ranks = len(RANK_ORDER)
    for node in tree.postorder(include_self=True):
        if node.istip():
            # set to True if we have consensus information
            node.NumTips = node.Consensus != ([None] * n_ranks)
        else:
            node.NumTips = reduce(add, [c.NumTips for c in node.Children])

def pick_names(tree, verbose=False):
    """Picks RankSafe names, sets RankNames on tree
    
    NOTE: tree is modified in place
    """
    if verbose:
        print "picking names..."

    for node in tree.nontips(include_self=True):
        names = []
        count = 0

        # set names at ranksafe nodes, stop if we've set a name and descendent
        # rank names are not safe.
        for rank, is_safe in enumerate(node.RankSafe):
            if is_safe:
                # place best name
                count += 1
                relfreq = node.ConsensusRelFreq[rank]
                names.append(sorted(relfreq.items(), key=itemgetter(1))[-1][0])
            else:
                # if we've had one or more useless rank, set remaining to None
                if count >= 1:
                    left = len(node.RankSafe) - len(names)
                    for i in range(left):
                        names.append(None)
                    break
                else:
                    names.append(None)
        node.RankNames = names

def fmeasure(precision, recall):
    """Returns the fmeasure (or F1-score)

    http://en.wikipedia.org/wiki/F1_score
    """
    return 2.0 * ((precision * recall) / (precision + recall))

def fpoint5measure(precision, recall):
    """Returns the f0.5measure (or F0.5-score)"""
    beta = 0.5
    betasqrd = beta**2

    tmp = (precision * recall) / ((betasqrd * precision) + recall)
    return (1 + betasqrd) * tmp

def f2measure(precision, recall):
    """Returns the f2measure (or F2-score)"""
    return (1.0+(2**2)) * ((precision * recall) / ((2**2 * precision) + recall))

def min_tips(nodes):
    """For a list of nodes, return the node with the fewest tips

    implemented naively...
    """
    scores = []
    for n in nodes:
        if n is None:
            scores.append(99999999999)
        else:
            scores.append(len(n.tips()))
    return nodes[argmin(scores)]

def name_node_score_fold(tree, score_f=fmeasure, tiebreak_f=min_tips, \
        verbose=False):
    """Compute name scores for internal nodes, pick the 'best'
    
    For this method, we traverse the tree once building up a dict of scores 
    for names and nodes, we can then pick the 'best' node out of the dict
    to avoid horrible lookups in the tree
    """
    if verbose:
        print "Starting name_node_score_fold..."
    name_node_score = dict([(i, {}) for i in range(len(RANK_ORDER))])
    n_ranks = len(RANK_ORDER)

    for node in tree.nontips(include_self=True):
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
    for rank, names in name_node_score.items():
        for name, node_scores in names.items():
            node_scores_sorted = sorted(node_scores, key=itemgetter(1))[::-1]
            nodes, scores = unzip(node_scores_sorted)
            scores = array(scores)

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
                for node,score in node_scores_sorted:
                    if node == node_to_keep:
                        continue
                    else:
                        node.RankNames[rank] = None
            else:
                for node,score in node_scores_sorted[1:]:
                    node.RankNames[rank] = None
def score_tree(tree, verbose=False):
    """Scores the tree based on RankNameScores and tip coverage

    if the node has a name in RankNames, multiply score in RankNameScores @ idx
    by the number of tips that are "valid".. Sum these scores, divide sum by 
    total number of tips covered (a tip can be counted multiple times)
    """
    total_score = 0.0
    tip_count = 0

    if verbose:
        print "Scoring tree..."

    for n in tree.nontips(include_self=True):
        for idx,name in enumerate(n.RankNames):
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

    for node in tree.nontips(include_self=True):
        node.Name = None
        node.Rank = None

        if empty_ranknames == node.RankNames:
            continue

        # take the "deepest" name
        for idx, name in enumerate(node.RankNames[::-1]):
            if name is None:
                continue
            node.Name = name
            node.Rank = n_ranks - (idx + 1) # adjust for 0-based index
            break

def make_consensus_tree(cons_split, check_for_rank=True):
    """Returns a mapping by rank for names to their parent names and counts"""
    god_node = TreeNode(Name=None)
    god_node.Rank = None
    
    base = cons_split[0]
    cur_node = god_node

    # create a base path in the tree 
    for rank, name in enumerate(base):
        new_node = TreeNode(Name=name)
        new_node.Rank = rank
        cur_node.append(new_node)
        cur_node = new_node

    # setup the initial childlookup structure so taht we don't have to
    # always iterate over .Children
    for n in god_node.traverse(include_self=True):
        if n.istip():
            n.ChildLookup = {}
            continue
        n.ChildLookup = {n.Children[0].Name:n.Children[0]}

    # for every consensus string, start at the "god" node
    for con in cons_split:
        cur_node = god_node

        # for each name, see if we've seen it, if not, add that puppy on
        for rank, name in enumerate(con):
            if name in cur_node.ChildLookup:
                cur_node = cur_node.ChildLookup[name]
            else:
                #print "adding to %s to %s" % (name, cur_node.Name)
                new_node = TreeNode(Name=name)
                new_node.Rank = rank
                new_node.ChildLookup = {}
                cur_node.append(new_node)
                cur_node.ChildLookup[name] = new_node
                cur_node = new_node

    # build an assist lookup dict
    lookup = {}
    for node in god_node.traverse():
        if node.Name is None:
            continue
        if check_for_rank and '__' in node.Name and \
                node.Name.split('__')[1] == '':
            continue
        lookup[node.Name] = node

    return god_node, lookup

def get_nearest_named_ancestor(node):
    """Returns node of nearest .Name'd ancestor

    returns None if there does not exist a named ancestor
    """
    curr = node.Parent
    while curr is not None and curr.Name is None:
        curr = curr.Parent
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
    n_ranks = len(RANK_ORDER) - 1
    for node in tree.nontips(include_self=True):
        
        if node.Name is not None:
            node.BackFillNames = [node.Name]
        else:
            node.BackFillNames = []
            continue

        if node.Parent is None:
            continue
        if node.Rank is None:
            continue
        if node.istip():
            continue
        # if node is kingdom, more on
        if node.Rank == 0:
            continue

        # find nearest ranked parent
        named_ancestor = get_nearest_named_ancestor(node)

        if named_ancestor is None:
            if verbose:
                print "Unable to find a named parent for %s" % (node.Name)
            continue

        levels = node.Rank - named_ancestor.Rank

        # this is indicative of a problem with the consensus strings
        if levels < 0:
            print node.Name, node.Rank, named_ancestor.Rank
            print '\t',node.RankSafe
            print '\t',node.RankNames
            print '\t',named_ancestor.RankSafe
            print '\t',named_ancestor.RankNames
            node.BackFillNames = []
            continue
        elif levels == 1:
            continue

        # walk consensus tree for missing names
        names = walk_consensus_tree(consensus_lookup, node.Name, levels)
        node.BackFillNames = names


def walk_consensus_tree(lookup, name, levels, reverse=True, verbose=False):
    """Walk up the consensus tree for n levels, return names
    
    if reverse is True, names are [::-1]
    """
    node = lookup[name]
    names = [name]
    curr = node.Parent

    for i in range(1,levels):
        if curr.Rank is None:
            # at root...
            break

        if curr is None:
            if verbose:
                print "Possible malformed consensus tree! See node %s" % name
            continue
        if curr.Name is None:
            if verbose:
                print "Gap in consensus tree! See node %s" % name
            names.append('%s__' % RANK_ORDER[curr.Rank])
        else:
            names.append(curr.Name)
        curr = curr.Parent

    if reverse:
        names = names[::-1]

    return names

def commonname_promotion(tree):
    """Promote names if possible from BackFillNames"""
    for node in tree.preorder(include_self=True):
        queue = [c for c in node.Children[:] if not c.istip()]
        backfill_nodes = []
        push_name = True

        # collect nearest descendants backfilled names
        while queue:
            curr = queue.pop()
            if len(curr.BackFillNames) == 0:
                queue.extend([c for c in curr.Children[:] if not c.istip()])
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

    #set_names_dict = {}
    # set the .Name attribute on the tree based on .BackFillNames
    for node in tree.preorder(include_self=True):
        if node.istip():
            continue

        if node.BackFillNames:
            unique = []
            for name in node.BackFillNames:
                unique.append(TaxaName.getTaxaName(name))
            node.BackFillNames = unique
        if len(node.BackFillNames) > 1:
            node.Name = '; '.join(node.BackFillNames)
        elif len(node.BackFillNames) == 1:
            node.Name = TaxaName.getTaxaName(node.BackFillNames[0])
        else:
            node.Name = None

class TaxaName(object):
    """Support object to provide unique taxa names"""
    _names = {}

    def __init__(self):
        raise NotImplementedError, "Object not meant for instantiation"

    @classmethod
    def getTaxaName(cls, request):
        return request
        #if '__' in request:
        #    base = request.split('__')
        #    if base[1] == '':
        #        return request
        #if request not in cls._names:
        #    cls._names[request] = 0
        #    return request
        #else:
        #    cls._names[request] += 1
        #    return '_'.join([request,str(cls._names[request])])

def make_names_unique(tree, append_suffix=True, verbose=False):
    """Appends on a unique number if multiple of the same names exist

    ordered by number of tips, ie, _1 has more tips that _2
    
    expects .BackFillNames to be set
    """
    if verbose:
        print "Assigning unique tags to duplicate names..."

    # build up a dict of the names and how many tips descend    
    name_lookup = {}
    for node in tree.nontips(include_self=True):
        if node.Name is None:
            continue
        else:
            for idx, name in enumerate(node.BackFillNames):
                if name not in name_lookup:
                    name_lookup[name] = []
                name_info = ((node.TipStop - node.TipStart),idx,node)
                name_lookup[name].append(name_info)
                #name_lookup[name].append(((node.TipStop - node.TipStart), idx, node))

    # assign unique numbers based on the number of tips that descend
    for name, scores_and_nodes in name_lookup.items():
        sorted_scores = sorted(scores_and_nodes)[::-1]
        for count, (score, idx, node) in enumerate(sorted_scores):
            # only assign a number of we have more than 1
            if count > 0:
                if node.BackFillNames[idx].split('__')[1] != '':
                    if append_suffix:
                        unique_name = '_'.join([node.BackFillNames[idx],str(count)])
                        node.BackFillNames[idx] = unique_name
                    #node.BackFillNames[idx] = '_'.join([node.BackFillNames[idx], '%d' % count])

    # should probably be refactored, but assign .Name based on .BackFillNames
    for node in tree.nontips(include_self=True):
        if len(node.BackFillNames) == 0:
            node.Name = None
        elif len(node.BackFillNames) == 1:
            node.Name = node.BackFillNames[0]
        else:
            node.Name = '; '.join(node.BackFillNames)

def pull_consensus_strings(tree, verbose=False):
    """Pulls consensus strings off of tree

    assumes .Name is set
    """
    if verbose:
        print "Pulling consensus strings..."

    constrings = []
    rank_order_rev = dict([(r,i) for i,r in enumerate(RANK_ORDER)])
    # start at the tip and travel up
    for tip in tree.tips():
        consensus_string = ['%s__' % r for r in RANK_ORDER]
        tipid = tip.Name
        n = tip.Parent

        # walk up the tree filling in the consensus string
        while n.Parent:
            if n.Name:
                if '; ' in n.Name:
                    names = n.Name.split('; ')
                    for name in names:
                        rank_idx = rank_order_rev[name[0]]
                        consensus_string[rank_idx] = name
                else:
                    rank_idx = rank_order_rev[n.Name[0]]
                    consensus_string[rank_idx] = n.Name
            n = n.Parent

        # if there is a name at the root we need to make sure we grab it
        if n.Name:
            if '; ' in n.Name:
                names = n.Name.split('; ')
                for name in names:
                    rank_idx = rank_order_rev[name[0]]
                    consensus_string[rank_idx] = name
            else:
                rank_idx = rank_order_rev[n.Name[0]]
                consensus_string[rank_idx] = n.Name

        # join strings with tip id
        constrings.append('\t'.join([tipid, '; '.join(consensus_string)]))
    return constrings

def save_bootstraps(tree, verbose=False):
    """Retains .Bootstrap if set in .Name"""
    if verbose:
        print "Attempting to retain bootstrap values"
    for n in tree.nontips(include_self=True):
        if n.Bootstrap is not None:
            if n.Name is None:
                n.Name = str(n.Bootstrap)
            else:
                if ';' in n.Name:
                    n.Name = ':'.join([str(n.Bootstrap), "%s" % n.Name])
                else:
                    n.Name = ':'.join([str(n.Bootstrap), n.Name])

def getpath(foo):
    while foo.Parent:
        print foo.Name
        if hasattr(foo, 'RankNames'):
            print foo.RankNames
        foo = foo.Parent
    print foo.Name
    if hasattr(foo, 'RankNames'):
        print foo.RankNames

def is_float(s):
    """Returns True if the value in string s is a float"""
    if s is None:
        return False

    if '.' not in s:
        return False

    try:
        test = float(s)
        return True
    except:
        return False


def validate_all_paths(tree):
    """Walk each path in the tree and make sure there aren't any conflicts"""
    # helper getpath method
    def getpath_f(n):
        path = []
        n = n.Parent
        while n.Parent:
            if n.Name is not None:
                if ':' in n.Name:
                    names = n.Name.split(':',1)
                    if not is_float(names[1]):
                        path.append(names[1])
                else:
                    if not is_float(n.Name):
                        path.append(n.Name)
            n = n.Parent
        if n.Name is not None:
            path.append(n.Name)

        clean_path = []
        for p in path:
            if '; ' in p:
                clean_path.extend(p.split('; ')[::-1])
            else:
                clean_path.append(p)
        return clean_path

    rank_order_rev = dict([(r,i) for i,r in enumerate(RANK_ORDER)])

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

def nlevel_workflow(tree, contree_lookup, verbose=False):
    counts = collect_names_at_ranks_counts(tree,verbose)
    decorate_ntips(tree)
    print "DECORATING WITH NTIPS, SHOULD USE RANGENODE METHODS"
    min_count = 2
    decorate_name_relative_freqs(tree, counts, min_count, verbose)
    set_ranksafe(tree, verbose)
    pick_names(tree, verbose)
    name_node_score_fold(tree, verbose)
    set_preliminary_name_and_rank(tree)
    backfill_names_gap(tree, contree_lookup, verbose)
    commonname_promotion(tree)
    make_names_unique(tree, verbose)
    return tree

from sys import argv
def main(args=argv):
    parser = OptionParser(option_list=options)
    opts, args = parser.parse_args(args=args)
    tipname_map = load_consensus_map(open(opts.consensus_map), 
                                     opts.append_rank,
                                     opts.verbose) 
    tree = load_tree(open(opts.tree), tipname_map, opts.verbose) 
    counts = collect_names_at_ranks_counts(tree, opts.verbose)
    decorate_ntips(tree)
    min_count = 2
    decorate_name_relative_freqs(tree, counts, min_count, opts.verbose)
    set_ranksafe(tree, opts.verbose)
    pick_names(tree, opts.verbose)
    name_node_score_fold(tree, verbose=opts.verbose)
    
    if opts.verbose:
        print "SCORE: ", score_tree(tree, opts.verbose)
    set_preliminary_name_and_rank(tree)

    contree, contree_lookup = make_consensus_tree(tipname_map.values())
    backfill_names_gap(tree, contree_lookup, opts.verbose)
    commonname_promotion(tree)
    make_names_unique(tree, append_suffix=False, verbose=opts.verbose)

    
    constrings = pull_consensus_strings(tree, opts.verbose)

    f = open(opts.output + '-consensus-strings', 'w')
    f.write('\n'.join(constrings))
    f.close()

    save_bootstraps(tree, opts.verbose)
    f = open(opts.output, 'w')
    f.write(tree.getNewick(with_distances=True))
    f.close()

                
if __name__ == '__main__':
    main()


