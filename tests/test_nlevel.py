#!/usr/bin/env python

"""nlevel tests

While many of the tests utilize similar trees and input data, the overlap
is not necessarily 100%. Many of these inputs are written with specific tests
in mind.
"""

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

from unittest import TestCase, main
from t2t.nlevel import (load_consensus_map, collect_names_at_ranks_counts,
                        load_tree, decorate_name_relative_freqs, decorate_name_counts,
                        set_ranksafe,
                        pick_names, has_badname, get_nearest_named_ancestor,
                        walk_consensus_tree, make_consensus_tree,
                        backfill_names_gap, commonname_promotion,
                        decorate_ntips, decorate_ntips_rank,
                        name_node_score_fold, backfill_from_secondary,
                        validate_all_paths, score_tree,
                        promote_to_multifurcation,
                        recover_from_polyphyletic_sibling,
                        _named_siblings, polyphyletic_unique,
                        normalize_species_binomial,
                        correct_species_binomial)

from skbio import TreeNode
import sys

if sys.version_info[0] == 2:
    from StringIO import StringIO
else:
    from io import StringIO

class NLevelTests(TestCase):

    def setUp(self):
        pass

    def test_score_tree(self):
        """Determine's the tree's fmeasure score"""
        # set RankNames and RankNameScores
        # if name in RankNames, check score, look at tips, etc
        t_str = StringIO(u"(((a,b),(c,d))e,(f,g),h)i;")
        t = TreeNode.read(t_str)
        t.RankNames = ['i', None, None, None]  # 1.0 * 6
        t.RankNameScores = [1.0, None, None, None]
        t.children[0].RankNames = [None, 'e', 'foo', None]  # 0.5 * 3, 0.6 * 3
        t.children[0].RankNameScores = [None, 0.5, 0.6, None]
        t.children[0].children[0].RankNames = [None] * 7
        t.children[0].children[1].RankNames = [None] * 7
        t.children[1].RankNames = [None] * 7
        t.children[1].RankNameScores = [None] * 7
        tips = list(t.tips())
        tips[0].Consensus = [None] * 7
        tips[1].Consensus = [1, 3, None, None]
        tips[2].Consensus = [2, 4, 5, None]
        tips[3].Consensus = [None, 1, None, None]
        tips[4].Consensus = [None, 1, None, None]
        tips[5].Consensus = [2, None, 3, None]
        tips[6].Consensus = [None, 4, None, None]
        decorate_ntips(t)
        exp = ((1.0 * 6) + (0.5 * 3) + (0.6 * 3)) / (6 + 3 + 3)
        obs = score_tree(t)
        self.assertEqual(obs, exp)

    def test_has_badname(self):
        """correctly determines if a string is bad"""
        data = "assadasd environmental sample asdasdd"
        self.assertTrue(has_badname(data))
        data = "asdasdsad dsasda dasd as"
        self.assertFalse(has_badname(data))

    def test_load_consensus_map(self):
        """correctly returns a consensus map"""
        data = ["foo\ta; b; c; d; e; f; g",
                 "bar\th; i; j; k; l; m; n",
                 "foobar\th; i; j; None; l; ; foo uncultured bar"]
        exp_noappend = {'foo': ['a', 'b', 'c', 'd', 'e', 'f', 'g'],
                        'bar': ['h', 'i', 'j', 'k', 'l', 'm', 'n'],
                        'foobar': ['h', 'i', 'j', None, 'l', None, None]}
        exp_append = {
            'foo': ['d__a', 'p__b', 'c__c', 'o__d', 'f__e', 'g__f', 's__g'],
            'bar': ['d__h', 'p__i', 'c__j', 'o__k', 'f__l', 'g__m', 's__n'],
            'foobar': ['d__h', 'p__i', 'c__j', 'o__', 'f__l', 'g__', 's__']}
        obs_noappend = load_consensus_map(data, False)
        obs_append = load_consensus_map(data, True)
        self.assertEqual(obs_noappend, exp_noappend)
        self.assertEqual(obs_append, exp_append)

    def test_collect_names_at_ranks_counts(self):
        """correctly returns total counts for names at ranks"""
        data = StringIO(u"((a,b)c,(d,(e,f)g)h,(i,j)k)l;")
        tipname_map = {'a': ['1', '2', '3', '4', '5', '6', '7'],
                       'b': ['1', '2', '3', None, '5', '6', '8'],
                       'd': ['1', '2', '3', 'a', '5', '6', '9'],
                       'e': ['1', '2', '3', None, '5', '6', '9'],
                       'i': ['1', '2', '3', 'a', '5', '6', '9'],
                       'j': ['1', '2', '3', '4', '5', '6', '9']}
        tree = load_tree(data, tipname_map)

        exp = {0: {'1': 6},
               1: {'2': 6},
               2: {'3': 6},
               3: {'4': 2, 'a': 2},
               4: {'5': 6},
               5: {'6': 6},
               6: {'7': 1, '8': 1, '9': 4}}

        obs = collect_names_at_ranks_counts(tree)
        self.assertEqual(obs, exp)

    def test_load_tree(self):
        """correctly loads and decorates tiplook info on a tree"""
        data = StringIO(u"((a,b)c,(d,(e,f)g)h,(i,j)k)l;")
        tipname_map = {'a': ['1', '2', '3', '4', '5', '6', '7'],
                       'b': ['1', '2', '3', '4', '5', '6', '8'],
                       'd': ['1', '2', '3', '4', '5', '6', '9'],
                       'e': ['1', '2', '3', '4', '5', '6', '10'],
                       'i': ['1', '2', '3', '4', '5', '6', '12'],
                       'j': ['1', '2', '3', '4', '5', '6', '13']}

        # exp in Name: (tipstart, tipstop, consensus)
        exp = {'a': (0, 0, ['1', '2', '3', '4', '5', '6', '7']),
               'b': (1, 1, ['1', '2', '3', '4', '5', '6', '8']),
               'c': (0, 1, [None] * 7),
               'd': (2, 2, ['1', '2', '3', '4', '5', '6', '9']),
               'e': (3, 3, ['1', '2', '3', '4', '5', '6', '10']),
               'f': (4, 4, [None] * 7),
               'g': (3, 4, [None] * 7),
               'h': (2, 4, [None] * 7),
               'i': (5, 5, ['1', '2', '3', '4', '5', '6', '12']),
               'j': (6, 6, ['1', '2', '3', '4', '5', '6', '13']),
               'k': (5, 6, [None] * 7),
               'l': (0, 6, [None] * 7)}

        obstree = load_tree(data, tipname_map)
        obs = {}
        for node in obstree.traverse(include_self=True):
            obs[node.name] = (node.TipStart, node.TipStop, node.Consensus)

        self.assertEqual(obs, exp)

    def test_decorate_name_relative_freqs(self):
        """correctly decorate relative frequency information on a tree"""
        data = StringIO(u"((a,b)c,(d,(e,f)g)h,(i,j)k)l;")
        tipname_map = {'a': ['1', '2', '3', '4', '5', '6', '7'],
                       'b': ['1', '2', '3', '4', '5', '6', '8'],
                       'd': ['1', '2', '3', '4', '5', '6', '8'],
                       'e': ['1', '2', '3', '4', 'a', '6', '7'],
                       'i': ['1', '2', '3', '4', 'a', None, '7'],
                       'j': ['1', '2', '3', '4', 'a', None, '8']}

        total_counts = {0: {'1': 10, 'foo': 5},
                        1: {'2': 6},
                        2: {'3': 12},
                        3: {'4': 6, 'bar': 5},
                        4: {'5': 6, 'a': 3},
                        5: {'6': 6},
                        6: {'7': 3, '8': 3}}

        tree = load_tree(data, tipname_map)
        decorate_ntips(tree)
        decorate_name_relative_freqs(tree, total_counts, 1)

        exp_root = {0: {'1': .6},
                    1: {'2': 1.0},
                    2: {'3': .5},
                    3: {'4': 1.0},
                    4: {'5': .5, 'a': 1.0},
                    5: {'6': 4.0 / 6},
                    6: {'7': 1.0, '8': 1.0}}

        self.assertEqual(tree.ConsensusRelFreq, exp_root)

    def test_decorate_name_counts(self):
        """correctly decorate relative frequency information on a tree"""
        data = StringIO(u"((a,b)c,(d,(e,f)g)h,(i,j)k)l;")
        tipname_map = {'a': ['1', '2', '3', '4', '5', '6', '7'],
                       'b': ['1', '2', '3', '4', '5', '6', '8'],
                       'd': ['1', '2', '3', '4', '5', '6', '8'],
                       'e': ['1', '2', '3', '4', 'a', '6', '7'],
                       'i': ['1', '2', '3', '4', 'a', None, '7'],
                       'j': ['1', '2', '3', '4', 'a', None, '8']}

        tree = load_tree(data, tipname_map)
        decorate_ntips(tree)
        decorate_name_counts(tree)

        exp_root = {0: {'1': 6},
                    1: {'2': 6},
                    2: {'3': 6},
                    3: {'4': 6},
                    4: {'5': 3, 'a': 3},
                    5: {'6': 4},
                    6: {'7': 3, '8': 3}}

        self.assertEqual(tree.TaxaCount, exp_root)

    def test_set_ranksafe(self):
        """correctly set ranksafe on tree"""
        data = StringIO(u"((a,b)c,(d,(e,f)g)h,(i,j)k)l;")
        tipname_map = {'a': ['1', '2', '3', '4', '5', '6', '7'],
                       'b': ['1', '2', '3', 'b', '5', '6', '8'],
                       'd': ['1', '2', '3', '4', '5', '6', '8'],
                       'e': ['1', '2', '3', 'b', 'a', '6', '7'],
                       'i': ['1', '2', '3', '4', 'a', '6', '7'],
                       'j': ['1', '2', '3', 'b', 'a', '6', '8']}

        total_counts = {0: {'1': 10, 'foo': 5},
                        1: {'2': 20},
                        2: {'3': 10},
                        3: {'4': 4, 'b': 5},
                        4: {'5': 7, 'a': 5},
                        5: {'6': 6},
                        6: {'7': 3, '8': 3}}

        tree = load_tree(data, tipname_map)
        decorate_ntips(tree)
        decorate_name_relative_freqs(tree, total_counts, 1)
        set_ranksafe(tree)

        exp_root = [True, False, True, False, True, True, False]
        self.assertEqual(tree.RankSafe, exp_root)

    def test_name_node_score_fold(self):
        """hate taxonomy"""
        data = StringIO(u"((a,b)c,(d,(e,f)g)h,(i,j)k)l;")
        tipname_map = {'a': ['1', '2', '3', '4', '5', '6', '8'],
                       'b': ['1', '2', '3', '4', '5', '6', '8'],
                       'd': ['1', '2', '3', 'f', 'e', 'c', '9'],
                       'e': ['1', '2', '3', 'f', 'e', 'c', '9'],
                       'i': ['1', '2', '3', 'g', 'a', 'h', '11'],
                       'j': ['1', '2', '3', 'g', 'a', 'h', '12']}

        total_counts = {0: {'1': 10, 'foo': 5},
                        1: {'2': 10},
                        2: {'3': 10},
                        3: {'4': 3, 'f': 5, 'g': 5},
                        4: {'5': 7, 'a': 5, 'e': 4},
                        5: {'6': 3, 'c': 3, 'd': 2, 'h': 3},
                        6: {'8': 3, '9': 2, '10': 2, '11': 2, '12': 2}}

        tree = load_tree(data, tipname_map)
        decorate_ntips(tree)
        decorate_name_relative_freqs(tree, total_counts, 1)
        set_ranksafe(tree)
        pick_names(tree)
        name_node_score_fold(tree)
        exp_root = ['1', '2', '3', None, None, None, None]
        expc0 = [None, None, None, '4', None, None, None]
        expc1 = [None, None, None, None, 'e', 'c', '9']
        expc2 = [None, None, None, None, None, 'h', None]
        expc1c1 = [None] * 7

        self.assertEqual(tree.RankNames, exp_root)
        self.assertEqual(tree.children[0].RankNames, expc0)
        self.assertEqual(tree.children[1].RankNames, expc1)
        self.assertEqual(tree.children[2].RankNames, expc2)
        self.assertEqual(tree.children[1].children[1].RankNames, expc1c1)

    def test_validate_all_paths(self):
        """complains correctly about badpaths"""
        data = StringIO(u"(((((1,2)s__,(3,4)s__)g__)p__),((5,6)f__)f__,((7,8)c__)o__);")
        t = load_tree(data, {})
        exp = [n for n in t.tips() if n.name in ['5', '6', '7', '8']]
        obs = validate_all_paths(t)
        self.assertEqual(obs, exp)

    def test_best_name_freqs_for_nodes(self):
        """correctly gets the frequencies per name per node"""
        data = StringIO(u"((a,b)c,(d,(e,f)g)h,(i,j)k)l;")
        tipname_map = {'a': ['1', '2', '3', '4', '5', '6', '7'],
                       'b': ['1', '2', '3', 'b', '5', '6', '8'],
                       'd': ['1', '2', '3', '4', '5', '6', '7'],
                       'e': ['1', '2', '3', 'b', 'a', 'foo', '7'],
                       'i': ['1', '2', '3', '4', 'a', 'foo', '8'],
                       'j': ['1', '2', '3', 'b', 'a', 'foo', '8']}

        total_counts = {0: {'1': 10, 'foo': 5},
                        1: {'2': 10},
                        2: {'3': 10},
                        3: {'4': 4, 'b': 5},
                        4: {'5': 7, 'a': 5},
                        5: {'6': 3, 'foo': 3},
                        6: {'7': 3, '8': 4}}

        tree = load_tree(data, tipname_map)
        decorate_ntips(tree)
        decorate_name_relative_freqs(tree, total_counts, 1)
        set_ranksafe(tree)
        pick_names(tree)

        # result = best_name_freqs_for_nodes(tree)
        cnode = tree.children[0]
        hnode = tree.children[1]
        knode = tree.children[2]
        # exp = {0:{'1':[(0.6,tree)]}
        #        1:{'2':[(0.6,tree)]},
        #        2:{'3':[(0.6,tree)]},
        #        3:{'4':[(0.75,tree)],'b':[(0.6,tree)]},
        #        4:{'a':[(0.6, tree),[(0.5,knode)]

    def test_pick_names(self):
        """correctly pick names to retain on a tree"""
        data = StringIO(u"((a,b)c,(d,(e,f)g)h,(i,j)k)l;")
        tipname_map = {'a': ['1', '2', '3', '4', '5', '6', '7'],
                       'b': ['1', '2', '3', 'b', '5', '6', '8'],
                       'd': ['1', '2', '3', '4', '5', '6', '7'],
                       'e': ['1', '2', '3', 'b', 'a', 'foo', '7'],
                       'i': ['1', '2', '3', '4', 'a', 'foo', '8'],
                       'j': ['1', '2', '3', 'b', 'a', 'foo', '8']}

        total_counts = {0: {'1': 10, 'foo': 5},
                        1: {'2': 10},
                        2: {'3': 10},
                        3: {'4': 4, 'b': 5},
                        4: {'5': 7, 'a': 5},
                        5: {'6': 3, 'foo': 3},
                        6: {'7': 3, '8': 4}}

        tree = load_tree(data, tipname_map)
        decorate_ntips(tree)
        decorate_name_relative_freqs(tree, total_counts, 1)
        set_ranksafe(tree)
        pick_names(tree)
        exp_root = ['1', '2', '3', None, None, None, None]
        self.assertEqual(tree.RankNames, exp_root)

        expc0 = [None, None, None, None, None, '6', None]
        expc1 = [None, None, None, None, None, None, '7']
        expc2 = [None, None, None, None, None, 'foo', '8']
        expc1c1 = [None] * 7

        self.assertEqual(tree.children[0].RankNames, expc0)
        self.assertEqual(tree.children[1].RankNames, expc1)
        self.assertEqual(tree.children[2].RankNames, expc2)
        self.assertEqual(tree.children[1].children[1].RankNames, expc1c1)

    def test_walk_consensus_tree(self):
        """correctly walk consensus tree"""
        data = [['a', 'b', 'c', 'd', 'e', 'f', 'g'],
                 ['a', 'b', 'c', None, None, 'x', 'y'],
                 ['h', 'i', 'j', 'k', 'l', 'm', 'n'],
                 ['h', 'i', 'j', 'k', 'l', 'm', 'q'],
                 ['h', 'i', 'j', 'k', 'l', 'm', 'n']]
        _, lookup = make_consensus_tree(data, check_for_rank=False)
        exp1 = ['n', 'm', 'l']
        exp2 = ['a']
        exp3 = ['x', 'f__', 'o__', 'c']
        obs1 = walk_consensus_tree(lookup, 'n', 3, reverse=False)
        obs2 = walk_consensus_tree(lookup, 'a', 3, reverse=False)
        obs3 = walk_consensus_tree(lookup, 'x', 4, reverse=False)
        self.assertEqual(obs1, exp1)
        self.assertEqual(obs2, exp2)
        self.assertEqual(obs3, exp3)

    def test_make_consensus_tree(self):
        """correctly gets parent consensus counts"""
        data = [['a', 'b', 'c', 'd', 'e', 'f', 'g'],
                 ['a', 'b', 'c', None, None, 'x', 'y'],
                 ['h', 'i', 'j', 'k', 'l', 'm', 'n'],
                 ['h', 'i', 'j', 'k', 'l', 'm', 'q'],
                 ['h', 'i', 'j', 'k', 'l', 'm', 'n']]
        exp_str = "(((((((g)f)e)d,(((y)x)))c)b)a,((((((n,q)m)l)k)j)i)h);"

        obs_root, lookup = make_consensus_tree(data, check_for_rank=False)

        self.assertEqual(str(obs_root).strip(), exp_str)
        self.assertNotIn(None, lookup)

    def test_make_consensus_tree_withtips(self):
        """correctly constructs the taxonomy tree with tip info"""
        data = [['a', 'b', 'c', 'd', 'e', 'f', 'g'],
                 ['a', 'b', 'c', None, None, 'x', 'y'],
                 ['h', 'i', 'j', 'k', 'l', 'm', 'n'],
                 ['h', 'i', 'j', 'k', 'l', 'm', 'q'],
                 ['h', 'i', 'j', 'k', 'l', 'm', 'n']]
        input_ids = ['1', '2', '3', '4', '5']
        exp_str = ("((((((((1)g)f)e)d,((((2)y)x)))c)b)a,(((((((3,5)n,"
                   "(4)q)m)l)k)j)i)h);")

        obs_root, lookup = make_consensus_tree(
            data, check_for_rank=False, tips=input_ids)

        self.assertEqual(str(obs_root).strip(), exp_str)
        self.assertNotIn(None, lookup)

    def test_decorate_ntips(self):
        """correctly decorate the tree with the NumTips param"""
        data = StringIO(u"(((a,b)c,(d,e,f)g)h,(i,j)k)l;")
        tree = TreeNode.read(data)
        tips = dict([(tip.name, tip) for tip in tree.tips()])
        tips['a'].Consensus = [1, 2, 3, 4, 5, 6, 7]
        tips['b'].Consensus = [None, None, None, 5, None, None, None]
        tips['d'].Consensus = [1, 2, 3, 4, 5, 6, 8]
        tips['e'].Consensus = [None, None, None, None, None, None, None]
        tips['f'].Consensus = [1, 2, 3, 4, 5, 6, 8]
        tips['i'].Consensus = [1, 2, 3, 4, 5, 6, 8]
        tips['j'].Consensus = [1, 2, 3, 4, 5, 6, 8]
        decorate_ntips(tree)
        self.assertEqual(tree.NumTips, 6)
        self.assertEqual(tree.children[0].NumTips, 4)
        self.assertEqual(tree.children[1].NumTips, 2)
        self.assertEqual(tree.children[0].children[0].NumTips, 2)
        self.assertEqual(tree.children[0].children[1].NumTips, 2)

    def test_decorate_ntips_rank(self):
        """correctly decorate the tree with the NumTipsRank param"""
        data = StringIO(u"(((a,b)c,(d,e,f)g)h,(i,j)k)l;")
        tree = TreeNode.read(data)
        tips = dict([(tip.name, tip) for tip in tree.tips()])
        tips['a'].Consensus = [1, 2, 3, 4, 5, 6, 7]
        tips['b'].Consensus = [None, None, None, 5, None, None, None]
        tips['d'].Consensus = [1, 2, 3, 4, 5, 6, 8]
        tips['e'].Consensus = [None, None, None, None, None, None, None]
        tips['f'].Consensus = [1, 2, 3, 4, 5, 6, 8]
        tips['i'].Consensus = [1, 2, 3, 4, 5, 6, 8]
        tips['j'].Consensus = [1, 2, 3, 4, 5, 6, 8]
        decorate_ntips_rank(tree)

        self.assertEqual(tree.NumTipsRank[0], 5)
        self.assertEqual(tree.NumTipsRank[1], 5)
        self.assertEqual(tree.NumTipsRank[2], 5)
        self.assertEqual(tree.NumTipsRank[3], 6)
        self.assertEqual(tree.NumTipsRank[4], 5)
        self.assertEqual(tree.NumTipsRank[5], 5)
        self.assertEqual(tree.NumTipsRank[6], 5)

    def test_get_nearest_named_ancestor(self):
        """correctly get the nearest named ancestor"""
        t = TreeNode.read(StringIO(u"(((s1,s2)g1,s3))root;"))
        t2 = TreeNode.read(StringIO(u"(((s1,s2)g1,s3));"))
        exp_t = t
        exp_t2 = None
        obs_t = get_nearest_named_ancestor(t.find('s3'))
        obs_t2 = get_nearest_named_ancestor(t2.find('s3'))
        self.assertEqual(obs_t, exp_t)
        self.assertEqual(obs_t2, exp_t2)

    def test_backfill_names_gap(self):
        """correctly backfill names"""
        consensus_tree = TreeNode.read(StringIO(u"(((s1,s2)g1,(s3,s4)g2,(s5,s6)g3)f1)o1;"))
        rank_lookup = {'s': 6, 'g': 5, 'f': 4, 'o': 3, 'c': 2, 'p': 1, 'k': 0}
        for n in consensus_tree.traverse(include_self=True):
            n.Rank = rank_lookup[n.name[0]]
        data = StringIO(u"((((1)s1,(2)s2),((3)s3,(4)s5)))o1;")
        lookup = dict([(n.name, n)
                      for n in consensus_tree.traverse(include_self=True)])
        # exp = "((((1)s1,(2)s2)g1,((3)'g2; s3',(4)'g3; s5')))'o1; f1'"
        t = TreeNode.read(data)
        t.Rank = 3
        t.children[0].Rank = None
        t.children[0].children[0].Rank = None
        t.children[0].children[1].Rank = None
        t.children[0].children[0].children[0].Rank = 6
        t.children[0].children[0].children[1].Rank = 6
        t.children[0].children[1].children[0].Rank = 6
        t.children[0].children[1].children[1].Rank = 6

        backfill_names_gap(t, lookup)

        self.assertEqual(t.BackFillNames, ['o1'])
        self.assertEqual(t.children[0].BackFillNames, [])
        self.assertEqual(t.children[0].children[0].BackFillNames, [])
        self.assertEqual(t.children[0].children[1].BackFillNames, [])
        self.assertEqual(t.children[0].children[0]
                         .children[0].BackFillNames, ['f1', 'g1', 's1'])
        self.assertEqual(t.children[0].children[0]
                         .children[1].BackFillNames, ['f1', 'g1', 's2'])
        self.assertEqual(t.children[0].children[1]
                         .children[0].BackFillNames, ['f1', 'g2', 's3'])
        self.assertEqual(t.children[0].children[1]
                         .children[1].BackFillNames, ['f1', 'g3', 's5'])

    def test_backfill_names_dangling(self):
        """correctly fill in dangling missing ranks"""
        consensus_tree = TreeNode.read(StringIO(u"(((s1,s2)g1,(s3,s4)g2,(s5,s6)g3)f1)o1;"))
        data = "((((1),(2)),((3),(4))))'o1; f1';"
        lookup = dict([(n.name, n)
                      for n in consensus_tree.traverse(include_self=True)])
        # exp = "((((1),(2)),((3),(4))))'o1; f1';"

    def test_polyphyletic_unique(self):
        tests = [(['s__a', 's__b'], {'s__a', 's__b'}),
                 (['s__a', 's__b', 's__c'], {'s__a', 's__b', 's__c'}),
                 (['s__a'], {'s__a', }),
                 (['s__a', 's__a_AA'], {'s__a_AA', }),
                 (['s__a', 's__a_AA', 's__a_AB'], {'s__a', 's__a_AA', 's__a_AB'})]
        for in_, exp in tests:
            obs = polyphyletic_unique(in_)
            self.assertEqual(obs, exp)

    def test_correct_species_binomial(self):
        t = TreeNode.read(["(((1,2)'s__foo bar',(3,4)'s__foo baz')g__foo_A,(5,6)g__baz,(((7)'s__baz thing',8)g__baz))root;"],  # noqa
                          convert_underscores=False)
        exp = TreeNode.read(["(((1,2)'s__foo_A bar',(3,4)'s__foo_A baz')g__foo_A,(5,6)g__baz,(((7)'s__baz thing',8)g__baz))root;"],  # noqa
                            convert_underscores=False)
        obs = correct_species_binomial(t)
        self.assertEqual(str(obs), str(exp))

    def test_normalize_species_binomial(self):
        tests = [(('g__foo', 's__foo bar'), 's__foo bar'),
                 (('g__foo_A', 's__foo bar'), 's__foo_A bar'),
                 (('g__foo_A', 's__baz bar'), 's__baz bar'),  # do not create new species names!  # noqa
                 (('g__foo_A', 's__foo_A bar'), 's__foo_A bar'),
                 (('g__foo_A', 's__foo_A bar_A'), 's__foo_A bar_A'),
                 (('g__foo', 's__foo bar_A'), 's__foo bar_A')]
        for in_, exp in tests:
            obs = normalize_species_binomial(*in_)
            self.assertEqual(obs, exp)

        tests = [('g__foo_A', 's__foo_B bar'), ]
        for in_ in tests:
            with self.assertRaises(ValueError):
                normalize_species_binomial(*in_)

    def test_named_siblings(self):
        t = TreeNode.read(['(((1,2)s__A,(3,4)s__B)g__A,(5,6)g__B,((7,8)g__C))root;'],  # noqa
                          convert_underscores=False)
        t.assign_ids()
        exp = [t.find('g__A'), t.find('g__C')]
        obs = _named_siblings(t.find('g__B'))
        self.assertEqual(obs, exp)

    def test_recover_from_polyphyletic_sibling(self):
        # map species and genus including on nested name
        # do not map a substring (g__baz -> g__bazx)
        t = TreeNode.read(["((((1,2)s__foo,(3,4)s__foo_A)'g__bar; f__baz',(5)'s__x; g__bar_A'),(6,7)g__baz,(8,9)g__bazx)root;"],  # noqa
                          convert_underscores=False)
        exp = TreeNode.read(["((((1,2)s__foo_A,(3,4)s__foo_A)'g__bar_A; f__baz',(5)'s__x; g__bar_A'),(6,7)g__baz,(8,9)g__bazx)root;"],  # noqa
                            convert_underscores=False)
        for n in t.traverse(include_self=False):
            n.length = 0
        for n in exp.traverse(include_self=False):
            n.length = 0
        obs = recover_from_polyphyletic_sibling(t)
        self.assertEqual(str(obs), str(exp))

        # two options so do not map
        t = TreeNode.read(['((1,2)s__foo:0.1,(3,4)s__foo_A:0.2,(5,6)s__foo_B:0.3)root;'],
                          convert_underscores=False)
        exp = TreeNode.read(['((1,2)s__foo_A:0.1,(3,4)s__foo_A:0.2,(5,6)s__foo_B:0.3)root;'],
                            convert_underscores=False)
        obs = recover_from_polyphyletic_sibling(t)
        self.assertEqual(str(obs), str(exp))

    def test_backfill_from_secondary(self):
        secondary_taxonomy = TreeNode.read(["((((l1,l2)s1,(l3)s2)g1,((l4)s3,(l5)s4)g2,((l6)s5,(l7)s6)g3)f1)o1;"])  # noqa
        rank_lookup = {'s': 6, 'g': 5, 'f': 4, 'o': 3, 'c': 2, 'p': 1, 'd': 0}
        for n in secondary_taxonomy.non_tips():
            n.Rank = rank_lookup[n.name[0]]

        # setup a bunch of state on the tree
        tree = TreeNode.read(["(((X1,X2)s1,(l1,l2))f1,(((l3)),(l4,X3))f2)c1;"])
        for n in tree.non_tips(include_self=True):
            n.BackFillNames = []
        tree.Rank = 2
        tree.children[0].Rank = 4  # ((X1,X2),(l1,l2))
        tree.children[0].BackFillNames = ['oX', 'f1']
        tree.children[0].children[0].Rank = 6  # (X1,X2)
        tree.children[0].children[0].BackFillNames = ['gX', 's1']
        tree.children[0].children[1].Rank = None  # (l1,l2)
        tree.children[1].Rank = 4  # (((l3)),(l4,X3))
        tree.children[1].BackFillNames = ['f2']
        tree.children[1].children[0].Rank = None  # ((l3))
        tree.children[1].children[0].children[0].Rank = None  # (l3)
        tree.children[1].children[1].Rank = None  # (l4,X3)
        backfill_from_secondary(tree, secondary_taxonomy)

        n = tree.find('l1').parent
        self.assertEqual(n.BackFillNames, ['g1', 's1'])
        n = tree.find('l3').parent
        self.assertEqual(n.BackFillNames, ['g1', 's2'])
        n = tree.find('l4').parent
        self.assertEqual(n.BackFillNames, ['g2', 's3'])

    def test_commonname_promotion(self):
        """correctly promote names if possible"""
        consensus_tree = TreeNode.read(StringIO(u"(((s1,s2)g1,(s3,s4)g2,(s5,s6)g3)f1)o1;"))
        rank_lookup = {'s': 6, 'g': 5, 'f': 4, 'o': 3, 'c': 2, 'p': 1, 'k': 0}
        for n in consensus_tree.traverse(include_self=True):
            n.Rank = rank_lookup[n.name[0]]
        data = StringIO(u"((((1)s1,(2)s2),((3)s3,(4)s5)))o1;")
        lookup = dict([(n.name, n)
                      for n in consensus_tree.traverse(include_self=True)])
        exp = "((((1)s1,(2)s2)g1,((3)'g2; s3',(4)'g3; s5')))'o1; f1';"
        t = TreeNode.read(data)
        t.Rank = 3
        t.children[0].Rank = None
        t.children[0].children[0].Rank = None
        t.children[0].children[1].Rank = None
        t.children[0].children[0].children[0].Rank = 6
        t.children[0].children[0].children[1].Rank = 6
        t.children[0].children[1].children[0].Rank = 6
        t.children[0].children[1].children[1].Rank = 6
        backfill_names_gap(t, lookup)
        commonname_promotion(t)

        self.assertEqual(str(t).rstrip(), exp)

    def test_promote_to_multifurcation(self):
        tree = TreeNode.read(["""(((a,b)s__1,(c,d)s__2)g__1,
                                  ((e,f)s__3,(g,h)));"""],
                             convert_underscores=False)
        exp = TreeNode.read(["""(((a,b)s__1,(c,d)s__2)g__1,
                                 ((e,f),(g,h))s__3);"""],
                             convert_underscores=False)
        fragset = {'g', 'h'}
        obs = promote_to_multifurcation(tree, fragset)
        self.assertEqual(str(obs), str(exp))


if __name__ == '__main__':
    main()
