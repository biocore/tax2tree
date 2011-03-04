#!/usr/bin/env python

"""nlevel tests

While many of the tests utilize similar trees and input data, the overlap
is not necessarily 100%. Many of these inputs are written with specific tests
in mind.
"""

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from t2t.nlevel import load_consensus_map, collect_names_at_ranks_counts, \
        load_tree, decorate_name_relative_freqs, set_ranksafe, pick_names, \
        has_badname, get_nearest_named_ancestor, walk_consensus_tree,\
        make_consensus_tree, backfill_names_gap, commonname_promotion, \
        decorate_ntips, name_node_score_fold, validate_all_paths, \
        score_tree
from cogent.parse.tree import DndParser
from cogent.core.tree import TreeNode

class NLevelTests(TestCase):
    def setUp(self):
        pass

    def test_score_tree(self):
        """Determine's the tree's fmeasure score"""
        # set RankNames and RankNameScores
        # if name in RankNames, check score, look at tips, etc
        t_str = "(((a,b),(c,d))e,(f,g),h)i;"
        t = DndParser(t_str)
        t.RankNames = ['i',None,None,None] # 1.0 * 6
        t.RankNameScores = [1.0,None,None,None]
        t.Children[0].RankNames = [None,'e','foo',None] # 0.5 * 3, 0.6 * 3
        t.Children[0].RankNameScores = [None, 0.5, 0.6, None]
        t.Children[0].Children[0].RankNames = [None] * 7
        t.Children[0].Children[1].RankNames = [None] * 7
        t.Children[1].RankNames = [None] * 7
        t.Children[1].RankNameScores = [None] * 7
        tips = t.tips()
        tips[0].Consensus = [None] * 7
        tips[1].Consensus = [1,3,None,None]
        tips[2].Consensus = [2,4,5,None]
        tips[3].Consensus = [None,1,None,None]
        tips[4].Consensus = [None,1,None,None]
        tips[5].Consensus = [2,None,3,None]
        tips[6].Consensus = [None,4,None,None]
        decorate_ntips(t)
        exp = ((1.0 * 6) + (0.5 * 3) + (0.6 * 3)) / (6 + 3 + 3)
        obs = score_tree(t)
        self.assertEqual(obs, exp)

    def test_has_badname(self):
        """correctly determines if a string is bad"""
        input = "assadasd environmental sample asdasdd"
        self.assertTrue(has_badname(input))
        input = "asdasdsad dsasda dasd as"
        self.assertFalse(has_badname(input))

    def test_load_consensus_map(self):
        """correctly returns a consensus map"""
        input = ["foo\ta; b; c; d; e; f; g",
                 "bar\th; i; j; k; l; m; n",
                 "foobar\th; i; j; None; l; ; foo uncultured bar"]
        exp_noappend = {'foo':['a','b','c','d','e','f','g'],
                        'bar':['h','i','j','k','l','m','n'],
                        'foobar':['h','i','j',None,'l',None, None]}
        exp_append = {'foo':['k__a','p__b','c__c','o__d','f__e','g__f','s__g'],
                      'bar':['k__h','p__i','c__j','o__k','f__l','g__m','s__n'],
                      'foobar':['k__h','p__i','c__j','o__','f__l','g__','s__']}
        obs_noappend = load_consensus_map(input, False)
        obs_append = load_consensus_map(input, True)
        self.assertEqual(obs_noappend, exp_noappend)
        self.assertEqual(obs_append, exp_append)

    def test_collect_names_at_ranks_counts(self):
        """correctly returns total counts for names at ranks"""
        input = "((a,b)c,(d,(e,f)g)h,(i,j)k)l;"
        tipname_map = {'a':['1','2','3','4','5','6','7'],
                       'b':['1','2','3',None,'5','6','8'],
                       'd':['1','2','3','a','5','6','9'],
                       'e':['1','2','3',None,'5','6','9'],
                       'i':['1','2','3','a','5','6','9'],
                       'j':['1','2','3','4','5','6','9']}
        tree = load_tree(input, tipname_map)
        
        exp = {0:{'1':6}, 
               1:{'2':6}, 
               2:{'3':6}, 
               3:{'4':2, 'a':2}, 
               4:{'5':6},
               5:{'6':6},
               6:{'7':1,'8':1,'9':4}}

        obs = collect_names_at_ranks_counts(tree)
        self.assertEqual(obs, exp)

    def test_load_tree(self):
        """correctly loads and decorates tiplook info on a tree"""
        input = "((a,b)c,(d,(e,f)g)h,(i,j)k)l;"
        tipname_map = {'a':['1','2','3','4','5','6','7'],
                       'b':['1','2','3','4','5','6','8'],
                       'd':['1','2','3','4','5','6','9'],
                       'e':['1','2','3','4','5','6','10'],
                       'i':['1','2','3','4','5','6','12'],
                       'j':['1','2','3','4','5','6','13']}

        # exp in Name: (tipstart, tipstop, consensus)
        exp = {'a':(0,0,['1','2','3','4','5','6','7']),
               'b':(1,1,['1','2','3','4','5','6','8']),
               'c':(0,1,[None] * 7), 
               'd':(2,2,['1','2','3','4','5','6','9']), 
               'e':(3,3,['1','2','3','4','5','6','10']), 
               'f':(4,4,[None] * 7),
               'g':(3,4,[None] * 7), 
               'h':(2,4,[None] * 7), 
               'i':(5,5,['1','2','3','4','5','6','12']), 
               'j':(6,6,['1','2','3','4','5','6','13']), 
               'k':(5,6,[None] * 7),
               'l':(0,6,[None] * 7)}

        obstree = load_tree(input, tipname_map, verbose=False)
        obs = {}
        for node in obstree.traverse(include_self=True):
            obs[node.Name] = (node.TipStart, node.TipStop, node.Consensus)

        self.assertEqual(obs, exp)

    def test_decorate_name_relative_freqs(self):
        """correctly decorate relative frequency information on a tree"""
        input = "((a,b)c,(d,(e,f)g)h,(i,j)k)l;"
        tipname_map = {'a':['1','2','3','4','5','6','7'],
                       'b':['1','2','3','4','5','6','8'],
                       'd':['1','2','3','4','5','6','8'],
                       'e':['1','2','3','4','a','6','7'],
                       'i':['1','2','3','4','a',None,'7'],
                       'j':['1','2','3','4','a',None,'8']}

        total_counts = {0:{'1':10, 'foo':5},
                1:{'2':6},
                2:{'3':12},
                3:{'4':6,'bar':5},
                4:{'5':6,'a':3},
                5:{'6':6},
                6:{'7':3,'8':3}}

        tree = load_tree(input, tipname_map)
        decorate_ntips(tree)
        decorate_name_relative_freqs(tree, total_counts, 1)

        exp_root = {0:{'1':.6},
                    1:{'2':1.0},
                    2:{'3':.5},
                    3:{'4':1.0},
                    4:{'5':.5, 'a':1.0},
                    5:{'6':4.0/6},
                    6:{'7':1.0,'8':1.0}}

        self.assertFloatEqual(tree.ConsensusRelFreq, exp_root)

    def test_set_ranksafe(self):
        """correctly set ranksafe on tree"""
        input = "((a,b)c,(d,(e,f)g)h,(i,j)k)l;"
        tipname_map = {'a':['1','2','3','4','5','6','7'],
                       'b':['1','2','3','b','5','6','8'],
                       'd':['1','2','3','4','5','6','8'],
                       'e':['1','2','3','b','a','6','7'],
                       'i':['1','2','3','4','a','6','7'],
                       'j':['1','2','3','b','a','6','8']}

        total_counts = {0:{'1':10, 'foo':5},
                1:{'2':20},
                2:{'3':10},
                3:{'4':4,'b':5},
                4:{'5':7,'a':5},
                5:{'6':6},
                6:{'7':3,'8':3}}

        tree = load_tree(input, tipname_map)
        decorate_ntips(tree)
        decorate_name_relative_freqs(tree, total_counts, 1)
        set_ranksafe(tree)

        #exp_root = ['Yes','Majority','Yes','No','Yes','Yes','No']
        exp_root = [True,False,True,False,True,True,False]
        self.assertEqual(tree.RankSafe, exp_root)

    def test_name_node_score_fold(self):
        """fucking hate taxonomy"""
        input = "((a,b)c,(d,(e,f)g)h,(i,j)k)l;"
        tipname_map = {'a':['1','2','3','4','5','6','8'],
                       'b':['1','2','3','4','5','6','8'],
                       'd':['1','2','3','f','e','c','9'],
                       'e':['1','2','3','f','e','c','9'],
                       'i':['1','2','3','g','a','h','11'],
                       'j':['1','2','3','g','a','h','12']}

        total_counts = {0:{'1':10, 'foo':5},
                1:{'2':10},
                2:{'3':10},
                3:{'4':3,'f':5,'g':5},
                4:{'5':7,'a':5,'e':4},
                5:{'6':3,'c':3,'d':2,'h':3},
                6:{'8':3,'9':2,'10':2,'11':2,'12':2}}

        tree = load_tree(input,tipname_map)
        decorate_ntips(tree)
        decorate_name_relative_freqs(tree, total_counts, 1)
        set_ranksafe(tree)
        pick_names(tree)
        name_node_score_fold(tree)
        exp_root = ['1','2','3',None,None,None,None]
        expc0 = [None, None, None, '4',None,None,None]
        expc1 = [None, None, None, None,'e','c','9']
        expc2 = [None, None, None, None, None,'h',None]
        expc1c1 = [None] * 7

        self.assertEqual(tree.RankNames, exp_root)
        self.assertEqual(tree.Children[0].RankNames, expc0)
        self.assertEqual(tree.Children[1].RankNames, expc1)
        self.assertEqual(tree.Children[2].RankNames, expc2)
        self.assertEqual(tree.Children[1].Children[1].RankNames, expc1c1)

    def test_validate_all_paths(self):
        """complains correctly about badpaths"""
        input = "(((((1,2)s__,(3,4)s__)g__)p__),((5,6)f__)f__,((7,8)c__)o__);"
        t = load_tree(input, {})
        exp = [t.getNodeMatchingName('5'),
               t.getNodeMatchingName('6'),
               t.getNodeMatchingName('7'),
               t.getNodeMatchingName('8')]
        obs = validate_all_paths(t)
        self.assertEqual(obs, exp)

    def test_best_name_freqs_for_nodes(self):
        """correctly gets the frequencies per name per node"""
        input = "((a,b)c,(d,(e,f)g)h,(i,j)k)l;"
        tipname_map = {'a':['1','2','3','4','5','6','7'],
                       'b':['1','2','3','b','5','6','8'],
                       'd':['1','2','3','4','5','6','7'],
                       'e':['1','2','3','b','a','foo','7'],
                       'i':['1','2','3','4','a','foo','8'],
                       'j':['1','2','3','b','a','foo','8']}

        total_counts = {0:{'1':10, 'foo':5},
                1:{'2':10},
                2:{'3':10},
                3:{'4':4,'b':5},
                4:{'5':7,'a':5},
                5:{'6':3,'foo':3},
                6:{'7':3,'8':4}}

        tree = load_tree(input,tipname_map)
        decorate_ntips(tree)
        decorate_name_relative_freqs(tree, total_counts, 1)
        set_ranksafe(tree)
        pick_names(tree)


        #result = best_name_freqs_for_nodes(tree)

        cnode = tree.Children[0]
        hnode = tree.Children[1]
        knode = tree.Children[2]
        #exp = {0:{'1':[(0.6,tree)]}
        #        1:{'2':[(0.6,tree)]},
        #        2:{'3':[(0.6,tree)]},
        #        3:{'4':[(0.75,tree)],'b':[(0.6,tree)]},
        #        4:{'a':[(0.6, tree),[(0.5,knode)]

    def test_pick_names(self):
        """correctly pick names to retain on a tree"""
        input = "((a,b)c,(d,(e,f)g)h,(i,j)k)l;"
        tipname_map = {'a':['1','2','3','4','5','6','7'],
                       'b':['1','2','3','b','5','6','8'],
                       'd':['1','2','3','4','5','6','7'],
                       'e':['1','2','3','b','a','foo','7'],
                       'i':['1','2','3','4','a','foo','8'],
                       'j':['1','2','3','b','a','foo','8']}

        total_counts = {0:{'1':10, 'foo':5},
                1:{'2':10},
                2:{'3':10},
                3:{'4':4,'b':5},
                4:{'5':7,'a':5},
                5:{'6':3,'foo':3},
                6:{'7':3,'8':4}}

        tree = load_tree(input,tipname_map)
        decorate_ntips(tree)
        decorate_name_relative_freqs(tree, total_counts, 1)
        set_ranksafe(tree)
        pick_names(tree)
        exp_root = ['1','2','3',None,None,None,None]
        self.assertEqual(tree.RankNames, exp_root)

        expc0 = [None,None,None,None,None,'6',None]
        expc1 = [None,None,None,None,None,None,'7']
        expc2 = [None,None,None,None,None,'foo','8']
        expc1c1 = [None] * 7

        self.assertEqual(tree.Children[0].RankNames, expc0)
        self.assertEqual(tree.Children[1].RankNames, expc1)
        self.assertEqual(tree.Children[2].RankNames, expc2)
        self.assertEqual(tree.Children[1].Children[1].RankNames, expc1c1)

    def test_walk_consensus_tree(self):
        """correctly walk consensus tree"""
        input = [['a','b','c','d','e','f','g'],
                 ['a','b','c',None,None,'x','y'],
                 ['h','i','j','k','l','m','n'],
                 ['h','i','j','k','l','m','q'],
                 ['h','i','j','k','l','m','n']]
        root, lookup = make_consensus_tree(input, check_for_rank=False)
        exp1 = ['n','m','l']
        exp2 = ['a']
        exp3 = ['x','f__','o__','c']
        obs1 = walk_consensus_tree(lookup, 'n', 3, reverse=False)
        obs2 = walk_consensus_tree(lookup, 'a', 3, reverse=False)
        obs3 = walk_consensus_tree(lookup, 'x', 4, reverse=False)
        self.assertEqual(obs1, exp1)
        self.assertEqual(obs2, exp2)
        self.assertEqual(obs3, exp3)

    def test_make_consensus_tree(self):
        """correctly gets parent consensus counts"""
        input = [['a','b','c','d','e','f','g'],
                 ['a','b','c',None,None,'x','y'],
                 ['h','i','j','k','l','m','n'],
                 ['h','i','j','k','l','m','q'],
                 ['h','i','j','k','l','m','n']]
        exp_str = "(((((((g)f)e)d,(((y)x)))c)b)a,((((((n,q)m)l)k)j)i)h);"
        
        obs_root, lookup = make_consensus_tree(input, check_for_rank=False)

        self.assertEqual(obs_root.getNewick(with_distances=False), exp_str)
        self.assertNotContains(None, lookup)

    def test_decorate_ntips(self):
        """correctly decorate the tree with the NumTips param"""
        input = "(((a,b)c,(d,e,f)g)h,(i,j)k)l;"
        tree = DndParser(input)
        tips = dict([(tip.Name, tip) for tip in tree.tips()])
        tips['a'].Consensus = [1,2,3,4,5,6,7]
        tips['b'].Consensus = [None,None,None,5,None,None,None]
        tips['d'].Consensus = [1,2,3,4,5,6,8]
        tips['e'].Consensus = [None, None,None,None,None,None,None]
        tips['f'].Consensus = [1,2,3,4,5,6,8]
        tips['i'].Consensus = [1,2,3,4,5,6,8]
        tips['j'].Consensus = [1,2,3,4,5,6,8]
        decorate_ntips(tree)
        self.assertEqual(tree.NumTips, 6)
        self.assertEqual(tree.Children[0].NumTips, 4)
        self.assertEqual(tree.Children[1].NumTips, 2)
        self.assertEqual(tree.Children[0].Children[0].NumTips, 2)
        self.assertEqual(tree.Children[0].Children[1].NumTips, 2)

    def test_get_nearest_named_ancestor(self):
        """correctly get the nearest named ancestor"""
        t = DndParser("(((s1,s2)g1,s3))root;")
        t2 = DndParser("(((s1,s2)g1,s3));")
        exp_t = t
        exp_t2 = None
        obs_t = get_nearest_named_ancestor(t.getNodeMatchingName('s3'))
        obs_t2 = get_nearest_named_ancestor(t2.getNodeMatchingName('s3'))
        self.assertEqual(obs_t, exp_t)
        self.assertEqual(obs_t2, exp_t2)

    def test_backfill_names_gap(self):
        """correctly backfill names"""
        consensus_tree = DndParser("(((s1,s2)g1,(s3,s4)g2,(s5,s6)g3)f1)o1;")
        rank_lookup = {'s':6,'g':5,'f':4,'o':3,'c':2,'p':1,'k':0}
        for n in consensus_tree.traverse(include_self=True):
            n.Rank = rank_lookup[n.Name[0]]
        input = "((((1)s1,(2)s2),((3)s3,(4)s5)))o1;"
        lookup = dict([(n.Name, n) for n in consensus_tree.traverse(include_self=True)])
        #exp = "((((1)s1,(2)s2)g1,((3)'g2; s3',(4)'g3; s5')))'o1; f1'"
        t = DndParser(input)
        t.Rank = 3
        t.Children[0].Rank = None
        t.Children[0].Children[0].Rank = None
        t.Children[0].Children[1].Rank = None
        t.Children[0].Children[0].Children[0].Rank = 6
        t.Children[0].Children[0].Children[1].Rank = 6
        t.Children[0].Children[1].Children[0].Rank = 6
        t.Children[0].Children[1].Children[1].Rank = 6

        backfill_names_gap(t, lookup)

        self.assertEqual(t.BackFillNames, ['o1'])
        self.assertEqual(t.Children[0].BackFillNames, [])
        self.assertEqual(t.Children[0].Children[0].BackFillNames, [])
        self.assertEqual(t.Children[0].Children[1].BackFillNames, [])
        self.assertEqual(t.Children[0].Children[0].Children[0].BackFillNames, ['f1','g1','s1'])
        self.assertEqual(t.Children[0].Children[0].Children[1].BackFillNames, ['f1','g1','s2'])
        self.assertEqual(t.Children[0].Children[1].Children[0].BackFillNames, ['f1','g2','s3'])
        self.assertEqual(t.Children[0].Children[1].Children[1].BackFillNames, ['f1','g3','s5'])

    def test_backfill_names_dangling(self):
        """correctly fill in dangling missing ranks"""
        consensus_tree = DndParser("(((s1,s2)g1,(s3,s4)g2,(s5,s6)g3)f1)o1;")
        input = "((((1),(2)),((3),(4))))'o1; f1';"
        lookup = dict([(n.Name, n) for n in consensus_tree.traverse(include_self=True)])
        #exp = "((((1),(2)),((3),(4))))'o1; f1';"

    def test_commonname_promotion(self):
        """correctly promote names if possible"""
        consensus_tree = DndParser("(((s1,s2)g1,(s3,s4)g2,(s5,s6)g3)f1)o1;")
        rank_lookup = {'s':6,'g':5,'f':4,'o':3,'c':2,'p':1,'k':0}
        for n in consensus_tree.traverse(include_self=True):
            n.Rank = rank_lookup[n.Name[0]]
        input = "((((1)s1,(2)s2),((3)s3,(4)s5)))o1;"
        lookup = dict([(n.Name, n) for n in consensus_tree.traverse(include_self=True)])
        exp = "((((1)s1,(2)s2)g1,((3)'g2; s3',(4)'g3; s5')))'o1; f1';"
        t = DndParser(input)
        t.Rank = 3
        t.Children[0].Rank = None
        t.Children[0].Children[0].Rank = None
        t.Children[0].Children[1].Rank = None
        t.Children[0].Children[0].Children[0].Rank = 6
        t.Children[0].Children[0].Children[1].Rank = 6
        t.Children[0].Children[1].Children[0].Rank = 6
        t.Children[0].Children[1].Children[1].Rank = 6
        backfill_names_gap(t, lookup)
        commonname_promotion(t)
        self.assertEqual(t.getNewick(with_distances=False), exp)

if __name__ == '__main__':

    main()
