#!/usr/bin/env python

"""consistency tests

While many of the tests utilize similar trees and input data, the overlap
is not necessarily 100%. Many of these inputs are written with specific tests
in mind.
"""

__author__ = "Donovan Park"
__copyright__ = "Copyright 2014, The tax2tree project"
__credits__ = ["Donovan Park"]
__license__ = "BSD"
__version__ = "1.0"
__maintainer__ = "Donovan Park"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

from unittest import TestCase, main

import t2t.nlevel as nl
from t2t.consistency import Consistency

from io import StringIO

class ConsistencyTests(TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        nl.set_rank_order(['d', 'p', 'c', 'o', 'f', 'g', 's'])

    def test_consistency_missing(self):
        """Test consistency of taxa in tree with missing taxa"""

        seed_con = 'f__Lachnospiraceae; g__Bacteroides; s__'
        nl.determine_rank_order(seed_con)
        tipname_map = {'a': ['f__Lachnospiraceae', 'g__Bacteroides', None],
                       'c': ['f__Lachnospiraceae', 'g__Bacteroides', 's__Bacteroides pectinophilus'],
                       'b': ['f__Lachnospiraceae', 'g__Bacteroides', None], 'e': [None, None, None],
                       'd': ['f__Lachnospiraceae', 'g__Bacteroides', 's__Bacteroides pectinophilus'],
                       'g': [None, None, None], 'f': ['f__Lachnospiraceae', 'g__Lachnospira', None],
                       'h': ['f__Lachnospiraceae', 'g__Lachnospira', 's__Bacteroides pectinophilus']}
        tree = nl.load_tree(StringIO(u'(((a,b),(c,d)),((e,f),(g,h)));'), tipname_map)

        counts = nl.collect_names_at_ranks_counts(tree)
        nl.decorate_ntips_rank(tree)
        nl.decorate_name_counts(tree)

        # determine taxonomic consistency of rooted tree
        #expected_consistency_index
        c = Consistency(counts, len(nl.RANK_ORDER))
        consistency_index = c.calculate(tree, rooted=True)

        self.assertAlmostEqual(consistency_index[0]['f__Lachnospiraceae'], 1.0)
        self.assertAlmostEqual(consistency_index[1]['g__Bacteroides'], 1.0)
        self.assertAlmostEqual(consistency_index[1]['g__Lachnospira'], 1.0)
        self.assertAlmostEqual(consistency_index[2]['s__Bacteroides pectinophilus'], 1.0)

        #determine consistency of unrooted tree
        consistency_index = c.calculate(tree, rooted=False)

        self.assertAlmostEqual(consistency_index[0]['f__Lachnospiraceae'], 1.0)
        self.assertAlmostEqual(consistency_index[1]['g__Bacteroides'], 1.0)
        self.assertAlmostEqual(consistency_index[1]['g__Lachnospira'], 1.0)
        self.assertAlmostEqual(consistency_index[2]['s__Bacteroides pectinophilus'], 1.0)

    def test_consistency_unrooted(self):
        """Test consistency of taxa with a taxa that is only monophyletic in unrooted tree"""

        seed_con = 'f__Lachnospiraceae; g__Bacteroides; s__'
        nl.determine_rank_order(seed_con)
        tipname_map = {'a': ['f__Lachnospiraceae', 'g__Bacteroides', 's__Bacteroides pectinophilus'],
                       'b': ['f__Lachnospiraceae', 'g__Bacteroides', 's__Bacteroides pectinophilus'],
                       'c': ['f__Lachnospiraceae', 'g__Bacteroides', 's__Bacteroides pectinophilus'],
                       'd': ['f__Lachnospiraceae', 'g__Bacteroides', 's__Bacteroides acidifaciens'],
                       'e': ['f__Lachnospiraceae', 'g__Bacteroides', 's__Bacteroides acidifaciens']}

        tree = nl.load_tree(StringIO(u'((a,b),(c,(d,e)));'), tipname_map)

        counts = nl.collect_names_at_ranks_counts(tree)
        nl.decorate_ntips_rank(tree)
        nl.decorate_name_counts(tree)

        # determine taxonomic consistency of rooted tree
        #expected_consistency_index
        c = Consistency(counts, len(nl.RANK_ORDER))
        consistency_index = c.calculate(tree, rooted=True)

        self.assertAlmostEqual(consistency_index[0]['f__Lachnospiraceae'], 1.0)
        self.assertAlmostEqual(consistency_index[1]['g__Bacteroides'], 1.0)
        self.assertAlmostEqual(consistency_index[2]['s__Bacteroides pectinophilus'], 0.66666666)
        self.assertAlmostEqual(consistency_index[2]['s__Bacteroides acidifaciens'], 1.0)

        #determine consistency of unrooted tree
        consistency_index = c.calculate(tree, rooted=False)

        self.assertAlmostEqual(consistency_index[0]['f__Lachnospiraceae'], 1.0)
        self.assertAlmostEqual(consistency_index[1]['g__Bacteroides'], 1.0)
        self.assertAlmostEqual(consistency_index[2]['s__Bacteroides pectinophilus'], 1.0)
        self.assertAlmostEqual(consistency_index[2]['s__Bacteroides acidifaciens'], 1.0)

if __name__ == '__main__':
    main()
