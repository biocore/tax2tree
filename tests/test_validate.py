#!/usr/bin/env python

from unittest import TestCase, main

from skbio import TreeNode

from t2t.validate import (check_parse, check_n_levels, check_gap,
                          check_prefixes, ParseError, cache_tipnames,
                          get_polyphyletic, find_gap)


class VerifyTaxonomy(TestCase):
    def setUp(self):
        pass

    def test_check_parse(self):
        """returns valid parsed or raises"""
        exp = ("1", ("k__a", "p__b", "c__c", "o__d", "f__e", "g__f", "s__g"))
        obs = check_parse(good_string)
        self.assertEqual(obs, exp)

        self.assertRaises(ParseError, check_parse, bad_string)

    def test_check_n_levels(self):
        """requires N levels, or unclassified"""
        id_, parsed = check_parse(good_string)
        self.assertTrue(check_n_levels(parsed, 7))

        self.assertFalse(check_n_levels(parsed, 8))

        id_, parsed = check_parse(good_unclassified)
        self.assertTrue(check_n_levels(parsed, 7))

        id_, parsed = check_parse(bad_unclassified1)
        self.assertFalse(check_n_levels(parsed, 7))

        id_, parsed = check_parse(bad_unclassified2)
        self.assertFalse(check_n_levels(parsed, 7))

        id_, parsed = check_parse(bad_unclassified3)
        self.assertFalse(check_n_levels(parsed, 7))

        id_, parsed = check_parse(bad_unclassified4)
        self.assertFalse(check_n_levels(parsed, 7))

    def test_check_gap(self):
        """check if a gap exists in a string"""
        id_, parsed = check_parse(good_string)
        self.assertTrue(check_gap(parsed))

        id_, parsed = check_parse(good_trailing)
        self.assertTrue(check_gap(parsed))

        id_, parsed = check_parse(gap)
        self.assertFalse(check_gap(parsed))

    def test_find_gap(self):
        good_string_idx = -1
        gap_idx = 2
        trailing_idx = -1

        id_, parsed = check_parse(good_string)
        self.assertEqual(find_gap(parsed), good_string_idx)

        id_, parsed = check_parse(good_trailing)
        self.assertEqual(find_gap(parsed), trailing_idx)

        id_, parsed = check_parse(gap)
        self.assertEqual(find_gap(parsed), gap_idx)

    def test_check_prefixes(self):
        """Verify the expected prefixes are present"""
        prefixes = ['k', 'p', 'c', 'o', 'f', 'g', 's']
        id_, parsed = check_parse(good_string)
        self.assertTrue(check_prefixes(parsed, prefixes))

        id_, parsed = check_parse(bad_prefix)
        self.assertFalse(check_prefixes(parsed, prefixes))

    def test_cache_tipnames(self):
        """caches tipnames"""
        t = TreeNode.from_newick("((a,b)c,(d,e)f)g;")
        cache_tipnames(t)

        self.assertEqual(t.tip_names, ['a', 'b', 'd', 'e'])
        self.assertEqual(t.children[0].tip_names, ['a', 'b'])
        self.assertEqual(t.children[1].tip_names, ['d', 'e'])

    def test_get_polyphyletic(self):
        """get polyphyletic groups"""
        cons = {'a': ['K', 'X1', 'X'],
                'b': ['K', 'X1', 'X'],
                'd': ['K', 'X1', 'Y'],
                'e': ['K', 'X1', 'Y'],
                'g': ['K', 'X2', 'Z'],
                'h': ['K', 'X2', 'Z'],
                'i': ['K', 'X2', 'X'],
                'j': ['K', 'X2', 'X']}

        obs_poly = get_polyphyletic(cons)

        self.assertEqual(len(obs_poly), 6)
        self.assertEqual(sorted(obs_poly[('X', 2)].keys()), ['X1', 'X2'])
        self.assertEqual(sorted(obs_poly[('Y', 2)].keys()), ['X1'])
        self.assertEqual(sorted(obs_poly[('Z', 2)].keys()), ['X2'])
        self.assertEqual(sorted(obs_poly[('X1', 1)].keys()), ['K'])
        self.assertEqual(sorted(obs_poly[('X2', 1)].keys()), ['K'])
        self.assertEqual(sorted(obs_poly[('K', 0)].keys()), [None])

good_string = "1	k__a; p__b; c__c; o__d; f__e; g__f; s__g"

# space instead of tab
bad_string = "2 k__a; p__b; c__c; o__d; f__e; g__f; s__g"

gap = "20	k__a; p__b; c__; o__d; f__e; g__f; s__g"
good_trailing = "30	k__a; p__b; c__c; o__d; f__e; g__; s__"
good_unclassified = "30	k__; p__; c__; o__; f__; g__; s__"
bad_unclassified1 = "40	Unclassified"
bad_unclassified2 = "50	k__x; unclassified"
bad_unclassified3 = "60	k__x; "
bad_unclassified4 = "70	k__x"
bad_nlevels = "80	k__a; p__b; c__c"
bad_prefix = "1	k__a; p__b; c__c; q__d; f__e; g__f; s__g"

if __name__ == '__main__':
    main()
