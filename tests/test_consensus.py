#!/usr/bin/env python

from t2t.consensus import get_consensus_stats, taxa_score, hash_cons, \
    taxa_score_hash, merge_taxa_strings_and_scores
from unittest import TestCase, main
from numpy import array, array_equal


class ConsensusTests(TestCase):

    def setUp(self):
        pass

    def test_merge_taxa_strings_and_scores(self):
        """merge taxa scores and the strings"""
        master = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c2', None, 'f__f2', 'g__g1', 's__s2'],
            'c': [None, None, None, None, None, None, None],
            'd': ['k__k2', 'p__x', 'c__c1', 'o__o10', 'f__f2', 'g__g2',
                  's__s3']}
        rep1 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'c': [None, None, None, None, None, None, None],
            'd': ['k__k2', '', 'c__c1', 'f__f2', 'f__f2', 'g__g2', 's__s3']}
        rep2 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c2', None, 'f__f2', 'g__g1', 's__s9'],
            'c': [None, None, None, None, None, None, None]}
        rep3 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c2', None, 'f__f2', 'g__g1', 's__s2'],
            'd': ['k__k2', '', 'c__c1', 'f__f2', 'f__f2', 'g__g2', 's__s3']}
        rep4 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c10', None, 'f__f2', 'g__g1', 's__s2'],
            'c': [None, None, None, None, None, None, None],
            'd': ['k__k1', '', 'c__c1', 'f__f2', 'f__f2', 'g__g2', 's__s3']}
        scores = taxa_score(master, [rep1, rep2, rep3, rep4])
        exp = {
            'a': [('k__k1', 1.0), ('p__p1', 1.0), ('c__c1', 1.0),
                  ('o__o1', 1.0), ('f__f1', 1.0), ('g__g1', 1.0),
                  ('s__s1', 1.0)],
            'b': [('k__k1', 1.0), ('p__p1', 1.0), ('c__c2', 0.75), (None, 1.0),
                  ('f__f2', 1.0), ('g__g1', 1.0), ('s__s2', 0.75)],
            'c': [(None, 1.0), (None, 1.0), (None, 1.0), (None, 1.0),
                  (None, 1.0), (None, 1.0), (None, 1.0)],
            'd': [('k__k2', 0.75), ('p__x', 0.25), ('c__c1', 1.0),
                  ('o__o10', 0.25), ('f__f2', 1.0), ('g__g2', 1.0),
                  ('s__s3', 1.0)]}
        obs = merge_taxa_strings_and_scores(master, scores)
        self.assertEqual(obs, exp)

    def test_taxa_score(self):
        """score taxa strings in the face of replicates"""
        master = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c2', None, 'f__f2', 'g__g1', 's__s2'],
            'c': [None, None, None, None, None, None, None],
            'd': ['k__k2', 'p__x', 'c__c1', 'o__o10', 'f__f2', 'g__g2',
                  's__s3']}
        rep1 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'c': [None, None, None, None, None, None, None],
            'd': ['k__k2', '', 'c__c1', 'f__f2', 'f__f2', 'g__g2', 's__s3']}
        rep2 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c2', None, 'f__f2', 'g__g1', 's__s9'],
            'c': [None, None, None, None, None, None, None]}
        rep3 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c2', None, 'f__f2', 'g__g1', 's__s2'],
            'd': ['k__k2', '', 'c__c1', 'f__f2', 'f__f2', 'g__g2', 's__s3']}
        rep4 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c10', None, 'f__f2', 'g__g1', 's__s2'],
            'c': [None, None, None, None, None, None, None],
            'd': ['k__k1', '', 'c__c1', 'f__f2', 'f__f2', 'g__g2', 's__s3']}
        # if malformed, as in rep3 rep4 for d, count as contridiction
        exp = {'a': array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
               'b': array([1.0, 1.0, 0.75, 1.0, 1.0, 1.0, 0.75]),
               'c': array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
               'd': array([0.75, 0.25, 1.0, 0.25, 1.0, 1.0, 1.0])}
        obs = taxa_score(master, [rep1, rep2, rep3, rep4])
        self.assertEqual(obs.keys(), exp.keys())
        for k in exp:
            self.assertTrue(array_equal(obs[k], exp[k]))

    def test_taxa_score_hash(self):
        """test hash based consensus scoring"""
        master = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c2', None, 'f__f2', 'g__g1', 's__s2'],
            'c': [None, None, None, None, None, None, None],
            'd': ['k__k2', 'p__x', 'c__c1', 'o__o10', 'f__f2', 'g__g2',
                  's__s3']}
        rep1 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'c': [None, None, None, None, None, None, None],
            'd': ['k__k2', 'p__', 'c__c1', 'f__f2', 'f__f2', 'g__g2', 's__s3']}
        rep2 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c2', None, 'f__f2', 'g__g1', 's__s9'],
            'c': [None, None, None, None, None, None, None]}
        rep3 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c2', None, 'f__f2', 'g__g1', 's__s2'],
            'd': ['k__k2', 'p__', 'c__c1', 'f__f2', 'f__f2', 'g__g2', 's__s3']}
        rep4 = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c10', None, 'f__f2', 'g__g1', 's__s2'],
            'c': [None, None, None, None, None, None, None],
            'd': ['k__k1', 'p__', 'c__c1', 'f__f2', 'f__f2', 'g__g2', 's__s3']}
        # if malformed, as in rep3 rep4 for d, count as contridiction
        exp = {'a': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
               'b': [1.0, 1.0, 0.75, 1.0, 1.0, 1.0, 0.75],
               'c': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
               'd': [0.75, 0.25, 1.0, 0.25, 1.0, 1.0, 1.0]}
        obs = taxa_score_hash(master, [rep1, rep2, rep3, rep4])
        self.assertEqual(obs.keys(), exp.keys())
        for k in exp:
            self.assertTrue(array_equal(obs[k], exp[k]))

    def test_hash_cons(self):
        """test turning consensus strings into hashes"""
        input = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c2', None, 'f__f2', 'g__g1', 's__s2'],
            'c': [None, None, None, None, None, None, None],
            'd': ['k__k2', 'p__x', 'c__c1', 'o__o10', 'f__f2', 'g__g2',
                  's__s3']}
        exp = array([map(hash, input['a']),
                     map(hash, input['d']),
                     map(hash, input['c']),
                     map(hash, input['b'])])
        obs = hash_cons(input, ['a', 'd', 'c', 'b'], 7)
        self.assertTrue(array_equal(obs, exp))

        exp = array([map(hash, input['a']),
                     map(hash, input['d']),
                     [0, 0, 0, 0, 0, 0, 0],
                     map(hash, input['c']),
                     map(hash, input['b'])])
        obs = hash_cons(input, ['a', 'd', 'e', 'c', 'b'], 7)
        self.assertTrue(array_equal(obs, exp))

    def test_get_consensus_stats(self):
        """Produces the correct stats"""
        input = {
            'a': ['k__k1', 'p__p1', 'c__c1', 'o__o1', 'f__f1', 'g__g1',
                  's__s1'],
            'b': ['k__k1', 'p__p1', 'c__c2', None, 'f__f2', 'g__g1', 's__s2'],
            'c': [None, None, None, None, None, None, None],
            'd': ['k__k2', '', 'c__c1', 'f__f2', 'f__f2', 'g__g2', 's__s3']}
        exp_nseqs = {
            'k': (3, 1), 'p': (2, 2), 'c': (3, 1), 'o': (1, 3), 'f': (3, 1),
            'g': (3, 1), 's': (3, 1)}
        exp_names = {'k': set(['k__k1', 'k__k2']), 'p': set(['p__p1']),
                     'c': set(['c__c1', 'c__c2']), 'o': set(['o__o1']),
                     'f': set(['f__f1', 'f__f2']), 'g': set(['g__g1',
                                                             'g__g2']),
                     's': set(['s__s1', 's__s2', 's__s3'])}
        #exp_names = {'k':2,'p':1,'c':2,'o':1,'f':2,'g':2,'s':3}
        obs_nseqs, obs_names = get_consensus_stats(input)
        self.assertEqual(obs_nseqs, exp_nseqs)
        self.assertEqual(obs_names, exp_names)

if __name__ == '__main__':
    main()
