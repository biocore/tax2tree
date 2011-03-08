#!/usr/bin/env python

from t2t.consensus import get_consensus_stats
from cogent.util.unit_test import TestCase, main

class ConsensusTests(TestCase):
    def setUp(self):
        pass

    def test_get_consensus_stats(self):
        """Produces the correct stats"""
        input = {'a':['k__k1','p__p1','c__c1','o__o1','f__f1','g__g1','s__s1'],
                 'b':['k__k1','p__p1','c__c2',None,'f__f2','g__g1','s__s2'],
                 'c':[None,None,None,None,None,None,None],
                 'd':['k__k2','','c__c1','f__f2','f__f2','g__g2','s__s3']}
        exp_nseqs = {'k':(3,1),'p':(2,2),'c':(3,1),'o':(1,3),'f':(3,1),
                     'g':(3,1),'s':(3,1)}
        exp_names = {'k':set(['k__k1','k__k2']),'p':set(['p__p1']),
                     'c':set(['c__c1','c__c2']),'o':set(['o__o1']),
                     'f':set(['f__f1','f__f2']),'g':set(['g__g1','g__g2']),
                     's':set(['s__s1','s__s2','s__s3'])}
        #exp_names = {'k':2,'p':1,'c':2,'o':1,'f':2,'g':2,'s':3}
        obs_nseqs, obs_names = get_consensus_stats(input)
        self.assertEqual(obs_nseqs, exp_nseqs)
        self.assertEqual(obs_names, exp_names)

if __name__ == '__main__':
    main()
