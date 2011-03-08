#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from t2t.util import combine_alignments, reroot
from cogent.parse.tree import DndParser

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class UtilTests(TestCase):
    def setUp(self):
        pass

    def test_combine_alignments(self):
        """Combine alignments, raise if intersecting ids"""
        lines1 = ['>a','AATTGGCC','>b','AATTAATT']
        lines2 = ['>c','AATTAGCC','>d','AATTGATT']
        exp = {'a':'AATTGGCC','b':'AATTAATT', 
               'c':'AATTAGCC','d':'AATTGATT'}
        obs = combine_alignments(lines1, lines2)
        self.assertEqual(obs, exp)

        lines1 = ['>a','AATTGGCC','>b','AATTAATT']
        lines2 = ['>a','AATTAACC','>C','AATTGATT']
        self.assertRaises(ValueError, combine_alignments, lines1, lines2)
        
    def test_reroot(self):
        """Should correctly reroot a tree"""
        t = DndParser("(((a,b)c,(d,e)f)g,(h,i)j);")
        tips = ['a','b']
        for n in t.traverse():
            n.Length = 1.0
        
        # note, g is lost because it has a single descendent and gets pruned off
        exp = "((a:1.0,b:1.0)c:0.5,((d:1.0,e:1.0)f:1.0,(h:1.0,i:1.0)j:2.0):0.5);"
        obs = reroot(t, tips)
        self.assertEqual(obs.getNewick(with_distances=True), exp)

if __name__ == '__main__':
    main()

