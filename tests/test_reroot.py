#!/usr/bin/env python

from cogent.parse.tree import DndParser
from t2t.reroot import reroot
from cogent.util.unit_test import TestCase, main

class RerootTests(TestCase):
    def setUp(self):
        pass

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

