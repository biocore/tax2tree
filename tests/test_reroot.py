#!/usr/bin/env python

from skbio.core.tree import TreeNode
from t2t.reroot import reroot
from unittest import TestCase, main

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class RerootTests(TestCase):
    def setUp(self):
        pass

    def test_reroot(self):
        """Should correctly reroot a tree"""
        t = TreeNode.from_newick("(((a,b)c,(d,e)f)g,(h,i)j);")

        tips = ['a','b']
        for n in t.traverse():
            n.length = 1.0

        # note, g is lost because it has a single descendent and gets pruned off
        exp = "((a:1.0,b:1.0)c:0.5,((d:1.0,e:1.0)f:1.0,(h:1.0,i:1.0)j:2.0):0.5)root;"
        obs = reroot(t, tips)
        self.assertEqual(obs.to_newick(with_distances=True), exp)

if __name__ == '__main__':
    main()

