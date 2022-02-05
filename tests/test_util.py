#!/usr/bin/env python

from unittest import TestCase, main
from t2t.util import reroot, unzip
from skbio import TreeNode

import sys

if sys.version_info[0] == 2:
    from StringIO import StringIO
else:
    from io import StringIO

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The tax2tree project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"


class UtilTests(TestCase):

    def setUp(self):
        pass

    def test_reroot(self):
        """Should correctly reroot a tree"""
        t = TreeNode.read(StringIO(u"(((a,b)c,(d,e)f)g,(h,i)j);"))
        tips = ['a', 'b']
        for n in t.traverse():
            n.Length = 1.0

        # note, g is lost because it has a single descendent and gets pruned
        # off
        exp = ("((a:1.0,b:1.0)c:0.5,((d:1.0,e:1.0)f:1.0,(h:1.0,i:1.0)"
               "j:2.0):0.5);")
        exp = "((a,b)c,((d,e)f,(h,i)j):0.0);"
        obs = reroot(t, tips)

        self.assertEqual(str(obs).rstrip(), exp)

    def test_reroot_edges(self):
        t = TreeNode.read(StringIO(u"(((a,b)c,(d,e)f)g,(h,i)j);"))
        tips = ['a', 'b']
        for n in t.traverse():
            n.Length = 1.0
        t.find('a').edge_num = 1
        t.find('b').edge_num = 2
        t.find('c').edge_num = 3
        t.find('f').edge_num = 4
        obs = reroot(t, tips)
        self.assertEqual(obs.find('a').edge_num, 1)
        self.assertEqual(obs.find('b').edge_num, 2)
        self.assertEqual(obs.find('c').edge_num, 3)
        self.assertEqual(obs.find('f').edge_num, 4)

    def test_unzip(self):
        """unzip(items) should be the inverse of zip(*items)"""
        chars = [list('abcde'), list('ghijk')]
        numbers = [[1, 2, 3, 4, 5], [0, 0, 0, 0, 0]]
        strings = [["abcde", "fghij", "klmno"], ['xxxxx'] * 3]

        lists = [chars, numbers, strings]
        zipped = [zip(*i) for i in lists]
        unzipped = [unzip(i) for i in zipped]

        for u, l in zip(unzipped, lists):
            self.assertEqual(u, l)

if __name__ == '__main__':
    main()
