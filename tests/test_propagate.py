import unittest
import skbio
from t2t.propagate import index_backbone, indexed_to_name


class PropagateTests(unittest.TestCase):
    def setUp(self):
        self.backbone = "((a,b)c,((d,e),(g,h)i)j)k;"

    def test_index_backbone(self):
        t = skbio.TreeNode.read([self.backbone])

        index_backbone(t)
        self.assertEqual(t.find('c').covered_tips,
                         frozenset(['a', 'b']))
        self.assertEqual(t.find('i').covered_tips,
                         frozenset(['g', 'h']))
        self.assertEqual(t.find('j').covered_tips,
                         frozenset(['d', 'e', 'g', 'h']))
        self.assertEqual(t.find('k').covered_tips,
                         frozenset(['a', 'b', 'd', 'e', 'g', 'h']))

        # the clade (d,e) does not have a named node so it should not
        # have names
        for n in t.traverse(include_self=True):
            if n.is_tip():
                self.assertEqual(n.covered_tips, frozenset())
            else:
                if n.name is None:
                    self.assertEqual(n.covered_tips, frozenset())
                else:
                    self.assertNotEqual(n.covered_tips, frozenset())

    def test_indexed_to_name(self):
        t = skbio.TreeNode.read([self.backbone])

        index_backbone(t)
        obs = indexed_to_name(t)
        exp = {frozenset(['a', 'b']): 'c',
               frozenset(['g', 'h']): 'i',
               frozenset(['d', 'e', 'g', 'h']): 'j',
               frozenset(['a', 'b', 'd', 'e', 'g', 'h']): 'k'}
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()
