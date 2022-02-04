import unittest
import skbio
from t2t.transfer import index_backbone, indexed_to_name


class PropagateTests(unittest.TestCase):
    def setUp(self):
        self.backbone = "((a,b)c,((d,e),(g,h)i)j)k;"
        self.with_placements = "((a,(b,X1)),(((X2,d),e),(g,h),X3),X4)k;"
        self.inconsistent = "((e,(b,X1)),(((X2,d),a),(g,h),X3),X4)k;"

    def test_transfer(self):
        b = skbio.TreeNode.read([self.backbone])
        wp = bp.parse_newick(self.with_placements)
        res = transfer(b, wp)
        res = bp.to_skbio_treenode(res)
        self.assertEqual(res.find('a').parent.name, 'c')
        self.assertEqual(res.find('b').parent.parent.name, 'c')
        self.assertEqual(res.find('X1').parent.name, 'c')
        self.assertEqual(res.find('d').parent.name, None)
        self.assertEqual(res.find('e').parent.name, None)
        self.assertEqual(res.find('g').parent.name, 'i')
        self.assertEqual(res.find('h').parent.name, 'i')
        self.assertEqual(res.find('X3').parent.name, 'j')
        self.assertEqual(res.find('X4').parent.name, 'k')

    def test_transfer_inconsistent(self):
        b = skbio.TreeNode.read([self.backbone])
        wp = bp.parse_newick(self.inconsistent)
        with self.assertRaises(KeyError):
            transfer(b, wp)

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
