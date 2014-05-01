#!/usr/bin/env python

from t2t.remap import parse_otu_map
from unittest import TestCase, main


class RemapTests(TestCase):
    def test_parse_otu_map(self):
        mapping = ("1\t2\t3\t4\n"
                   "2\t5\t6\n"
                   "3\t7")
        exp = {"2": ["2", "3", "4"],
               "5": ["5", "6"],
               "7": ["7"]}
        obs = parse_otu_map(mapping.splitlines())
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
