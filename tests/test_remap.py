#!/usr/bin/env python

from t2t.remap import parse_otu_map, members_to_rep, remap_taxonomy
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

    def test_members_to_rep(self):
        otus = {"2": ["2", "3", "4"],
                "5": ["5", "6"],
                "7": ["7"]}
        exp = {'2': '2', '3': '2', '4': '2', '5': '5', '6': '5',
               '7': '7'}
        obs = members_to_rep(otus)
        self.assertEqual(obs, exp)

    def test_remap_taxonomy(self):
        mapping= {"2": ["2", "3", "4"],
                  "5": ["5", "6"],
                  "7": ["7"]}

        tax = {"3": ["a", "b", "c"],
               "7": ["x", "y", "z"],
               "8": ['foo', 'bar']}

        exp = {"2": ["a", "b", "c"],
               "3": ["a", "b", "c"],
               "4": ["a", "b", "c"],
               "7": ["x", "y", "z"],
               "8": ['foo', 'bar']}

        obs = remap_taxonomy(mapping, tax)
        self.assertEqual(obs, exp)

if __name__ == '__main__':
    main()
