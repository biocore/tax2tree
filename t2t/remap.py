#!/usr/bin/env python

def parse_otu_map(lines):
    """Returns {rep: [members]}, members include rep"""
    res = {}
    for l in lines:
        fields = l.strip().split('\t')
        c_id = fields[0]
        rep = fields[1]
        members = fields[1:]
        res[rep] = members
    return res

