#!/usr/bin/env python


def parse_otu_map(lines):
    """Returns {rep: [members]}, members include rep"""
    res = {}
    for l in lines:
        fields = l.strip().split('\t')
        rep = fields[1]
        members = fields[1:]
        res[rep] = members
    return res


def members_to_rep(otus):
    """Provides a lookup for a cluster member to its rep"""
    res = {}
    for rep, members in otus.iteritems():
        for member in members:
            res[member] = rep
    return res


def remap_taxonomy(mapping, taxa):
    """Remaps the taxonomy over the OTU clusters"""
    res = {}
    reps = members_to_rep(mapping)

    for tax_id, tax_str in taxa.iteritems():
        rep = reps.get(tax_id, None)

        if rep is None:
            res[tax_id] = tax_str
            continue

        for m in mapping[rep]:
            if m in res:
                print "%s was: %s\nnow is %s" % (m, mapping[m], tax_str)
            res[m] = tax_str
    return res
