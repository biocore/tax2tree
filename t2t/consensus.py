#!/usr/bin/env python

"""TODO: place all consensus methods here

load_consensus_map
make_consensus_tree
etc...
"""
from nlevel import RANK_ORDER

def get_consensus_stats(consensus_map):
    """Returns consensus stats, expects rank prefix

    Returns a tuple of two dicts: 
    
    - sequence counts per level, (classified, unclassified)
    - contains name counts per level
    """
    n_seqs = {}
    n_names = {}

    cons = consensus_map.values()
    total_cons = len(cons)

    rank_names = [c[0] for c in cons[0]]
    for idx, rank in enumerate(rank_names):
        # collect all cons that are classified (ie more info than k__)
        klassed = [c[idx].lower() for c in cons if c[idx] and c[idx][0]==rank]

        n_classified = len(klassed)

        rank_names = set(klassed)

        n_seqs[rank] = (n_classified, total_cons - n_classified)
        n_names[rank] = rank_names

    return (n_seqs, n_names)

def pretty_print_consensus_stats(stats):
    seqs, names = stats
    print '\t'.join(['rank','num_classified','num_unclassified','num_names'])
    for k in RANK_ORDER:
        print '\t'.join(map(str, [k,seqs[k][0],seqs[k][1],len(names[k])]))
