#!/usr/bin/env python

"""TODO: place all consensus methods here

load_consensus_map
make_consensus_tree
etc...
"""
from t2t.nlevel import RANK_ORDER
from numpy import zeros, where, logical_or, int_


def taxa_score(master, reps):
    """Score taxa strings by contradictions observed in reps"""
    n_ranks = len(RANK_ORDER)
    master_ids = frozenset(master.keys())

    master_order = master.keys()
    master_rows = {k: idx for idx, k in enumerate(master_order)}
    master_cons = [master[k] for k in master_order]

    scores = zeros((len(master_ids), n_ranks), dtype=float)

    n_reps = 0
    for rep in reps:
        n_reps += 1

        for id_, con in rep.items():
            if id_ not in master_ids:
                raise KeyError("Unknown key %s in replicate" % id_)

            row = master_rows[id_]

            for rank, name in enumerate(con):
                if name == master_cons[row][rank]:
                    scores[row, rank] += 1

        # missing taxa are not considered contradictions
        missing_taxa = master_ids - frozenset(rep)
        for k in missing_taxa:
            row = master_rows[k]
            scores[row] += 1

    scores /= n_reps

    # slice and dice the scores
    return {k: scores[master_rows[k]] for k in master}


def merge_taxa_strings_and_scores(master, scores):
    """Merge taxa strings and their scores, return {id_:(taxa,score)}"""
    return {k: list(zip(v, scores[k])) for k, v in master.items()}


def taxa_score_hash(master, reps):
    """Score each taxonomy string based on contradictions observed in reps"""
    n_ranks = len(RANK_ORDER)

    master_order = master.keys()
    scores = zeros((len(master_order), n_ranks), dtype=float)
    master_hash = hash_cons(master, master_order, n_ranks)

    n_reps = 0
    for rep in reps:
        n_reps += 1
        rep_hash = hash_cons(rep, master_order, n_ranks)
        # where the taxons are equal to the master
        # or if the the taxon is not present in the replicate
        scores += where(logical_or(rep_hash == master_hash, rep_hash == 0),
                        1, 0)
    scores /= n_reps

    # slice and dice the scores
    return {k: scores[idx] for idx, k in enumerate(master_order)}


def hash_cons(cons, order, n_ranks):
    """Returns a numpy array of hash values for the cons

    NOTE: expects that cons are always specified even if the taxon name does
    not exist. In other words, the following are acceptable for missing fieids:
    [None, 'None', k__, p__, etc...]. It is _NOT_ okay to use the empty string
    at a field. The python hash method returns 0 on an empty string, but never
    otherwise and this method treats 0 specially.
    """
    hashes = zeros((len(order), n_ranks), dtype=int_)

    for idx, id_ in enumerate(order):
        try:
            hashes[idx] = [hash(c) for c in cons[id_]]
        except KeyError:
            pass
        # defaults to zero if the consensus string isn't represented
    return hashes


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

    rank_names = None
    for co in cons:
        if co is None:
            continue
        rank_names = [c[0] for c in co if c is not None]
        break
    for idx, rank in enumerate(rank_names):
        # collect all cons that are classified (ie more info than k__)
        klassed = [c[idx].lower()
                   for c
                   in cons if c[idx] and c[idx][0] == rank]

        n_classified = len(klassed)

        rank_names = set(klassed)

        n_seqs[rank] = (n_classified, total_cons - n_classified)
        n_names[rank] = rank_names

    return (n_seqs, n_names)


def pretty_print_consensus_stats(stats):
    seqs, names = stats
    print('\t'.join(['rank', 'num_classified', 'num_unclassified',
                     'num_names']))
    for k in RANK_ORDER:
        print('\t'.join(map(str, [k, seqs[k][0], seqs[k][1], len(names[k])])))
