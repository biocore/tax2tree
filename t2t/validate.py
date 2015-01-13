#!/usr/bin/env python

from collections import defaultdict
from operator import add

from t2t.nlevel import determine_rank_order, make_consensus_tree, load_consensus_map
from t2t.util import unzip


class ParseError(Exception):
    pass


def check_parse(line):
    """Make sure the line can parse into (id_, [some,thing])"""
    try:
        id_, initial = line.strip().split('\t')
    except ValueError:
        raise ParseError("Unable to split in tab")

    parsed = tuple([n.strip() for n in initial.split(";")])

    if not len(parsed):
        raise ParseError("Line appears to not have a taxonomy")

    return (id_, parsed)


def check_n_levels(parsed, n):
    """Make sure there are n levels in parsed"""
    if len(parsed) == n:
        return True
    else:
        return False


def find_gap(parsed):
    """Find a gap else -1"""
    reduced = [p.split('__')[-1] for p in parsed]

    end_found = False
    end_idx = -1

    for idx, n in enumerate(reduced):
        if n and end_found:
            break

        if not n:
            end_found = True
            end_idx = idx

    return end_idx if end_idx != len(parsed) - 1 else -1


def check_gap(parsed):
    return find_gap(parsed) == -1


def check_prefixes(parsed, expected_prefixes):
    """Make sure each rank has the expected prefix"""
    for name, prefix in zip(parsed, expected_prefixes):
        try:
            obs, level_name = name.split('__', 1)
        except ValueError:
            return False

        if obs != prefix:
            return False

    return True


def cache_tipnames(t):
    """cache tipnames on the internal nodes"""
    for n in t.postorder(include_self=True):
        if n.is_tip():
            n.tip_names = [n.name]
        else:
            n.tip_names = reduce(add, [c.tip_names for c in n.children])


def get_polyphyletic(cons):
    """get polyphyletic groups and a representative tip"""
    tips, taxonstrings = unzip(cons.items())
    tree, lookup = make_consensus_tree(taxonstrings, False, tips=tips)
    cache_tipnames(tree)

    names = {}
    for n in tree.non_tips():
        if n.name is None:
            continue
        if (n.name, n.Rank) not in names:
            names[(n.name, n.Rank)] = {}
        if n.parent is not None:
            names[(n.name, n.Rank)][n.parent.name] = n.tip_names[0]

    return names


def hierarchy_errors(tax_lines):
    """Get errors in the taxonomy hierarchy"""
    conmap = load_consensus_map(tax_lines, False)
    names = get_polyphyletic(conmap)
    errors = []

    for (name, rank), parents in names.iteritems():
        if len(parents) > 1:
            err = {'Taxon': name, 'Rank': rank, 'Parents': parents}
            errors.append(err)

    return errors


def flat_errors(tax_lines):
    """Flat file errors"""
    inc_prefix = 'Incorrect prefixes'
    inc_nlevel = 'Incorrect number of levels'
    inc_gap = 'Gaps in taxonomy'

    seed_con = tax_lines[0].strip().split('\t')[1]
    rank_order = determine_rank_order(seed_con)

    nlevels = len(rank_order)
    errors = defaultdict(list)
    errors_seen = defaultdict(set)

    for line in tax_lines:
        id_, parsed = check_parse(line)

        if not check_prefixes(parsed, rank_order):
            if parsed not in errors_seen[inc_prefix]:
                errors_seen[inc_prefix].add(parsed)
                errors[inc_prefix].append(id_)

        if not check_n_levels(parsed, nlevels):
            if parsed not in errors_seen[inc_nlevel]:
                errors_seen[inc_nlevel].add(parsed)
                errors[inc_nlevel].append(id_)

        if not check_gap(parsed):
            gap_idx = find_gap(parsed)
            taxon_following_gap = gap_idx + 1

            # another +1 as the slice is exclusive
            if parsed[:taxon_following_gap + 1] not in errors_seen[inc_gap]:
                errors_seen[inc_gap].add(parsed[:taxon_following_gap + 1])
                errors['Gaps in taxonomy'].append(id_)

    return errors
