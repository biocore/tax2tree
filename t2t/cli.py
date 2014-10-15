from skbio import TreeNode

import t2t.nlevel as nl
import t2t.validate as val


def fetch(tree):
    t = TreeNode.from_newick(open(tree))
    ranks = set(nl.RANK_ORDER)
    res = []
    error = True

    failed_nodes = []
    for n in t.non_tips():
        if n.name and n.name[0] not in ranks:
            failed_nodes.append(n)

    if failed_nodes:
        for n in failed_nodes:
            tip = list(n.tips())[0]
            tip_id = tip.name
            path = [a.name for a in tip.ancestors() if a.name is not None]
            res.append("Unknown rank: %s" % n.name)
            res.append("\tA tip ID from the clade: %s" % tip_id)
            res.append("\tCurrent lineage in tree: %s" % '; '.join(path[::-1]))

    else:
        res = nl.pull_consensus_strings(t)
        error = False

    return res, error


def validate(lines, limit, flat_errors, hierarchy_errors):
    res = []
    if flat_errors:
        flat = val.flat_errors(lines)
        for err_type in sorted(flat):
            ids = ','.join(flat[err_type][:10])

            if len(flat[err_type]) > limit:
                ellipse = '...'
            else:
                ellipse = ''

            res.append(err_type)
            res.append('\t%s%s' % (ids, ellipse))

    if hierarchy_errors:
        hier = val.hierarchy_errors(lines)
        if hier:
            res.append("Multiple parents")
        for err in hier:
            res.append("\t%s" % err['Taxon'])
            for parent in err['Parents']:
                res.append("\t\t%s, %s" % (parent, err['Parents'][parent]))

    return res, False
