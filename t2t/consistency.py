#!/usr/bin/env python

__author__ = "Donovan Park"
__copyright__ = "Copyright 2014, The tax2tree project"
__credits__ = ["Donovan Park"]
__license__ = "BSD"
__version__ = "1.0"
__maintainer__ = "Donovan Park"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

from collections import defaultdict

from numpy import mean


class Consistency(object):
    """Calculates the consistency of taxonomic groups within a reference tree.

    Consistency is a measure between 0 and 1 indicating how close a taxonomic
    group is to being monophyletic.
    """
    def __init__(self, taxa_counts, n_ranks):
        """Initialize class.

        Parameters
        ----------
        taxa_counts : dict of dict
            Maps taxonomic rank and taxon names to counts, [RANK][name] ->
            count
        n_ranks : Init
            Indicates number of ranks in taxonomic hierarchy
        """
        self.taxa_counts = taxa_counts
        self.n_ranks = n_ranks

    def calculate(self, tree, rooted):
        """Return taxonomic consistency of each taxa in tree.

        Consistency is at every node and the highest consistency
        reported for each taxa.

        Parameters
        ----------
        tree : TreeNode
        rooted : Boolean
            Indicates if tree should be treated as rooted
        """

        # determine total number of informative tips in tree for each rank
        total_informative_tips = defaultdict(int)
        for n in tree.tips():
            for rank in xrange(self.n_ranks):
                total_informative_tips[rank] += n.NumTipsRank[rank]

        # determine highest consistency node for each taxa
        consistency_index = {i: defaultdict(int) for i in xrange(self.n_ranks)}
        for n in tree.traverse(include_self=True):
            for rank in xrange(self.n_ranks):
                for name, total_taxa_cnt in self.taxa_counts[rank].iteritems():
                    node_taxa_count = n.TaxaCount[rank].get(name, 0)
                    incongruent_taxa = n.NumTipsRank[rank] - node_taxa_count

                    c = float(node_taxa_count) / (total_taxa_cnt +
                                                  incongruent_taxa)
                    if c > consistency_index[rank][name]:
                        consistency_index[rank][name] = c

                    if not rooted:
                        # consider consistency of taxa in other subtree since
                        # the tree is unrooted
                        node_taxa_count = total_taxa_cnt - node_taxa_count
                        incongruent_taxa = total_informative_tips[rank] - \
                            n.NumTipsRank[rank] - \
                            node_taxa_count
                        c = float(node_taxa_count) / (total_taxa_cnt +
                                                      incongruent_taxa)
                        if c > consistency_index[rank][name]:
                            consistency_index[rank][name] = c

        return consistency_index

    def write_taxon_consistency(self, output_file, consistency_index):
        """Write consistency of each taxon to file.

        Parameters
        ----------
        output_file : TreeNode
        consistency_index : dict of dict
          Taxonomic consistency returned by Consistency.calculate()
        """

        fout = open(output_file, 'w')
        fout.write('Taxon\tCount\tConsistency\n')
        for rank in xrange(self.n_ranks):
            for name, consistency in consistency_index[rank].iteritems():
                fout.write('%s\t%d\t%.3f\n' % (name,
                                               self.taxa_counts[rank][name],
                                               consistency))

        fout.close()

    def write_rank_consistency(self, output_file, consistency_index, min_taxa,
                               rank_order):
        """Write average consistency of each rank to file.

        Parameters
        ----------
        output_file : TreeNode
        consistency_index : dict of dict
          Taxonomic consistency returned by Consistency.calculate()
        min_taxa: minimum taxa in group for inclusion in calculated average
        """

        fout = open(output_file, 'w')
        fout.write('Rank #\tRank prefix\t# taxon\tAverage consistency\n')
        for rank in xrange(self.n_ranks):
            val = []
            for name, consistency in consistency_index[rank].iteritems():
                if self.taxa_counts[rank][name] >= min_taxa:
                    val.append(consistency)

            if len(val) > 0:
                fout.write('%s\t%s\t%d\t%.3f\n' % (rank, rank_order[rank],
                                                   len(val), mean(val)))
            else:
                fout.write('%s\t%s\t%d\t%s\n' % (rank, rank_order[rank],
                                                 len(val), 'NA'))

        fout.close()
