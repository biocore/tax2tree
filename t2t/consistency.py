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

class Consistency(object):
    def __init__(self, taxa_counts, n_ranks):
        self.taxa_counts = taxa_counts
        self.n_ranks = n_ranks
                         
    def calculate(self, tree, rooted):
        """Return taxonomic consistency of each taxa in tree.
        
        Parameters
        ----------
        tree : TreeNode
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
                for name, total_taxa_count in self.taxa_counts[rank].iteritems():  
                    node_taxa_count = n.TaxaCount[rank].get(name, 0)
                    incongruent_taxa = n.NumTipsRank[rank] - node_taxa_count
                    
                    c = float(node_taxa_count) / (total_taxa_count + incongruent_taxa)
                    if c > consistency_index[rank][name]:
                        consistency_index[rank][name] = c
                    
                    if not rooted:
                        # consider consistency of taxa in other subtree since the tree is unrooted
                        node_taxa_count = total_taxa_count - node_taxa_count
                        incongruent_taxa = total_informative_tips[rank] - n.NumTipsRank[rank] - node_taxa_count
                        c = float(node_taxa_count) / (total_taxa_count + incongruent_taxa)
                        if c > consistency_index[rank][name]:
                            consistency_index[rank][name] = c
                            
        return consistency_index

    def write(self, output_file, consistency_index):
        """Write consistency to file.
                
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
                fout.write('%s\t%d\t%.3f\n' % (name, self.taxa_counts[rank][name], consistency))
                
        fout.close()