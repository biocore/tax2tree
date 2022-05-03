tax2tree
========

[![Build Status](https://github.com/biocore/tax2tree/actions/workflows/python-package-conda.yml/badge.svg) 

tax2tree assists in decorating an existing taxonomy onto a phylogenetic tree
with overlapping tip names.

An example test tree and input consensus set can be found in tests/data. Once
tax2tree is installed and paths are setup appropriately (see INSTALL), you can
run the following to reproduce the output seen in Figure 1 in 
[McDonald et al. 2011](http://www.ncbi.nlm.nih.gov/pubmed/22134646).

    $ t2t decorate -t tests/data/test_consistency.ntree -m tests/data/test_consistency.cons -o test_output

Additional methods
==================

tax2tree provides support for decorating a backbone within [.jplace](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009). 

    $ t2t decorate --placement <path> ...

When decorating a taxonomy where many lineages are represented by individual tips, you can add "nameholders" to the tree. This works by creating a node with zero branch length as the parent of each tip, providing an internal node for a given label. The classic tax2tree algorithm cannot place labels on tips, and this provides a workaround.

    $ t2t decorate --add-nameholder ...

Following multifurcation placement resolution, it is possible a single node, or multifurcation will have been created on the parent edge of a named internal node. That new internal node is reasonably likely to be of the same clade. To account for this, it is possible to perform a single node promotion to the new internal node. To avoid overpromotion, it is necessary to specify the set of tips which correspond to the placements so that multifurcations can be detected

    $ t2t promote-multifurcation ...

Contact information
-------------------
Daniel McDonald (mcdonadt@colorado.edu)
