tax2tree
========

[![Build Status](https://travis-ci.org/biocore/tax2tree.png?branch=master)](https://travis-ci.org/biocore/tax2tree) [![Coverage Status](https://coveralls.io/repos/biocore/tax2tree/badge.png)](https://coveralls.io/r/biocore/tax2tree)

tax2tree assists in decorating an existing taxonomy onto a phylogenetic tree
with overlapping tip names.

An example test tree and input consensus set can be found in tests/data. Once
tax2tree is installed and paths are setup appropriately (see INSTALL), you can
run the following to reproduce the output seen in Figure 1 in 
[McDonald et al. 2011](http://www.ncbi.nlm.nih.gov/pubmed/22134646).

    $ t2t -t tests/data/test.ntree -m tests/data/tests.con -o test_output

Contact information
-------------------
Daniel McDonald (mcdonadt@colorado.edu)
