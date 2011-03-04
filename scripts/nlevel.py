#!/usr/bin/env python

"""Run the nlevel workflow"""

from optparse import make_option
from cogent.util.misc import parse_command_line_parameters
from sys import exit, stdout

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2011"
__credits__ = ["Daniel McDonald"] 
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Production"

script_info = {}
script_info['brief_description'] = "Map a taxonomy onto a tree"
script_info['script_description'] = """Take an input taxonomy and attempt to place the names on the tree in a reasonable fashion."""

script_info['script_usage'] = []
script_info['script_usage'].append(("""Example:""","""python nlevel.py --consensus-map=taxonomy.txt --tree=tree.ntree --output=my_result.ntree"""))

script_info['output_description'] = """The output is both a tree and its associated taxonomy strings.'"""

script_info['required_options']=[make_option('--tree','-t',dest='tree',\
                    help='Input tree'),
              make_option('--consensus-map','-m',dest='consensus_map',\
                    help='Input consensus map'),
              make_option('--output','-o', dest='output',\
                    help='Output file name')]

script_info['optional_options']=[make_option('--verbose',dest='verbose',
                  action='store_true',default=False)]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if not opts.tree:
        print option_parser.usage()
        print
        print "Tree must be specified!"
        exit(1)

    if not opts.consensus_map:
        print option_parser.usage()
        print
        print "Consensus map must be specified!"
        exit(1)

    if not opts.output:
        print option_parser.usage()
        print
        print "Output must be specified!"
        exit(1)

    append_rank = False
    tipname_map = load_consensus_map(open(opts.consensus_map),
                                     append_rank,
                                     opts.verbose)
    tree = load_tree(open(opts.tree), tipname_map, opts.verbose)
    counts = collect_names_at_ranks_counts(tree, opts.verbose)
    decorate_ntips(tree)
    decorate_name_relative_freqs(tree, counts, 3, opts.verbose)
    set_ranksafe(tree, opts.verbose)
    pick_names(tree, opts.verbose)
    name_node_score_fold(tree, verbose=opts.verbose)

    if opts.verbose:
        print "SCORE: ", score_tree(tree, opts.verbose)
    set_preliminary_name_and_rank(tree)

    contree, contree_lookup = make_consensus_tree(tipname_map.values())
    backfill_names_gap(tree, contree_lookup, opts.verbose)
    commonname_promotion(tree)
    make_names_unique(tree, opts.verbose)


    constrings = pull_consensus_strings(tree, opts.verbose)

    f = open(opts.output + '-consensus-strings', 'w')
    f.write('\n'.join(constrings))
    f.close()

    save_bootstraps(tree, opts.verbose)
    f = open(opts.output, 'w')
    f.write(tree.getNewick(with_distances=True))
    f.close()

   
if __name__ == '__main__':
    main()
