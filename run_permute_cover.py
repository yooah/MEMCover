#!/usr/bin/env python

########################################################################
# permute the given mutation profile
#
# python run_permute_cover.py fid pnum ptype -mut mut_file
#
# e.g.  to read "mut_data.txt" and run 100 * |E| edge swaps
# 	to create a permuted mutational profile
#   using TR (type restricted) method
#	and create "tr_permuted_cover_1.txt"
#
# >> python run_permute_cover.py 1 100  tr -mut=mut_data.txt
#
# use "to" instead of "tr" to run TO_permutation (type oblivious)
#   and create "to_permuted_cover_1.txt"
# >> python run_permute_cover.py 1 100 to -mut=mut_data.txt
#
# modify INPUT/OUTPUT FILES as well as config.py file to change options.
# see README.txt for input file format
########################################################################

import argparse
import logging

import config
from config import cancers, nsamples
from mut_ex import permute_mut_data
import utils.io as io

# read arguments
parser = argparse.ArgumentParser()
parser.add_argument("fid", help="permutation file id", type=int)
parser.add_argument("pnum", help="number of edge swaps", type=int)
parser.add_argument("ptype", help="permutation type (tr or to)", type=str)
parser.add_argument("-mut", "--mut_file", help="mutation file name", type=str)

args = parser.parse_args()
fid, pnum = args.fid, args.pnum,

logging.info("%d-th permutation.. (permute %d times)" % (fid, pnum))


if args.ptype == "to":
	prefix = "to_"     # TO permutation
	logging.info("permute all types...")
elif args.ptype == "tr":
	prefix = "tr_"        # TR permutation
	logging.info("permute within types...")

# INPUT_FILES
if args.mut_file is None:
	mut_file = config.mut_file # use default
else:
	mut_file = config.data_dir + args.mut_file

# OUTPUT_FILES
output_file = config.permute_dir+prefix+"permuted_cover_"+str(fid)+".txt"

# read data
(genes, samples, mut_dic) = io.read_mut_matrix(mut_file)
logging.info("\n****** reading data files...\n")

# cancer types and permutation methods
if args.ptype == "to":  # for TO permutation
	type_idx_dic = dict([("all", range(len(samples)))])
	cancers = ["all"]
elif args.ptype == "tr":  # for TR permutation
	# read subtype data (sample -> type)
	sample_type_dic = io.read_dic(config.subtype_file)
	# construct dictionary each cancer type mapped to a list of sample names
	type_idx_dic = {}
	for cancer in cancers:
		type_idx_dic[cancer] = list(filter(lambda i: sample_type_dic[samples[i]] == cancer, range(len(samples))))

# construct a bipartite graph for each cancer type
mut_graphs = permute_mut_data.construct_mut_graph_per_type(mut_dic, cancers, type_idx_dic)

# permute each bipartite graph separately
logging.info("\n****** permuting the graph...\n")
permuted_graphs = {}
for cancer in cancers:
	subG = mut_graphs[cancer]
	permuted_graphs[cancer] = permute_mut_data.permute_mut_graph(subG, genes, type_idx_dic[cancer], pnum)

# merge muted graphs and construct the permuted dic
logging.info("construct mut_dic from permuted graphs\n")
permuted_mut_dic = permute_mut_data.construct_mut_dic_from_graphs(permuted_graphs, genes, nsamples)

# write the coverage of genes in permutation alteration in compact format (unweighted)
io.write_mut_list(permuted_mut_dic, output_file)
