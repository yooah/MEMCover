#!/usr/bin/env python

########################################################################
#
# >> run_module_cover k ewth
#
# * positional arguments:
#   k                     number of cover for each patient
#   ewth                  edge weight threshold
#
#
# e.g.,
# to run module cover with k=15 (coverage), edge weight threshold=0.2 and other optional arguments
# >> python run_module_cover.py 15 0.2 -mut=mut_data.txt -net=human_net.net -out=combined
#
# * optional arguments:
#   -h, --help            show this help message and exit
#   -mut MUT_FILE, --mut_file MUT_FILE
#                         mutation file name (default in config file)
#   -net NET_FILE, --net_file NET_FILE
#                         network file name (default in config file)
#   -mw MUTW, --mutw MUTW
#                         somatic mutation rate weight relative to cnv
#   -nw NW_FILE, --nw_file NW_FILE
#                         node weight file name
#   -out OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
#                        specify output filename prefix
#
# See README file for input file formats.
#
# * results will be stored in
#   "results/combined_itr_15_0.2.txt": each gene in a row in the order added
#   "results/combined_modules_15_0.2.txt": one module in a row
#
#
# For the details of the algorithms and parameters, see
# Yoo-Ah Kim, Dong-Yeon Cho, Phuong Dao, and Teresa M. Przytycka
# MEMCover: integrated analysis of mutual exclusivity and functional network
# reveals dysregulated pathways across multiple cancer types
# Bioinformatics 2015 31: i284-i292
#
#
########################################################################

# import packages
import argparse
import logging
import networkx as nx

import module_cover.module_cover as module_cover
import utils.io as io
import config

# Read arguments
parser = argparse.ArgumentParser()
parser.add_argument("k", help="number of cover for each patient", type=int)
parser.add_argument("ewth", help="edge weight threshold", type=float)
parser.add_argument("-mut", "--mut_file", help="mutation file name", type=str)
parser.add_argument("-net", "--net_file", help="network file name", type=str)
parser.add_argument("-mw", "--mutw", help="mut weight relative to cnv", type=str)
parser.add_argument("-nw", "--nw_file", help="node weight file name", type=str)
parser.add_argument("-out", "--output_prefix", help="specify output filename prefix", type=str)

args = parser.parse_args()

# essential arguments
k, ewth = args.k, args.ewth

# optional arguments
if args.mutw is None:
    mutw = 1
else:
    mutw = args.mutw

# INPUT FILES
if args.mut_file is None:
    mut_file = config.mut_file # use default
else:
    mut_file = config.data_dir + args.mut_file
if args.net_file is None:
    net_file = config.hn_file # use default
else:
    net_file = config.data_dir + args.net_file

nw_file = args.nw_file

# OUTPUT FILES
mod_dir = "results/"
logging.debug("... results will be stored in "+mod_dir+" directory....\n")
if args.output_prefix is None:
    output = mod_dir + "module_cover"
else:
    output = mod_dir + args.output_prefix

# read mutation file
logging.info("... read "+mut_file+"....\n")
genes, samples, mut_dic = io.read_mut_matrix(mut_file, mutw, nw_file)

# read network file
logging.info("... read graph and edge weight....\n")
G = nx.read_edgelist(net_file, nodetype=str, data=(('weight',float),))
# create score dictionary
temp_dic =nx.get_edge_attributes(G, 'weight')
score_dic = {}
for e in temp_dic:  # make sure to have entries for both direction
    if e[0] not in score_dic:
        score_dic[e[0]] = {}
    score_dic[e[0]][e[1]] = temp_dic[e]
    if e[1] not in score_dic:
        score_dic[e[1]] = {}
    score_dic[e[1]][e[0]] = temp_dic[e]


# run module cover
logging.info("... output will be stored in "+output+"....\n")
# results in each iteration of module cover will be written in this file
results_output=output + "_itr_"+str(k)+"_"+str(ewth)+".txt"

M, total_cost = module_cover.greedy_module_cover(mut_dic, G, k, score_dic, ewth, results_output, True)

# write each module in a row
module_file = open(output + "_modules_"+str(k)+"_"+str(ewth)+".txt", 'w')
for i in range(len(M)):
    module_file.write("%d\t%s\n" %(i, ",".join(M[i])))

module_file.close()


