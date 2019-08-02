###### post processing of module cover 
# python postproc_module_cover.py 

import sys
import os
import argparse
import networkx as nx
import pandas
import importlib

from config import *
import module_cover.module_cover2 as module_cover2
import utils.io as io
importlib.reload(io)

# parser = argparse.ArgumentParser()
# parser.add_argument("pth", help="edge weight percent threhold", type=int)
# parser.add_argument("mutw", help="mut weight relative to cnv ", type=int)
# args = parser.parse_args()
# pth = args.pth

# arguments
wth = 0.18540455679 # weight threhold corresponding to the top 40%
weight_type = "norm_combined"

# input filenames
module_file = "results/norm_modules.txt"
scoresfile = "data/all_edge_score_norm_0.5_0.5.eda" # edge scores

# output filenames
merged_module_file = "results/merged_modules.txt"
overlapped_module_file = "results/overlapped_modules.txt"

### read edge scores
sys.stdout.write("... read graph and edge weight....\n")
G = io.read_net(hn_file, 100)
(score_dics, slabels) = io.read_edge_attrs(scoresfile) ## edge scores from human net
score_dic = score_dics[slabels.index(weight_type)]

# read modules
wdic = {}
(ogM, M_dic) = io.read_module_file(module_file)

sys.stdout.write("merging modules ....\n")
(mgM, cost_dic) = module_cover2.merge_modules(ogM, wth, G, wdic, [score_dic], [1], wth)
io.write_genes_in_modules(mgM, merged_module_file)

sys.stdout.write("overlapping modules ....\n")
ovM = module_cover2.overlap_modules(mgM, 0, G, wdic, [score_dic], [1], wth)
io.write_genes_in_modules(ovM, overlapped_module_file)