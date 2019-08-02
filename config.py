#!/usr/bin/env python
# define global variables

import os
# 11 cancer types used in the pan-cancer analysis
cancers = ["BLCA", "BRCA", "CRC", "GBM", "HNSC", "KIRC", "LAML", "LUAD", "LUSC", "OV", "UCEC"]

# DIRECTORIES
cur_dir = os.getcwd()
results_dir = cur_dir+"/results/"

data_dir = "./data/"
permute_dir = results_dir
me_dir = results_dir


hn_file = data_dir+"human_net.net" # HumanNet (gene interaction)
mut_file = data_dir+"mut_data.txt" # mutation file

# subtype information
# format: patient label in the first column followed by cancer type
# the labels used in the mutation file and cancer types from config.cancers
subtype_file = data_dir+"subtype.txt"  # cancer type for each sample

# constants
nsamples = 3182  # total number of samples
