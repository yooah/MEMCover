--------------------------------------
Introduction
--------------------------------------

This README file includes the description of pan_me package,
which take a mutation profiles from multipe types of cancers as an input
and compute pairwise mutual exclusivity patterns of genes
 
For the full description of the method, see

[1] YA Kim, et al. MEMCover: Integrated Analysis of Mutual Exclusivity and
Functional Network Reveals Dysregulated Pathways Across Multiple Cancer
Types, ISMB, July, 2015 (Bioinformatics 2015)

Last revised on 08/02/19

********* Contents ***********
- How to use
- Directories/Files List

--------------------------------------
How to use
--------------------------------------
0. To check/modify parameters, see config.py

1. create permutation instances

# the following command reads mut_data.txt in data directory,
# creates a permuted mutation profile
# by swapping edges 100*|E| times
# default method is "Type Resticted/Separated" permutation
# the results are stored in "tr_permuted_cover_0.txt" in results directory
# change OUTPUT_FILES for a different filename

>> python run_permute_cover.py 0 100 tr -mut=mut_data.txt

# use "to" option  to have "Type olivious" permutation

>> python run_permute_cover.py 0 100 to -mut=mut_data.txt

# See below for mut_data.txt format

# if there is only one cancer type, tr and to will give the same results
# be sure to change config.py to define cancer types, number of samples etc.

2. compute empirical p-values via permutation test

#  to read 10 TR permutation files starting from permutation file index 0
#  compute ranks for all pairs given in efile

>> python run_ptest.py 0 10 tr -ef human_net.net

3. Run Module Cover

# run module cover with edge weight threshold = 0.2, k=15 (coverage)
# using "mut_data.txt" for mutation profile
# "human_net.net" for network (include edge weight)
# "combined" for output prefix

>> python run_module_cover.py 15 0.2 -mut=mut_data.txt -net=human_net.net -out=combined

# * positional arguments:
#   k                     number of cover for each patient
#   ewth                  edge weight threshold, 
#			  to choose the threshold, check th edge weight distribution of your data 
#	                  and try different values in the range of edge weights
#			  in general, smaller threshold will create bigger modules
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
# 			  for example, if mutw is 2, the wegith of somatic mutation (M) is 2 whereas CNV (C) weight is 1
#			  if both (B), the weight is 3 (1 for CNV and 2 for mutation)
# 			  see the format of mut_data.txt below.
#			   
#   -nw NW_FILE, --nw_file NW_FILE
#                         node weight file name
#   -out OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
#                        specify output filename prefix
#
#
# * results will be stored in
#   "results/combined_itr_15_0.2.txt": each gene in a row in the order added
#   "results/combined_modules_15_0.2.txt": one module in a row


>> python postproc_module_cover.py

# read modules from "results/norm_modules.txt"
# perform postprocessing (merging + overlapping)
# edge scores from "data/all_edge_score_norm_all_type_0.5_0.5.eda" is used
# with weight type = "norm_combined" and wth = 0.18540455679 (corresponding to the top 40%)

# run postprocessing
--------------------------------------
Directories/Files
--------------------------------------

./
    README.txt
    config.py
    run_permute_cover.py
    run_ptest.py
    run_module_cover.py

    data/
        human_net.net
        # interaction (edge) data file
        # format: Each row is an interaction
        # two genes in the first two columns
        # third column with weight

        mut_data.mat
        #   gene/sample mutation file
        #   format: the first row is the list of sample labels
        #   the first column is for gene names
        #   Each element is labeled as N(None), C(CNV), M(Somatic Mutation), or B(Both)
        #   Example:
        #				s1	s2	s3	s4 	S5  ...
        #		PTEN	N 	C	N	M	B   ...
        # 	    RB1		B	N	C	N	N

        subtype.txt

        all_edge_score_norm_0.5_0.5.eda
        # edge scores file for all edges


    mut_ex/
        mut_ex.py
        permute_mut_data.py

    utils/
        io.py
        misc.py

    results/
        to_permuted_cover_0.txt  # to permutation instance
        tr_permuted_cover_0.txt  # tr permutation instance
        norm_modules.txt         # modules found by module cover
        tr_human_net_me_rank_0_10000.txt # ME ranks for all humannet edges using 10,000 TR permutations
		merged_modules.txt       # modules after merging
		overlapped_modules.txt   # modules after overlapping
		hint_all_me_pvs.txt.gz   # HINT network edge scores
