#!/usr/bin/python
# compute mutual exclusivity rank with permutation test

import logging
import numpy as np


def comp_type_dic(original_gene_dic, type_idx_dic):
	"""
	divide the given mut_dic into each type
	:param original_gene_dic: gene-> set of covered samples
	:param type_idx_dic: cancer type -> set of samples
	:return:original_gene_type_dic: cancer type -> gene -> set of covered samples
	"""

	original_gene_type_dic = {}
	for can in type_idx_dic:
		temp_dic = {}
		can_idx = type_idx_dic[can]
		for gene in original_gene_dic:
			temp_dic[gene] = set(original_gene_dic[gene]).intersection(can_idx)
		original_gene_type_dic[can] = temp_dic

	return original_gene_type_dic


def comp_pair_cover(gene_mut_dic, pairs):
	"""
	given a mut_dic, compute the size of cover for each pair
	if either gene does not exist, cover size is 0
	:param gene_mut_dic: gene -> list of covered samples
	:param pairs: list of gene pairs
	:return:
		cover_size_list: cover sizes in the same order as pairs
	"""
	cover_size_list = []

	for pair in pairs:
		(x, y) = pair
		if x not in gene_mut_dic or y not in gene_mut_dic:
			cover = []
		else:
			cover = set(gene_mut_dic[x]).union(gene_mut_dic[y])
		cover_size_list.append(len(cover))

	return cover_size_list


def norm_cover_size(type_cover_list, coefs):
	"""
	compute normalized sum of each element
	:param type_cover_list: can -> list of cover_sizes
	:param coefs: can -> coef
	:return:
		[sum([type_cover_list[can][i]/coefs[can] for can in cancers]) for i in range(len(type_cover_list) ]
	"""
	temp_type_cover = []
	for can in coefs:
		coef = coefs[can]
		type_cover = type_cover_list[can]
		temp_type_cover.append([x*coef for x in type_cover])
	return [sum(elems) for elems in zip(*temp_type_cover)]


def update_rank(rank_list, original_cover_size_list, permuted_cover_size_list, ep=1e-5):
	"""
	update ranking based on new permuted instance
	compute a vector to indicate 1 if original is no bigger than permuted
	and add to the rank
	:param rank_list: list of rank
	:param original_cover_size_list: list of original cover sizes
	:param permuted_cover_size_list: list of permuted cover sizes
	:param ep: epsilon for numerical precision error (conservatively compute the rank)
	:return: rank_list: updated rank list (new)
	"""

	llen = len(original_cover_size_list) # length of the list
	if llen != len(permuted_cover_size_list):
		logging.warning("two lists have different sizes\n")
	temp_ranks = [int(original_cover_size_list[i] <= permuted_cover_size_list[i] + ep) for i in range(llen)]
	rank_list = [x+y for x, y in zip(temp_ranks, rank_list)]

	return rank_list


def comp_pv(all_ranks, ptype, pnum, cancers):
	"""
	compute p-values

	:param all_ranks: pandas data frame
	:param ptype: permutaion type (tr or to)
	:param pnum: number of permutation instances
	:param cancers: cancer type

	:return: pvs: pvalues pandas data frame
	"""
	pvs = all_ranks.copy()
	pv_labels = ["gene1", "gene2", "raw_pv"]
	if ptype == "tr":
		pv_labels += (["norm_pv"]+[can+"_pv" for can in cancers])
	pvs.columns = pv_labels  # overwrite labels
	pvs.iloc[:, 2:] += 1
	pvs.iloc[:, 2:] /= float(pnum+1)

	return pvs


def comp_logp(pvs, ptype, cancers):
	"""
	compute log2 pvalues

	:param pvs: pvalues, pandas data frame
	:param ptype: permutaion type (tr or to)
	:param pnum: number of permutation instances
	:param cancers: cancer type

	:return: logps: pvalues pandas data frame
	"""
	logps = pvs.copy()
	logp_labels = ["gene1", "gene2", "raw_logp"]
	if ptype == "tr":
		logp_labels += (["norm_log"]+[can+"_logp" for can in cancers])
	logps.columns = logp_labels  # overwrite labels
	logps.iloc[:, 2:] = -np.log10(logps.iloc[:, 2:])

	return logps
