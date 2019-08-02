#!/usr/bin/env python
# module cover library 2

"""
module cover routines
optimized version using heap
"""

import sys
import csv
import os
import numpy as np
import networkx as nx
import utils.misc as misc
import time
import heapq as hq

### utility functions
def update_w_dic(wdic, set_a, set_b, score_dics, coefs):
	""" 
	update weight dic between set_a and set b
	w(x, y) = [score,mi]*coefs-th (see comp_w) for all pairts x in set_a, y in set_b
	score = score_dic[x][y] (0 if there is no edge)

	Parameters:
		wdic: weight dic (gene->(gene->weight))
		set_a, set_b: compute weights between the two
		score_dics:  list of interaction scores
		g: file handler
	Return wdic
	"""
	for x in set_a:
		if x not in wdic:
			wdic[x] = {}
		for y in set_b:
			if y not in wdic[x]:
				scores = []
				for sd in score_dics:
					if x in sd and y in sd[x]:
						score = sd[x][y]
					else:
						score = 0
					scores.append(score)
				wdic[x][y] = comp_w(scores, coefs)
			if y not in wdic:
				wdic[y] = {}
			wdic[y][x] = wdic[x][y]

	return wdic

def comp_mcost(aset, wdic, th):
	"""
	compute module cost
	definition:
	m > 1: m (1+th) - sum_x sum_ (y not x) wdic[x][y]/ (m-1)
	m == 1: 1
	"""
	m = len(aset)
	if m == 1:
		return 1
	## total coherence
	tcoh = sum([sum([wdic[g][h] for h in filter(lambda h: h != g, aset)]) for g in aset])

	return (m * (1 + th) - tcoh / (m - 1))

def comp_w(values, coefs, avg=0):
	""" 
	compute weight of an edge in module cover
	weighted sum of values minus avg
	Parameters:
		values: list of floats
		coefs: list of floats
		avg: float (not to be used, left just in case 09/29)
	Return:
		float	
	"""
	return sum([values[i] * coefs[i] for i in range(len(values))]) - avg

def find_best_gene_module_pair(modules, G, gm_dic, m_dic):
	""" compute the between cost of every pair of modules
	and pick a pair with the maximum cost > 0

	Parameters:
		modules
		G
		gm_dic: module -> gene -> mcost(m+g)
		m_dic: module cost
	Returns:
		g : max pair gene
		m : max pair module
		max_cost
	"""
	max_cost = -1
	max_pair = (None, None)
	sys.stderr.write("finding best (gene, module) pairs...\n")
	all_genes = []
	for i in range(len(modules)):
		m = modules[i]
		for g in gm_dic[i]:
			if g in m:
				continue
			cost = 1 + m_dic[i] - gm_dic[i][g]
			if max_cost < cost:
				max_cost = cost
				max_pair = (g, m)
	return (max_pair[0], max_pair[1], max_cost)

def find_best_module_pair(modules, mcost_dic):
	""" compute module cost for all moduels and
	all possible merged modules

	Parameters:
		modules
		mcost_dic mcost_dic[i1][i2] = mcost(m1 + m2)
		if i1=i2, most(mi)
	Returns:
		module1
		module2
		score (mcost(m1) + mcost(m2)- mcost(m1+m2)
	"""
	max_cost = -1
	max_pair = (None, None)
	sys.stderr.write("finding best module pairs...\n")
	for i in mcost_dic:
		m1 = modules[i]
		for j in mcost_dic[i]:
			if i == j:
				continue
			m2 = modules[j]
			cost = mcost_dic[i][i] + mcost_dic[j][j] - mcost_dic[i][j]
			if max_cost < cost:
				max_cost = cost
				max_pair = (m1, m2)
	return (max_pair[0], max_pair[1], max_cost)

## compute module cost (all modules and all possible merged modules)
def comp_mcost_all(modules, G, wdic, score_dics, coefs, th):
	"""
	mcost_dic[i][j] = mcost(mi+mj)
	Parameters:
		modules
		G: Graph (necessary for optimizations - costs are computed only if two modules have neighbors)
		wdic: weight dic
		score_dics:  list of interaction scores
		coef
		th
		g
	Returns:
		cost_dic: between module cost
	"""
	mcost_dic = {}
	for m in modules:
		wdic = update_w_dic(wdic, m, m, score_dics, coefs)
	for i in range(len(modules)):
		mcost_dic[i] = {}
		m1 = modules[i]
		mcost_dic[i][i] = comp_mcost(m1, wdic, th)
		for j in range(len(modules)):
			m2 = modules[j]
			if i == j or len(set(misc.neighbors(G, m1)).intersection(m2)) == 0:
				continue
			wdic = update_w_dic(wdic, m1, m2, score_dics, coefs)
			mcost_dic[i][j] = comp_mcost(set(m1).union(m2), wdic, th)
	return mcost_dic

def comp_gmcost_all(modules, G, wdic, th):
	"""
	gmcost_dic[i][g] = mcost(mi+g)
	mcost_dic[i] = mcost(mi)
	Parameters:
		modules
		G
		wdic (computed based on score_dics, coefs)
	Returns:
		gmcost_dic
		mcost_dic
	"""
	gmcost_dic = {}
	mcost_dic = {}
	all_genes = []
	for m in modules:
		all_genes.extend(m)
	for i in range(len(modules)):
		m1 = modules[i]
		mcost_dic[i] = comp_mcost(m1, wdic, th)

		gmcost_dic[i] = {}
		neighbors = set(all_genes).intersection(misc.neighbors(G, m1))
		if i == 57:
			print(neighbors)
		for g in neighbors:
			gmcost_dic[i][g] = comp_mcost(set(m1).union([g]), wdic, th)
	return (mcost_dic, gmcost_dic)

def comp_between_cost_all(modules, G, wdic, score_dics, coefs, g=None):
	""" compute the between cost of every pair of modules
	cost_dic[i][j] = sum(wdic[x][y] for all x,y)/len(modules(i))*len(modules(j)) (see comp_between_cost)
	i, j: module indices, x in modules[i], y in modules[j]
	Parameters:
		modules
		G: Graph (necessary for optimizations - costs are computed only if two modules have neighbors)
		wdic: weight dic
		score_dics:  list of interaction scores
		coef
		th
		g
	Returns:
		cost_dic: between module cost
	"""
	cost_dic = {}
	for i in range(len(modules)):
		if i % 100 == 0:
			print(i)
		m1 = modules[i]
		neighs = misc.neighbors(G, m1)
		cost_dic[i] = {}
		neigh_modules = filter(lambda j: len(set(neighs).intersection(modules[j])) > 0, range(i + 1, len(modules)))
		for j in neigh_modules:
			m2 = modules[j]
			wdic = update_w_dic(wdic, m1, m2, score_dics, coefs)
			cost_dic[i][j] = comp_between_cost(m1, m2, wdic)
	return cost_dic

def comp_between_cost(m1, m2, wdic):
	""" compute the cost between two modules
	"""
	total_cost = 0
	total_cost = sum([sum([wdic[x][y] for y in m2]) for x in m1])
	return total_cost / (len(m1) * len(m2))

### functions for post processing
def merge_modules(modules, alpha2, G, wdic, score_dics, coefs, th):
	""" merge modules after obtaining modules from module_cover
	merge two modules if the merged module has a lower cost
	mcost(m1) + mcost(m2) >= mcost(m1+m2) + alpha2
	Paramters:
		modules
		alpha2: theshold
		G graph (for optimization)
		wdic: weight dic gene -> gene -> weight (see update_wdic)
			computed using score_dics, coefs, th
	Returns:
		new_modules: merged modules
		cost_dic: cost between modules (after merging)
	"""
	## compute module cost (all modules and all possible merged modules)
	old_modules = modules
	modules = [list(m) for m in old_modules]  ## make a copy
	mcost_dic = comp_mcost_all(modules, G, wdic, score_dics, coefs, th)
	while True:
		## score =  (mcost(m1) + mcost(m2)) - mcost(m1+m2)
		(m1, m2, score) = find_best_module_pair(modules, mcost_dic)
		print(score)
		if score <= -alpha2:
			break
		else:
			i1 = modules.index(m1)
			i2 = modules.index(m2)
			sys.stderr.write("(%s) and (%s) are merged (cost %f)\n" % (",".join(m1), ",".join(m2), score))
			modules[i1].extend(m2)
			modules[i2] = []
			mcost_dic[i1][i1] = mcost_dic[i1][i2]  ## new mcost(mi)
			del mcost_dic[i1][i2]  ### mcost_dic[i][j] = mcost(mi+mj)
			del mcost_dic[i2][i2]  ## mcost_dic[i][i] = mcost(mi)
			for i in set(mcost_dic[i1]).union(mcost_dic[i2]):
				if i == i1 or i == i2:
					continue
				wdic = update_w_dic(wdic, modules[i], modules[i1], score_dics, coefs)
				mcost_dic[i][i1] = comp_mcost(set(modules[i]).union(modules[i1]), wdic, th)
				mcost_dic[i1][i] = mcost_dic[i][i1]
				if i2 in mcost_dic[i]:
					del mcost_dic[i][i2]
			del mcost_dic[i2]
	new_modules = list(filter(lambda x: len(x) > 0, modules))
	return (new_modules, mcost_dic)


def overlap_modules(modules, alpha3, G, wdic, score_dics, coefs, th):
	""" compute overlapping modules after obtaining modules from module_cover + merging
	choose a (gene, module) pair with maximum cost (>0) and add the gene to the module
	if gene is not in the module
	Paramters:
		modules
		alpha3: threshold
		G graph (for optimization)
		wdic: weight dic gene -> gene -> weight
			(assume every pair (to be considered) has an entry in wdic,
			which may not be true)
			computed using score_dics, coefs, th
	Returns:
		new_modules: overlapped modules
	"""
	new_modules = [list(m) for m in modules]
	(mcost_dic, gmcost_dic) = comp_gmcost_all(modules, G, wdic, th)

	while True:
		(g, m, score) = find_best_gene_module_pair(new_modules, G, gmcost_dic, mcost_dic)
		print(score)
		if score <= -alpha3:
			break
		else:
			mid = new_modules.index(m)
			sys.stderr.write("(%s) is added to (%s) (cost %f)\n" % (g, ",".join(m), score))
			wdic = update_w_dic(wdic, [g], misc.neighbors(G, m), score_dics, coefs)
			wdic = update_w_dic(wdic, m, misc.neighbors(G, [g]), score_dics, coefs)
			new_modules[mid].append(g)
			del gmcost_dic[mid][g]
			for g2 in gmcost_dic[mid]:
				wdic = update_w_dic(wdic, [g2], new_modules[mid], score_dics, coefs)
				gmcost_dic[mid][g2] = comp_mcost(set(new_modules[mid]).union([g2]), wdic, th)
	return new_modules




