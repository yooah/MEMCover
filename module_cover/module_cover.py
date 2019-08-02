#!/usr/bin/env python
# module cover library

import logging
import time
import heapq as hq

import utils.misc as misc

def comp_cost(M, gene, Mavg, sum_dic, th, best_mod=-1):
	"""
	compute the cost of adding a gene g to the current module set
	find module which incurs the minimum cost
	:param M:  a list of currently selected modules
	:param gene: gene to be considered
	:param Mavg: Mavg[i] - the current avg cost of module i
	:param sum_dic: sum_dic[z][i] - sum of cost from z to module i
	:param th:
	:param best_mod:
	:return:min_cost: the min cost of adding g
	:return min_cost_m: module index to which gene g is to be added
			if it is len(M), g is added as a separate module
	"""
	# cost of adding the node as a separate node
	MAX = 100
	min_cost_m = len(M)
	min_cost = 1
	min_avg = MAX
	if best_mod == -1:
		return min_cost, min_cost_m
	else:
		module = M[best_mod]
		avg_corr = sum_dic[gene][best_mod] / len(module)
		add_cost = 1 + (Mavg[best_mod] / len(module) - 2 * avg_corr + th)  # based on new weight definition
		if add_cost < min_cost or (add_cost == min_cost and min_avg < avg_corr):
			min_cost = add_cost
			min_cost_m = best_mod

	return min_cost, min_cost_m


def update_modules(max_id, modules, max_g, Mavg, max_cost, th):
	"""
	Update modules by adding max_g
	and Mavg

	:param max_id: module id to be added
	:param modules: current module set
	:param max_g: gene to be added
	for Mavg
	:param max_cost: cost difference
	:param Mavg: original Mavg
	:param th: module cost threshold

	:return: modules: updated modules
			 mavg: updated module cost for max_id
	"""
	if max_id < len(modules):
		modules[max_id].append(max_g)
		mavg = Mavg[max_id] + 1 - max_cost + th
	else:
		modules.append([max_g])
		mavg = th
	return (modules, mavg)


def update_sum_dic(module_nodes, neighbors, sum_dic, max_module, score_dic, max_g):
	"""
	update sum_dic
	sum_dic[x][i]: sum of weights between x and a node in module_nodes

	:param module_nodes: nodes in the new module
	:param neighbors: all neighbors of module_nodes
	:param sum_dic: old sum_dic
	:param max_module: new module index
	:param score_dic: edge score dic
	:param max_g: selected gene
	:return:

	"""
	for neigh in neighbors:
		# update sum_dic
		# if the neighbor already has an entry
		if neigh in sum_dic and max_module in sum_dic[neigh]:
			score = misc.get_val(score_dic, neigh, max_g)
			sum_dic[neigh][max_module] += score
		# if the neighbor is new (i.e., the neighbor of max_g)
		else:
			cost_list = []
			for x in module_nodes:
				score = misc.get_val(score_dic, neigh, x)
				cost_list.append(score)
			if neigh not in sum_dic:
				sum_dic[neigh] = {}
			sum_dic[neigh][max_module] = sum(cost_list)
	if max_g in sum_dic: # gene removed once selected
		del sum_dic[max_g]


def update_best_heap(best_heap, best_dic, neighbors, max_module, M,  sum_dic, Mavg, th):
	"""
	update best_dic and best_heap (using priority queue)
	with a module updated, best module for each gene is updated

	:param best_heap: a heap representation of best_dic
	:param best_dic: node -> module -> cost
	:param neighbors: all nodes to be updated
	:param max_module: module to be updated
	:param M: modules
	:param sum_dic: dict g -> module -> sum of edge scores between g amd module
	:param Mavg: dict average score within module
	:param th: edge score threshold
	:return:
	"""
	for neigh in neighbors:
		# the average score from neigh to module should be > th
		if sum_dic[neigh][max_module] < th * len(M[max_module]):
			if max_module in best_dic[neigh]:
				del best_dic[neigh][max_module]
				best_heap[neigh] = [(best_dic[neigh][x], x) for x in best_dic[neigh]]
				hq.heapify(best_heap[neigh])
			continue
		new = (Mavg[max_module] - 2 * sum_dic[neigh][max_module]) / len(M[max_module]) + th
		best_dic[neigh][max_module] = new
		if len(best_heap[neigh]) > 0 and best_heap[neigh][0][1] == max_module:
			hq.heapreplace(best_heap[neigh], (new, max_module))
		else:
			best_heap[neigh] = [(best_dic[neigh][x], x) for x in best_dic[neigh]]
			hq.heapify(best_heap[neigh])


def update_cover_info(selected_g, revised_mut_dic, sample_cover_list, uncovered_samples):
	"""
	update coverage information

	:param selected_g: selected gene
	:param revised_mut_dic: copy of mutation dic
			for covered sample, set value to be 0 so that it is not counted in the cover
			remove selected_g from the dict so that it is not selected again
			revised each iteration in module cover
	:param sample_cover_dic: number of times TO BE COVERED for each sample index
	:param uncovered_samples: set of samples not completely covered yet

	:return: new_covered: list of newly covered sample indices
	"""
	# decrease cover count for each sample
	for i in uncovered_samples:
		sample_cover_list[i] -= revised_mut_dic[selected_g][i]

	# a gene cannot be selected multiple times
	revised_mut_dic.pop(selected_g)

	new_covered = set(filter(lambda x: sample_cover_list[x] <= 0, uncovered_samples))
	for c in new_covered:
		for g in revised_mut_dic:
			revised_mut_dic[g][c] = 0

	return new_covered


# main functions
def greedy_module_cover(mut_dic, GNet, k, score_dic, th=0, output=None, stop=False, l=0):
	"""
	Find a module cover
	module weight = alpha + # nodes - sum of average weight for all nodes
	:param mut_dic: dict[g] = list of weights (w[g][i] > 0 iff g covering sample i)
	:param GNet: interaction network
	:param k: number of times a sample is covered
	:param score_dic: edge score dic score_dic[x][y]
	:param th: module cost threshold
	Optional
	:param output: file to write the results
	:param stop: boolean : stop module cover if the best module is singleton covering only one sample
	Outdated
	:param: l: number of outliers
	:return:
		M: list of selected modules
		total_cost
	"""

	# initializing...
	# samples and genes
	nsamples = len(mut_dic.values()[0])
	nodes = set(mut_dic).intersection(GNet)
	uncovered = set(range(nsamples))  # uncovered samples
	sample_cover_count = [k for x in uncovered]  # number of times covered for each sample

	revised_dic = dict([(g, list(mut_dic[g])) for g in mut_dic])  # deep copy of mut_dic
	# when partial cover is given..

	# selected modules, selected genes, total module cost
	M, selected, total_cost = [], [], 0

	# precomputed dics to make the module cost update simple
	# Mavg[i] - the avg cost of module i
	# sum_dic: sum_dic[z][i] - sum of edge scores from z to module i
	Mavg, sum_dic = {}, {}

	# For running time optimization
	# best_dic: dictionary for the sum cost for (node, module) pair
	# best_heap: a heap representation of best_dic
	# makes it easier to find the best module and update the cost
	best_dic, best_heap = {}, {}
	for g in GNet:
		best_dic[g] = {}
		best_heap[g] = []

	# for debugging
	prevtime = time.mktime(time.localtime())
	starttime = prevtime

	if output is not None:
		f = open(output, 'w')
		f.write("itr\tmodule_id\tselected_gene\tmost_covered\tleast_covers\tavg_covered\tbenefit\tcost\n")

	itr = len(selected)
	max_g = None

	while len(uncovered) > l:
		itr += 1
		curtime = time.mktime(time.localtime())
		logging.debug("%f sec taken " % (curtime - prevtime))
		logging.debug("(%f sec taken in total)\n" % (curtime - starttime))
		prevtime = curtime

		# find the best node to add
		max_ben_cost = 0
		for g in nodes:
			# when best module for g is known
			if len(best_heap[g]) > 0:
				bestg = best_heap[g][0][1]
			else:
				bestg = -1
			(cost, m) = comp_cost(M, g, Mavg, sum_dic, th, bestg)
			new_module = m
			benefit = sum(revised_dic[g])
			ben_cost = benefit / float(cost)
			if ben_cost < max_ben_cost:
				continue
			if ben_cost > max_ben_cost or (max_g is not None and sum(mut_dic[g]) > sum(mut_dic[max_g])):
				# if the gene is better than the current i
				# or ties with the best and its original benefit is better
				(max_g, max_ben_cost, max_ben, max_cost, max_module) = (g, ben_cost, benefit, cost, new_module)

		# STOP condition
		# bug fixed for the last node being added with benefit = 0
		if (max_g in selected) or (max_ben_cost == 0) or (stop is True and max_ben == 1 and max_cost == 1):
			logging.info("stop condition met..")
			break

		# add max_g to module set (M), update module cost, gene sets, total module cost, sample set
		(M, Mavg[max_module]) = update_modules(max_module, M, max_g, Mavg, max_cost, th)
		selected.append(max_g)  # max_g is selected
		nodes.remove(max_g)  # max_g is removed from available genes
		new_covered = update_cover_info(max_g, revised_dic, sample_cover_count, uncovered)  # revised_dic, cover_count
		uncovered = uncovered.difference(new_covered)  # remove covered samples
		total_cost += max_cost  # total cost

		# update sum_dic, best_dic, best_heap
		module_nodes = M[max_module]
		neighbors = misc.neighbors(GNet, module_nodes).difference(selected)
		update_sum_dic(module_nodes, neighbors, sum_dic, max_module, score_dic, max_g)
		update_best_heap(best_heap, best_dic, neighbors, max_module, M, sum_dic, Mavg, th)
		# write progress
		print_debug_message(itr, max_module, max_ben, max_cost, max_ben_cost, mut_dic[max_g], max_g)

		# write in the output file
		if output is not None:
			f.write("%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\n" % (itr, max_module, max_g,
					k-max(sample_cover_count), k-min(sample_cover_count),
					k-sum(sample_cover_count)/nsamples, max_ben, max_cost))

	if output is not None:
		f.close()

	return M, total_cost


# print function
def print_debug_message(itr, max_module, max_ben, max_cost, max_ben_cost, rel_g, max_g, interval=1):
	if itr % interval == 0:
		print str(itr) + ": " + max_g
		logging.debug("%d-th iteration----------------\n" % itr)
		logging.debug("module id that the selected node will be added to :  %d\n" % max_module)
		logging.debug("benefit is %d , cost is %d, max ben/cost is %d\n " % (max_ben, max_cost, max_ben_cost))
		logging.debug("original benefit %d\n" % sum(rel_g))
		logging.debug("best node is %s\n" % max_g)

