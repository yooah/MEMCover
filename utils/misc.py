#!/usr/bin/env python
# miscellaneous utility functions


def get_val(dic, g1, g2):
	"""
	:param dic: any 2-d dic
	:param g1: key 1
	:param g2: key 2
	:return: return dic[g1][g2] if exists, 0 otherwise
	"""
	if g1 in dic and g2 in dic[g1]:
		s = dic[g1][g2]
	elif g2 in dic and g1 in dic[g2]:
		s = dic[g2][g1]
	else:
		s = 0
	return s


def get_positives(in_list):
	""" given a list of values, return positive indices
	"""
	return filter(lambda i: in_list[i] > 0, range(len(in_list)))

def neighbors(G, aset):
	"""
	find neighbor nodes given a set of nodes

	:param G: nx.Graph
	:param aset: a set of nodes
	:return: neighbors: list
	"""
	neighs = set([])
	for x in aset:
		neighs.update(G.neighbors(x))
	return neighs




