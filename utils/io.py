#!/usr/bin/env python

# various I/O functions for processing TCGA data


import sys

import networkx as nx
import numpy as np
import pandas

import utils.misc as misc


def read_net(netfile, top="NA"):
    """ read net file
    Parameters:
		netfile: net file name

			gene1\tgene2\tweight(optional)

		top : top % of weights to be included
			use "NA" if no weights given
	Returns:
		nx.Graph
	"""
    lines = open(netfile).readlines()
    G = nx.Graph()
    tknss = [tuple(l.split()) for l in lines]
    if (top == "NA") or (len(list(filter(lambda tkns: len(tkns) < 3, tknss))) != 0): # unweighted
        tknss = [tuple(tkns[:2]) for tkns in tknss]
        G.add_edges_from(tknss)
    else: # weighted
        weights = [float(x[2]) for x in tknss]
        wth = np.percentile(weights, 100-top)  # weight threshold
        for tkns in tknss:
                if float(tkns[2]) >= wth:
                    G.add_edge(tkns[0], tkns[1], weight=float(tkns[2]))
    return G


def read_net_top(netfile, top=100):
	""" read net file and keep top % edges
	Parameters:
		netfile: string name of net file
		gene1\tgene2\tweight (optional)
		top : top % of weights to be included
			use "NA" if no weights given
	Returns:
		Graph
	"""
	lines = open(netfile).readlines()
	G = nx.Graph()

	## remove self loops
	tknss = [tuple(l.split()) for l in lines]

	if top == "NA" or len(tknss[0]) < 3:
		print("no weights..")
		G.add_edges_from(tknss)
		return G

	weights = [float(x[2]) for x in tknss]
	weights.sort(reverse=True)
	wth = weights[int((len(weights)-1)*top/100)]
	print(wth)
	for tkns in tknss:
		if float(tkns[2]) >= wth:
			G.add_edge(tkns[0], tkns[1], weight=float(tkns[2]))

	return G



def read_dic(filename):
	""" read a file with two columns and create a dic
		the first column as keys and the second as values
	:param filname
	:return: dict first column -> the second column
	"""
	lines = open(filename).readlines()
	data_dic = {}
	for l in lines:
		tkns = l.split()
		if len(tkns) > 2:
			sys.stderr.write("more than two columns")
			return {}
		data_dic[tkns[0]] = tkns[1]
	return data_dic


def read_mut_list(filename, sep=","):
	"""
	read a bipartite graph B(G, S) in the form of list of positive samples
		gene\tcs where csi for comma separated indices (default)
		G: genes, S: samples
	** USE THIS FORMAT FOR UNWEIGHTED CASES. COMPACT FORMAT
	example:
		PTEN	12,14,36,40
	:param filename
	:param sep default=","
	:return dict gene -> list of sample indices
	"""
	temp_dic = read_dic(filename)
	rel_dic_list = {}
	for g in temp_dic:
		rel_dic_list[g] = [int(x) for x in temp_dic[g].split(sep)]

	return rel_dic_list


def read_mut_matrix(filename, mw=3, mutsig_file=None):
	""" read a bipartite graph B(G, S) in the form of a labeled matrix
		G: genes, S: samples
		the first row is the list of sample labels
		the first column is for gene names
		Each element is labeled as N(None), C(CNV), M(Somatic Mutation), or B(Both)
		example:
					s1	s2	s3	s4 	S5  ...
			PTEN	N 	C	N	M	B   ...

		labels are converted to weights based on weight_dic
		weight_dic = dict([('N', 0), ('C', 1), ('M', mw), ('B', mw+1)])
		if mutsig_file is given, multiply mutsig score + 1 to the each scores

	:param filename
	:param mw (the relative weight of somatic mutation) default=3
	:param mutsig_file: gene weight file from mutsig
		ex.
			1 gene mutsig_score
			2 ABCA1 -0.0
			3 ABCA2 0.0862616260817

	:return genes: list of genes
		samples: list of samples
		data_dic: dict gene g -> the list converted from the label to edge weight e(g, s)
	"""
	lines = open(filename).readlines()
	genes = []
	samples = lines[0].split()[1:]
	data_dic = {}
	weight_dic = dict([('N', 0), ('C', 1), ('M', mw), ('B', mw+1)])
	for l in lines[1:]:
		tkns = l.split()
		data_dic[tkns[0]] = [float(weight_dic[x]) for x in tkns[1:]]
		genes.append(tkns[0])

	if mutsig_file is not None: # if mutsig file is given
		mutsig = pandas.read_table(mutsig_file, sep=" ", index_col=0).to_dict()['mutsig_score']
		mutsig_data_dic = {}
		for g in data_dic:
			# multiply (mutsig_score+1) for each gene
			mutsig_data_dic[g] = [(mutsig[g]+1)*x for x in data_dic[g]]
		data_dic = mutsig_data_dic

	return genes, samples, data_dic


def write_mut_matrix(genes, samples, data_dic, filename):
	"""
	write a bipartite graph B(G, S) in the form of a weighted matrix
		G: genes, S: samples
		the first row is the list of sample names
		the first column is for gene names
		example:
				sample1\tsample2\t ...
			PTEN\t0\t1\t\t0\t3..

	:param genes: list of genes
	:param samples: list of samples
	:param data_dic: dict gene -> weights of edges in B(G, S), 0 for no edge
	:param filename
	:return
	"""
	f = open(filename, 'w')
	f.write("gene\t%s\n" % "\t".join(samples))
	for g in genes:
		f.write("%s\t" %g)
		f.write("%s\n" % "\t".join([str(x) for x in data_dic[g]]))
	f.close()
	

def write_mut_list(rel_dic, filename, sep=","):
	"""
	write a bipartite graph B(G, S) in the form of list of positive samples
		gene\tcs where csi for comma separated indices (default)
		G: genes, S: samples
	** USE THIS FORMAT FOR UNWEIGHTED CASES. COMPACT FORMAT
	example:
		PTEN	12,14,36,40
	:param rel_dic: dict gene -> weights of edges in B(G, S), 0 for no edge
	:param filename
	:param sep default=","
	:return
	"""

	list_dic = {}
	for g in rel_dic:
		covers = misc.get_positives(rel_dic[g])
		if len(covers) == 0:
			continue
		list_dic[g] = sep.join([str(x) for x in covers])

	write_dic(list_dic, filename)


def write_dic(any_dic, filename):
	""" write a file with two columns and create a dic
		the first column as keys and the second as values
	:param any_dic
	:param filename
	:return
	"""

	f = open(filename, 'w')
	for x in any_dic:
		f.write("%s\t%s\n" % (str(x), str(any_dic[x])))
	f.close()

def create_type_idx(samples, types, sample_type_dic):
	"""
	create a dict mapping type -> list of sample indices for the type
	the indices given in the samples are used

	:param samples: list of samples
	:param types: list of cancer types
	:param sample_type_dic: sample->type dic
	:return:
	"""
	type_idx_dic = {}
	for ty in types:
		type_idx_dic[ty] = filter(lambda i: sample_type_dic[samples[i]] == ty, range(len(samples)))

	return type_idx_dic




def read_edge_attrs(filename):
	""" read score file between two genes
	file format:
		gene1 gen2 weight
		starting with the labels in the first row
	Returns: tuple of score dics
	"""
	lines = open(filename).readlines()

	if len(lines) == 0:  ## when the file is empty
		return {}, []

	labels = lines[0].split()[1:]
	score_dics = [{} for la in range(len(labels))]
	for l in lines[1:]:
		tkns = l.split("\t")
		(x, etype, y) = tkns[0].split()
		for i in range(len(labels)):
			score_dic = score_dics[i]
			if x not in score_dic:
				score_dic[x] = {}
			if tkns[i + 1] != "NA":
				score_dic[x][y] = float(tkns[i + 1])

	return score_dics, labels


def read_module_file(filename):
	""" read module file generated by greedy_min_cost_module_cover
	:
	:param	filename:
			file format: Each row is the info for each selected gene with the labels in the first row

				itr module_id gene_name min_cover max_cover avg_cover
	:return
		list of modules (a module is a list of genes)
		mod_dic : gene -> (itr, module_id,  min, max, avg_cover)
	"""
	lines = open(filename).readlines()
	modules = []
	mod_dic = {}
	for l in lines[1:]:
		tkns = l.split()
		module_id = int(tkns[1])
		itr = int(tkns[0])
		min = int(tkns[3])
		max = int(tkns[4])
		avg_cover = int(tkns[5])
		if len(modules) <= module_id:  # new module
			modules.append([])
		modules[module_id].append(tkns[2])
		mod_dic[tkns[2]] = (itr, module_id, min, max, avg_cover)

	return modules, mod_dic


def read_simple_module_file(filename, nM_dic=None):
	""" read module file generated by greedy_min_cost_module_cover

	:param filename:
			each row is one module with genes (comma separated)
	:return modules: list of modules (a module is a list of genes)
	:return mod_dic : gene -> (itr, module_id,  min, max, avg_cover, merged_module_id1, 2,...)
		a gene may belong to multiple modules
	"""
	lines = open(filename).readlines()
	modules = [l.strip().split(",") for l in lines]

	if nM_dic is not None:
		M_dic = {}
		for i in range(len(modules)):
			m = modules[i]
			for g in m:
				if g not in M_dic:
					M_dic[g] = nM_dic[g]
				M_dic[g] = M_dic[g] + (i,)

	else:
		M_dic = None

	return modules, M_dic


def write_genes_in_modules(M, filename):
	""" write genes in each module in one line (comma separated)
	M: modules
	filename
	"""
	f = open(filename, 'w')
	for m in M:
        	f.write("%s\n" %",".join(m))

	f.close()



