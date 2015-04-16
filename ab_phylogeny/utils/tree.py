#!/usr/bin/python
# filename: tree.py



###########################################################################
#
# Copyright (c) 2014 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################


import os
import subprocess as sp

import ete2



def make_tree(alignment, aa):
	'''
	Builds a tree file (using FastTree) from a sequence alignment in FASTA format

	Input
	path to a FASTA-formatted sequence alignment

	Output
	path to a Newick-formatted tree file
	'''
	tree = alignment.replace('_aligned.aln', '_tree.nw')
	return fast_tree(alignment, tree, aa)



def fast_tree(alignment, tree, aa):
	if aa:
		ft_cmd = 'fasttree -nt {} > {}'.format(alignment, tree)
	else:
		ft_cmd = 'fasttree {} > {}'.format(alignment, tree)
	ft = sp.Popen(ft_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
	stdout, stderr = ft.communicate()
	return tree


def make_figure(tree, timepoints):
	fig = tree.replace('_tree.nw', '_tree.pdf')
	orders = {tp.name: tp.order for tp in timepoints}
	colors = {tp.name: tp.color for tp in timepoints}
	# make the tree
	t = ete2.Tree(tree)
	# style the nodes based on timepoint
	for node in t.traverse():
		earliest = get_earliest_leaf(node.get_leaf_names(), orders)
		color = colors[earliest]
		style = ete2.NodeStyle()
		style['size'] = 0
		style['vt_line_width'] = 1.0
		style['hz_line_width'] = 1.0
		style['vt_line_color'] = color
		style['hz_line_color'] = color
		style['vt_line_type'] = 0
		style['hz_line_type'] = 0
		node.set_style(style)
	# style the full tree
	ts = ete2.TreeStyle()
	ts.show_leaf_names = False
	ts.scale = 500
	ts.branch_vertical_margin = 1.5
	# render the tree
	t.render(fig, tree_style=ts)


def get_earliest_leaf(leaves, order):
	counts = {}
	for leaf in leaves:
		tp = leaf.split('_')[0]
		counts[tp] = counts[tp] + 1 if tp in counts else 1
	total = sum(counts.values())
	if 'root' in counts:
		return 'root'
	timepoints = sorted(counts.keys(), key=lambda x: order[x])
	for tp in timepoints:
		if 100. * counts[tp] / total >= 10:
			return tp
