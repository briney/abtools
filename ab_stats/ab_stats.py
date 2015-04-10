#!/usr/bin/env python
# filename: ab_stats.py


###########################################################################
#
# Copyright (c) 2015 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################


from __future__ import print_function

import os
import sys
import uuid
import time
import shelve
import urllib
import sqlite3
import tempfile
import argparse
import subprocess as sp
import multiprocessing as mp
from StringIO import StringIO
from collections import Counter

import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import seaborn as sns

from pymongo import MongoClient



parser = argparse.ArgumentParser("Computes and plots basic repertoire information from one or more antibody NGS datasets.")
parser.add_argument('-o', '--output', dest='output', required=True,
					help="Output directory for the FASTA files. Required.")
parser.add_argument('-t', '--temp', dest='temp_dir', required=True,
					help="Directory for temporary files. Required.")
parser.add_argument('-d', '--database', dest='db', required=True,
					help="Name of the MongoDB database to query. Required.")
parser.add_argument('-c', '--collection', dest='collection', default=None,
					help="Name of the MongoDB collection to query. \
					If not provided, all collections in the given database will be processed iteratively.")
parser.add_argument('-i', '--ip', dest='ip', default='localhost',
					help="IP address for the MongoDB server.  Defaults to 'localhost'.")
parser.add_argument('-u', '--user', dest='user', default=None,
					help="Username for the MongoDB server. Not used if not provided.")
parser.add_argument('-p', '--password', dest='password', default=None,
					help="Password for the MongoDB server. Not used if not provided.")
parser.add_argument('-V', '--var', dest='var_plot', default=None,
					help="Plot the distribution of variable families or genes. \
					Options are 'fam', 'gene', or 'both'. \
					If not provided, no plot will be made.")
parser.add_argument('-D', '--div', dest='div_plot', default=None,
					help="Plot the distribution of diversity families or genes. \
					Options are 'fam', 'gene', or 'both'. \
					If not provided, no plot will be made.")
parser.add_argument('-J', '--join', dest='join_plot', default=None,
					help="Plot the distribution of joining families or genes. \
					Only option is 'gene'. \
					If not provided, no plot will be made.")
parser.add_argument('-3', '--cdr3', dest='cdr3_plot', default=None,
					help="Plot the distribution of CDR3 lengths, as nucleotide or amino acid. \
					Options are 'nt', 'aa', or 'both'. \
					If not provided, no plot will be made.")
parser.add_argument('-H', '--heatmap', dest='vj_heatmap', default=False, action='store_true',
					help="If set, generates a heatmap (or quiltplot, if you're a novice) of combined variable and joining gene use. \
					If not provided, no plot will be made.")
parser.add_argument('-I', '--isotype', dest='isotype_plot', default=False, action='store_true',
					help="Plot the isotype frequency. \
					Requires that isotypes have been identified using AbStar, or the 'isotype' field is present in the given MongoDB collection. \
					If not provided, no plot will be made.")
parser.add_argument('-C', '--chain', dest='chain', default='heavy', choices=['heavy', 'kappa', 'lambda'],
					help="Select the antibody chain to analyze. \
					Options are 'heavy', 'kappa', and 'lambda'. \
					Default is 'heavy'.")
parser.add_argument('--debug', dest='debug', action='store_true', default=False,
					help="If set, will run in debug mode.")
args = parser.parse_args()




def get_db():
	if args.user and args.password:
		password = urllib.quote_plus(password)
		uri = 'mongodb://{}:{}@{}'.format(args.user, password, args.ip)
		conn = MongoClient(uri)
	else:
		conn = MongoClient(args.ip, 27017)
	return conn[args.db]


def get_collections(db):
	if args.collection:
		return [args.collection, ]
	collections = db.collection_names(include_system_collections=False)
	return sorted(collections)



def query(db, collection):
	print('Getting sequences from MongoDB...', end='')
	sys.stdout.flush()
	c = db[collection]
	results = c.find({'chain': args.chain, 'prod': 'yes'},
					 {'_id': 0,
					  'v_gene.gene': 1, 'd_gene.gene': 1, 'j_gene.gene': 1,
					  'v_gene.fam': 1, 'd_gene.fam': 1,
					  'cdr3_len': 1, 'isotype': 1}).limit(5000)
	seqs = [r for r in results if '/OR' not in r['v_gene']['gene']]
	print('Done.\nFound {} sequences\n'.format(len(seqs)))
	return seqs


def aggregate(data, norm=True, sort_by='value'):
	'''
	Counts the number of occurances of each item in 'data'.

	Inputs
	data: a list of values.
	norm: normalize the resulting counts (as percent)
	sort_by: how to sort the retured data. Options are 'value' and 'count'.

	Output
	a non-redundant list of values (from 'data') and a list of counts.
	'''
	vdict = {}
	for d in data:
		vdict[d] = vdict[d] + 1 if d in vdict else 1
	vals = [(k, v) for k, v in vdict.iteritems()]
	if sort_by == 'value':
		vals.sort(key=lambda x: x[0])
	else:
		vals.sort(key=lambda x: x[1])
	xs = [v[0] for v in vals]
	if norm:
		raw_y = [v[1] for v in vals]
		total_y = sum(raw_y)
		ys = [100. * y / total_y for y in raw_y]
	else:
		ys = [v[1] for v in vals]
	return xs, ys




def cdr3_plot(seqs, collection, make_plot):
	if not make_plot:
		return None
	cdr3s = [s['cdr3_len'] for s in seqs if s['cdr3_len'] > 0]
	x, y = aggregate(cdr3s)
	color = sns.hls_palette(7)[4]
	plot_file = '{}_{}_cdr3_lengths.pdf'.format(collection, args.chain)
	make_barplot([str(i) for i in x], y, color, plot_file)


def germline_plot(seqs, gene, collection, level):
	if not level:
		return None
	level = ['fam', 'gene'] if level == 'both' else [level, ]
	for l in level:
		gtype = '{}_gene'.format(gene.lower())
		chain = args.chain[0].upper()
		data = [s[gtype][l].rstrip('D') for s in seqs]
		x, y = aggregate(data)
		x = [i.replace('IG{}{}'.format(chain, gene), '{}{}'.format(gene, chain)) for i in x]
		plot_file = '{}_{}{}_{}.pdf'.format(collection, gene, chain, l)
		colors = get_germline_plot_colors(x, l)
		make_barplot(x, y, colors, plot_file, l)


def get_germline_plot_colors(data, l):
	fams = [d.split('-')[0].replace('/OR', '') for d in data]
	nr_fams = sorted(list(set(fams)))
	num_colors = len(nr_fams)
	rgbs = sns.hls_palette(num_colors)
	rgb_dict = {i[0]: i[1] for i in zip(nr_fams, rgbs)}
	return [rgb_dict[f] for f in fams]



def isotypes():
	pass


def vj_heatmap():
	pass







def make_barplot(x, y, colors, ofile, l=None):
	# set bar locations and width
	ind = np.arange(len(x))
	width = 0.75
	# plot objects
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# axis limits and ticks
	ax.set_ylim(0, 1.05 * max(y))
	ax.set_xlim(-width / 2, len(ind))
	ax.set_xticks(ind + width / 2)
	xtick_names = ax.set_xticklabels(x)
	if l == 'gene':
		plt.setp(xtick_names, rotation=90, fontsize=7)
	# make the plot
	bar = ax.bar(ind, y, width, color=colors)
	plt.savefig(os.path.join(args.output, ofile))
	plt.close()







def print_collection_info(collection):
	print('')
	print('')
	print('-' * 25)
	print(collection)
	print('-' * 25)
	print('')








def main():
	sns.set_style('white')
	db = get_db()
	for collection in get_collections(db):
		print_collection_info(collection)
		seqs = query(db, collection)
		germline_plot(seqs, 'V', collection, level=args.var_plot)
		germline_plot(seqs, 'D', collection, level=args.div_plot)
		germline_plot(seqs, 'J', collection, level=args.join_plot)
		cdr3_plot(seqs, collection, args.cdr3_plot)




if __name__ == '__main__':
	main()
