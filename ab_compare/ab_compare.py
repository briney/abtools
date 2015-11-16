#!/usr/bin/python
# filename: ab_compare.py


#
# Copyright (c) 2015 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


from __future__ import print_function

import os
import sys
import math
import uuid
import time
import random
import urllib
import sqlite3
import tempfile
import argparse
import itertools
import subprocess as sp
import multiprocessing as mp
from StringIO import StringIO
from collections import Counter

import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from pymongo import MongoClient, ASCENDING



parser = argparse.ArgumentParser("Compares two or more antibody repertoires.")
parser.add_argument('-o', '--output', dest='output', required=True,
					help="Output directory for figure and data files. Required.")
parser.add_argument('-d', '--database', dest='db', required=True,
					help="Name of the MongoDB database to query. Required.")
parser.add_argument('-1', '--collection1', dest='collection1', default=None,
					help="Name of the first MongoDB collection to query for comparison. \
					If -1 and -2 are not provided, all collections in the given database will be processed iteratively.")
parser.add_argument('-2', '--collection2', dest='collection2', default=None,
					help="Name of the second MongoDB collection to query for comparison. \
					If not provided, the collection provided by -1 will be iteratively compared with all other collections in the database. \
					If both -1 and -2 are not provided, all collections in the given database will be processed iteratively.")
parser.add_argument('--collection_prefix', dest='collection_prefix', default=None,
					help="If supplied, will iteratively process only collections beginning with <collection_prefix>.\
					if '-1' is also provided, it will be iteratively compared with all collections beginning with <collection_prefix>.")
parser.add_argument('-i', '--ip', dest='ip', default='localhost',
					help="IP address for the MongoDB server.  Defaults to 'localhost'.")
parser.add_argument('-u', '--user', dest='user', default=None,
					help="Username for the MongoDB server. Not used if not provided.")
parser.add_argument('-p', '--password', dest='password', default=None,
					help="Password for the MongoDB server. Not used if not provided.")
parser.add_argument('-k', '--chunksize', dest='chunksize', type=int, default=100000,
					help="Number of sequences to use for each similarity calculation. \
					Default is 100,000.")
parser.add_argument('-I', '--iterations', dest='iterations', type=int, default=1000,
					help="Number of iterations of the similarity calculation to perform. \
					Default is 10000.")
parser.add_argument('-s', '--similarity_method', dest='method',
					choices=['marisita-horn', 'kullback-leibler', 'jensen-shannon', 'jaccard', 'bray-curtis', 'renkonen', 'cosine'],
					default='marisita-horn',
					help="The similarity calculation to use. \
					Options are 'marisita-horn', 'kullback-leibler' and 'jensen-shannon'. \
					Note that kullback-leibler is actually a divergence measure (lower values indicate greater similarity) \
					Default is 'marisita-horn'.")
parser.add_argument('-c', '--control', dest='control_similarity', default=False, action='store_true',
					help="Plot the control similarity of the two collections")
parser.add_argument('-C', '--chain', dest='chain', default='heavy', choices=['heavy', 'kappa', 'lambda'],
					help="Select the antibody chain to analyze. \
					Options are 'heavy', 'kappa', and 'lambda'. \
					Default is 'heavy'.")
parser.add_argument('--debug', dest='debug', action='store_true', default=False,
					help="If set, will run in debug mode.")
args = parser.parse_args()



# ================================================
#
#                   MONGO
#
# ================================================



def get_db():
	'''
	Gets a MongoDB database connection object.
	'''
	if args.user and args.password:
		password = urllib.quote_plus(password)
		uri = 'mongodb://{}:{}@{}'.format(args.user, password, args.ip)
		conn = MongoClient(uri)
	else:
		conn = MongoClient(args.ip, 27017)
	return conn[args.db]


def create_index(collection, fields):
	'''
	Creates an index in collection on fields.

	Inputs
	collection: the name of a MongoDB collection
	fields: a list of one or more fields
	'''
	index_fields = [(f, ASCENDING) for f in fields]
	db[collection].create_index(index_fields)


def get_collection_pairs():
	if args.collection1 and args.collection2:
		print_single_pair_info(args.collection1, args.collection2)
		return [(args.collection1, args.collection2), ]
	collections = db.collection_names(include_system_collections=False)
	if args.collection_prefix:
		collections = [c for c in collections if c.startswith(args.collection_prefix)]
	if args.collection1:
		collections.sort()
		pairs = zip([args.collection1] * len(collections), collections)
		print_multiple_pairs_info(pairs)
		return pairs
	pairs = [p for p in itertools.combinations(sorted(collections), 2)]
	print_multiple_pairs_info(pairs)
	return pairs


def index_collections(pairs):
	collections = sorted(list(set([c for tup in pairs for c in tup])))
	fields = ['chain', 'prod']
	print('\n\nCreating indexes...')
	for c in collections:
		print(c)
		create_index(c, fields)


def query(collection):
	print('\nGetting {} sequences from MongoDB...'.format(collection), end='')
	sys.stdout.flush()
	c = db[collection]
	results = c.find({'chain': args.chain, 'prod': 'yes'},
					 {'_id': 0, 'v_gene.gene': 1, 'cdr3_len': 1})
	seqs = [r for r in results if '/OR' not in r['v_gene']['gene']]
	print('Done.\nFound {} sequences'.format(len(seqs)))
	return seqs



# ================================================
#
#                   MUNGING
#
# ================================================



def random_sample_no_replacement(s1, s2, only_continuous, norm):
	chunksize = min(int(0.9 * len(s1)), int(0.9 * len(s2)), args.chunksize)
	c1_indexes = []
	c2_indexes = []
	# while len(c1_indexes) < chunksize:
	# 	index = random.randint(0, len(c1) - 1)
	# 	if index not in c1_indexes:
	# 		c1.indexes.append(index)
	# while len(c2_indexes) < chunksize:
	# 	index = random.randint(0, len(c2) - 1)
	# 	if index not in c2_indexes:
	# 		c2.indexes.append(index)
	# c1 = [s1[i] for i in c1_indexes]
	# c2 = [s2[i] for i in c2_indexes]
	c1 = random.sample(s1, chunksize)
	c2 = random.sample(s2, chunksize)
	agg1 = aggregate(c1)
	agg2 = aggregate(c2)
	df = pd.DataFrame({'s1': agg1, 's2': agg2})
	df.fillna(0, inplace=True)
	if only_continuous:
		df = df[(df['s1'] != 0) & (df['s2'] != 0)]
	if norm:
		return normalize(df['s1'], df['s2'])
	return df['s1'], df['s2']


def get_vgenes(c1, prev_data=None):
	data1 = query(c1)
	v1 = [d['v_gene']['gene'] for d in data1]
	if prev_data:
		size = min(len(prev_data), len(v1))
		v1 = v1[:size]
		prev_data = prev_data[:size]
		return prev_data, v1
	return v1


def aggregate(data):
	'''
	Counts the number of occurances of each item in 'data'.

	Input
	data: a list of values.

	Output
	a dict of bins and counts.
	'''
	vdict = {}
	for d in data:
		vdict[d] = vdict[d] + 1 if d in vdict else 1
	return vdict


def normalize(s1, s2):
	return s1 / np.sum(s1), s2 / np.sum(s2)





# ================================================
#
#                 SIMILARITIES
#
# ================================================



def simdif_method(sample1, sample2):
	'''
	Determines the appropriate similarity/divergence method to use
	based on user options.

	Inputs
	sample1: list of properties (i.e. Var genes) for sample 1
	sample2: list of properties (i.e. Var genes) for sample 2

	Outputs (from appropriate method function)
	median: median of similarity/divergence score distribution
	counts: from binned scores, counts for each bin
	bins: from binned scores, left-most boundary of each bin
	similarities: raw similarity scores.
	'''
	methods = {'marisita-horn': (mh_similarity, False, False),
 		       'kullback-leibler': (kl_divergence, True, True),
 		       'jensen-shannon': (js_similarity, False, True),
 		       'jaccard': (jaccard_similarity, False, False),
 		       'bray-curtis': (bc_similarity, False, True),
 		       'renkonen': (renkonen_similarity, False, True),
 		       'cosine': (cosine_similarity, False, False)}
 	method, continuous, normalize = methods[args.method]
 	s1, s2 = random_sample_no_replacement(sample1,
 										  sample2,
 										  continuous,
 										  normalize)
 	return method(s1, s2)


def calculate_similarities(s1, s2):
	if args.debug:
		similarities = []
		for i in xrange(args.iterations):
			similarities.append(simdif_method(s1, s2))
	else:
		async_results = []
		p = mp.Pool()
		for i in xrange(args.iterations):
			async_results.append(p.apply_async(simdif_method, (s1, s2)))
		monitor_mp_jobs(async_results, sys.stdout)
		similarities = [a.get() for a in async_results]
	median = np.median(similarities)
	counts, bins = np.histogram(similarities)
	return median, counts, bins, similarities


def calculate_control_similarities(s1, s2):
	pool = s1 + s2
	random.shuffle(pool)
	return calculate_similarities(pool, pool)


def mh_similarity(sample1, sample2):
	'''
	Calculates the Marista-Horn similarity for two samples.

	NOTE: sample1 and sample2 should be the same length, and
	the sum of each sample should be greater than 0.

	Inputs
	sample1: list of frequencies
	sample2: list of frequencies

	Output
	Marista-Horn similarity (float, between 0 and 1)
	'''
	X = sum(sample1)
	Y = sum(sample2)
	XY = X * Y
	sumXiYi = 0
	sumXiSq = 0
	sumYiSq = 0
	for x, y in zip(sample1, sample2):
		sumXiYi += x * y
		sumXiSq += x * x
		sumYiSq += y * y
	num = 2 * sumXiYi
	denom = (float(sumXiSq) / (X * X) + float(sumYiSq) / (Y * Y)) * XY
	return 1. * num / denom


def kl_divergence(s1, s2):
	'''
	Calculates the Kullback-Leibler divergence for two samples.

	NOTE: s1 and s2 should be the same length, and the sum of
	each should equal 1. Probabilities should be continuous for
	both s1 and s2.

	Inputs
	sample1: probability distribution
	sample2: probability distribution

	Output
	Kullbeck-Leibler similarity (float)
	'''
	n1, n2 = normalize(s1, s2)
	kl = 0
	for Pi, Qi in zip(n1, n2):
		kl += Pi * math.log(Pi / Qi)
	return kl


def js_similarity(s1, s2):
	'''
	Calculates the Jensen-Shannon similarity for two samples.

	NOTE: s1 and s2 should be the same length, and the sum of
	each should equal 1.

	Inputs
	sample1: probability distribution
	sample2: probability distribution

	Output
	Jensen-Shannon similarity (float, between 0 and 1)
	'''
	n1, n2 = normalize(s1, s2)
	Hsum = np.zeros(len(n1))
	sumH = 0
	for s in [n1, n2]:
		sumH += 0.5 * shannon_entropy(s)
		for i, val in enumerate(s):
			Hsum[i] += 0.5 * val
	return 1 - (shannon_entropy(Hsum) - sumH)


def shannon_entropy(prob_dist):
	'''
	Calculates the Shannon entropy for a single probability distribution.

	Input
	prob_dist: probability distribution, must sum to 1

	Output
	Shannon entropy (float)
	'''
	return -sum([p * math.log(p) for p in prob_dist if p != 0])


def jaccard_similarity(s1, s2):
	'''
	Calculates the Jaccard similarity for two samples.

	NOTE: sample1 and sample2 should be the same length, and
	the sum of each sample should be greater than 0.

	Inputs
	sample1: list of frequencies
	sample2: list of frequencies

	Output
	Jaccard similarity (float, between 0 and 1)
	'''
	num = 0
	denom = 0
	for xy in zip(s1, s2):
		num += min(xy)
		denom += max(xy)
	return 1. * num / denom


def renkonen_similarity(s1, s2):
	'''
	Calculates the Renkonen similarity (also known as the
	percentage similarity) for two samples.

	NOTE: s1 and s2 should be the same length, and
	the sum of each sample should equal 1.

	Inputs
	s1: probability distribution
	s2: probability distribution

	Output
	Renkonen similarity (float, between 0 and 1)
	'''
	return sum(min(vals) for vals in zip(s1, s2))


def bc_similarity(s1, s2):
	'''
	Calculates the Bray-Curtis similarity for two samples.

	NOTE: s1 and s2 should be the same length, and
	the sum of each sample should equal 1.

	Inputs
	s1: probability distribution
	s2: probability distribution

	Output
	Bray-Curtis similarity (float, between 0 and 1)
	'''
	return 2.0 * sum(min(vals) for vals in zip(s1, s2)) / (sum(s1) + sum(s2))


def cosine_similarity(s1, s2):
	'''
	Calculates the cosine (angular) similarity for two samples.

	NOTE: s1 and s2 should be the same length.

	Inputs
	s1: list of frequencies
	s2: list of frequencies

	Output
	Cosine similarity (float, between 0 and 1)
	'''
	num = sum([Pi * Qi for Pi, Qi in zip(s1, s2)])
	denomPi = mp.sqrt(sum([Pi * Pi for Pi in s1])) * mp.sqrt(sum([Qi * Qi for Qi in s2]))
	return float(num) / denom


def sd_similarity(s1, s2):
	'''
	Calculates the Brey-Curtis similarity for two samples.

	NOTE: s1 and s2 should be the same length, and
	the sum of each sample should equal 1..

	Inputs
	s1: list of frequencies
	s2: list of frequencies

	Output
	Brey-Curtis similarity (float, between 0 and 1)
	'''
	num = 0
	denom = 0
	for x, y in zip(s1, s2):
		num += min(x, y) > 0
		denom += min(x) > 0
		denom += min(y) > 0
	return 1. * num / denom


def bin_similarities(similarities):
	return np.histogram(similarities, bins=10)


def monitor_mp_jobs(results, log):
	finished = 0
	jobs = len(results)
	while finished < jobs:
		time.sleep(1)
		ready = [ar for ar in results if ar.ready()]
		finished = len(ready)
		update_progress(finished, jobs, log)
	log.write('\n\n')


def update_progress(finished, jobs, log, failed=None):
	pct = int(100. * finished / jobs)
	ticks = pct / 2
	spaces = 50 - ticks
	if failed:
		prog_bar = '\r({}/{}) |{}{}|  {}% ({}, {})'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct, finished - failed, failed)
	else:
		prog_bar = '\r({}/{}) |{}{}|  {}%'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct)
	sys.stdout.write(prog_bar)
	sys.stdout.flush()





# ================================================
#
#                OUTPUT FILES
#
# ================================================



def make_sim_plot(counts, bins, median, ofile):
	sns.set_style('white')
	bin_mdpt = (bins[1] - bins[0]) / 2
	x = [str(round(b + bin_mdpt, 4)) for b in bins[:-1]]
	y = counts
	# set bar locations and width
	ind = np.arange(len(x))
	width = 0.75
	# plot objects
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# plot aesthetics
	color = sns.hls_palette(7)[4]
	# axis limits, labels and ticks
	ax.set_ylim(0, 1.1 * max(y))
	ax.set_xlim(-width / 2, len(ind))
	ax.set_xticks(ind + width / 2)
	xtick_names = ax.set_xticklabels(x)
	ax.set_xlabel('Similarity index')
	ax.set_ylabel('Frequency')
	# make the plot
	bar = ax.bar(ind, y, width, color='lightgray')
	# plot the median line
	adj_median = ((median - bins[0]) / (bins[-1] - bins[0]) * len(x))
	plt.axvline(x=adj_median, ymin=0, ymax=1, linewidth=1.5, color='r', linestyle='dashed')
	text_xpos = adj_median + 0.2 if adj_median < 6.5 else adj_median - 2.2
	text_ypos = 1.04 * max(y)
	med = round(median, 4)
	med_str = str(med) if len(str(med)) == 6 else str(med) + '0' * (6 - len(str(med)))
	ax.text(text_xpos,
			text_ypos,
			'median: {}'.format(med_str),
			color='r',
			weight='semibold')
	# save the final figure
	plt.savefig(os.path.join(args.output, ofile))
	plt.close()


def write_data(s1, s2, median, counts, bins, similarities, ofile):
	bin_midpoints = [b + (bins[1] - bins[0]) / 2 for b in bins]
	ostring = 'samples: {} and {}'.format(s1, s2)
	ostring += '\n' + ('=' * len(ostring)) + '\n\n'
	ostring += 'method: {} {}\n\n'.format(args.method.title(),
										'divergence' if args.method == 'kullback-leibler' else 'similarity')
	ostring += 'median: {}\n\n'.format(median)
	ostring += 'iterations: {}\n'.format(len(similarities))
	ostring += 'sequences per iteration: {}\n\n'.format(args.chunksize)
	ostring += 'score bin (midpoint)\tcount\n'
	for b, c in zip(bin_midpoints, counts):
		ostring += '{}\t{}\n'.format(b, c)
	ostring += '\nraw scores:\n'
	ostring += '\n'.join([str(s) for s in similarities])
	open(ofile, 'w').write(ostring)


def write_output(s1, s2, median, counts, bins, similarities, control=False):
	ctrl = 'control_' if control else ''
	m = '{}'.format('divergences' if args.method == 'kullback-leibler' else 'similarities')
	fig_file = os.path.join(args.output, '{}_{}_{}_{}{}.pdf'.format(s1, s2, args.method, ctrl, m))
	data_file = os.path.join(args.output, '{}_{}_{}_{}data.txt'.format(s1, s2, args.method, ctrl))
	make_sim_plot(counts, bins, median, fig_file)
	write_data(s1, s2, median, counts, bins, similarities, data_file)


def update_scores(s1, s2, median, scores):
	if s1 not in scores:
		scores[s1] = {}
	scores[s1][s2] = median
	return scores




# ================================================
#
#                  PRINTING
#
# ================================================



def print_method():
	print('\n\nmethod: {} {}'.format(args.method.title(),
								   'divergence' if args.method == 'kullback-leibler' else 'similarity'))


def print_collection_info(collection):
	coll_string = '{0}{1}{0}'.format(' ' * 15, collection)
	print('')
	print('\n')
	print('=' * len(coll_string))
	print(coll_string)
	print('=' * len(coll_string))
	print('')


def print_single_pair_info(p1, p2):
	print('\n')
	print('Found two collections: {} and {}'.format(p1, p2))
	print('\n')
	print('Making a single comparison:')
	print('---------------------------')
	print('{} vs {}'.format(p1, p2))


def print_multiple_pairs_info(pairs):
	print('\n')
	collections = sorted(list(set([c for tup in pairs for c in tup])))
	print('Found {} collections: {}, and {}'.format(len(collections),
													', '.join(collections[:-1]),
													collections[-1]))
	print('\n')
	comp_string = 'Making {} comparisons:'.format(len(pairs))
	print(comp_string)
	print('-' * len(comp_string))
	for p1, p2 in pairs:
		print('{} vs {}'.format(p1, p2))


def print_pair_info(s1, s2):
	print('\n\n')
	compare_string = 'Comparing {} and {}'.format(s1, s2)
	print(compare_string)
	print('-' * len(compare_string))


def print_final_results(scores, control=False):
	'''
	Prints nicely formatted results to sys.stdout.

	Input
	2D dict, containing scores for pairs of samples
	'''
	if not scores:
		return 0
	print('')
	if control:
		print('Control scores:\n')
	else:
		print('Experimental scores:\n')
	for s1 in sorted(scores.keys()):
		print(s1)
		print('-' * len(s1))
		for s2 in sorted(scores[s1].keys()):
			print('{}\t{}'.format(s2, round(scores[s1][s2], 4)))
		print('')







def main():
	print_method()
	pairs = get_collection_pairs()
	index_collections(pairs)
	prev1 = None
	scores = {}
	cscores = {}
	for pair in pairs:
		s1, s2 = pair
		curr1 = s1
		if prev1 != curr1:
			print_collection_info(s1)
			s1_all_vgenes = get_vgenes(s1)
		print_pair_info(s1, s2)
		s1_vgenes, s2_vgenes = get_vgenes(s2, prev_data=s1_all_vgenes)
		print('\nCalculating similarities...')
		median, counts, bins, similarities = calculate_similarities(s1_vgenes, s2_vgenes)
		write_output(s1, s2, median, counts, bins, similarities)
		scores = update_scores(s1, s2, median, scores)
		if args.control_similarity:
			print('\nCalculating control similarities...')
			cmedian, ccounts, cbins, csimilarities = calculate_control_similarities(s1_vgenes, s2_vgenes)
			write_output(s1, s2, cmedian, ccounts, cbins, csimilarities)
			cscores = update_scores(s1, s2, cmedian, cscores)
		prev1 = s1
	print_final_results(scores)
	print_final_results(cscores, control=True)


if __name__ == '__main__':
	db = get_db()
	main()
