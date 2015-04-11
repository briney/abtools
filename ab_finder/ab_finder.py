#!/usr/bin/python
# filename: ab_finder.py


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
import time
import tempfile
import argparse
import subprocess as sp
import multiprocessing as mp

import numpy as np
import pandas as pd

from pymongo import MongoClient

from Bio import SeqIO

import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser("For a MongoDB collection, plots the germline divergence against the sequence identity to a given 'subject' sequence.")
parser.add_argument('-d', '--database', dest='db', required=True,
					help="Name of the MongoDB database to query. Required.")
parser.add_argument('-c', '--collection', dest='collection', default=None,
					help="Name of the MongoDB collection to query. If not provided, all collections in the given database will be processed iteratively.")
parser.add_argument('-o', '--output', dest='output_dir', default=None,
					help="Output directory figure files. If not provided, figures will not be generated. Directory will be created if it does not already exist.")
parser.add_argument('-t', '--temp', dest='temp_dir', required=True,
					help="Directory for temporary storage. Will be created if it does not already exist. Required.")
parser.add_argument('-l', '--log', dest='log', default=None,
					help="The log file, to which the blast_parse log info will be written. Default is stdout.")
parser.add_argument('-C', '--cluster', dest="cluster", default=False, action='store_true',
					help="Use if performing computation on a Celery cluster. If set, input files will be split into many subfiles and passed to a Celery queue. \
					If not set, input files will still be split, but will be distributed to local processors using multiprocessing.")
parser.add_argument('-i', '--ip', dest='ip', default='localhost',
					help="The IP address for the MongoDB server.  Defaults to 'localhost'.")
parser.add_argument('-p', '--port', dest='port', default=27017,
					help="The port for the MongoDB server.  Defaults to '27017'.")
parser.add_argument('-s', '--standard', dest='standard', required=True,
					help='Path to a file containing the standard sequence(s) for which identity/divergence will be calculated, in FASTA format. \
					All sequences in the standard file will iteratively processed. Required')
parser.add_argument('-q', '--chain', dest='chain', default='heavy', choices=['heavy', 'kappa', 'lambda', 'light'],
					help="The chain type of the subject sequence.  Options are 'heavy', 'kappa', 'lambda' and 'light'.  \
					Default is 'heavy'.")
parser.add_argument('-n', '--no_update', dest='update', action='store_false', default=True,
					help="Does not update the MongoDB with ab_compare info. Can save some time if the idenentity calculations aren't needed again.")
parser.add_argument('-N', '--nucleotide', dest='is_aa', action='store_false', default=True,
					help="Use nucleotide sequences for alignment. Default is amino acid sequences. Ensure standard format matches.")
parser.add_argument('-x', '--xmin', dest='x_min', type=int, default=-1,
					help="Minimum X-axis (germline divergence) value for the AbCompare plot. Default is -1.")
parser.add_argument('-X', '--xmax', dest='x_max', type=int, default=35,
					help="Maximum X-axis (germline divergence) value for the AbCompare plot. Default is 35.")
parser.add_argument('-y', '--ymin', dest='y_min', type=int, default=65,
					help="Minimum Y-axis (mAb identity) value for the AbCompare plot. Default is 65.")
parser.add_argument('-Y', '--ymax', dest='y_max', type=int, default=101,
					help="Maximum Y-axis (mAb identity) value for the AbCompare plot. Default is 101.")
parser.add_argument('-g', '--gridsize', dest='gridsize', type=int, default=0,
					help="Gridsize for the AbCompare plot. Default is 36 for amino acid sequences and 50 for nucleotide sequences.")
parser.add_argument('-D', '--debug', dest="debug", action='store_true', default=False,
					help="If set, will write all failed/exception sequences to file and should give more informative errors.")
args = parser.parse_args()



# ================================================
#
#             FILES AND DIRECTORIES
#
# ================================================



def make_directories():
	for d in [args.output_dir, args.temp_dir]:
		if d:
			_make_direc(d)


def _make_direc(d):
	if not os.path.exists(d):
		os.makedirs(d)
	if args.cluster:
		cmd = 'sudo chmod 777 {}'.format(d)
		p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
		stdout, stderr = p.communicate()


def get_standards():
	standards = []
	for s in SeqIO.parse(open(args.standard, 'r'), 'fasta'):
		standards.append(s)
	return standards


def get_chain():
	if args.chain == 'light':
		return ['kappa', 'lambda']
	return [args.chain, ]


def get_sequences(collection, temp_dir, log):
	files = []
	fastas = []
	chunksize = 1000
	seq_counter = 0
	total_seq_counter = 0
	query_results = query(collection)
	iden_field = 'aa_identity' if args.is_aa else 'nt_identity'
	vdj_field = 'vdj_aa' if args.is_aa else 'vdj_nt'
	for seq in query_results:
		fastas.append('>{}_{}\n{}'.format(seq['seq_id'], seq[iden_field]['v'], seq[vdj_field]))
		seq_counter += 1
		total_seq_counter += 1
		if seq_counter == chunksize:
			files.append(write_to_temp_file(fastas, temp_dir))
			fastas = []
			seq_counter = 0
	if fastas:
		files.append(write_to_temp_file(fastas, temp_dir))
	print_query_done()
	return files


def write_to_temp_file(fastas, temp_dir):
	tfile = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
	tfile.write('\n'.join(fastas))
	tfile.close()
	return tfile.name


def clean_up(files):
	for f in files:
		os.unlink(f)



# ================================================
#
#                   MONGO
#
# ================================================



def get_collections():
	if not args.collection:
		subjects = db.collection_names()
		subjects.remove('system.indexes')
		return sorted(subjects)
	return [args.collection, ]


def query(collection):
	coll = db[collection]
	chain = get_chain()
	ensure_index(collection, 'chain')
	print_query_info()
	iden_field = 'aa_identity.v' if args.is_aa else 'nt_identity.v'
	vdj_field = 'vdj_aa' if args.is_aa else 'vdj_nt'
	return coll.find({'chain': {'$in': chain}}, {'_id': 0, 'seq_id': 1, iden_field: 1, vdj_field: 1})


def ensure_index(collection, field):
	print("Indexing '{}' on {}...".format(field, collection))
	coll = db[collection]
	coll.ensure_index(field)


def chunker(l, n):
	'Generator that produces n-length chunks from iterable l.'
	for i in xrange(0, len(l), n):
		yield l[i:i + n]


def update_db(standard, scores, collection):
	print_update_info()
	start = time.time()
	standard = standard.replace('.', '_')
	g = scores.groupby('identity')
	p = mp.Pool(processes=250)
	async_results = []
	groups = regroup(g.groups)
	for group in groups:
		async_results.append(p.apply_async(update, args=(collection, group, standard)))
	monitor_update(async_results)
	p.close()
	p.join()
	run_time = time.time() - start
	print('Updating took {} seconds. ({} sequences per second)'.format(round(run_time, 2), round(len(scores) / run_time, 1)))
	print_done()


def update(collection, data, standard):
	coll = db[collection]
	score = data[0]
	ids = data[1]
	mab_id_field = 'mab_identity_aa' if args.is_aa else 'mab_identity_nt'
	coll.update({'seq_id': {'$in': ids}}, {'$set': {'{}.{}'.format(mab_id_field, standard.lower()): float(score)}}, multi=True)


def monitor_update(results):
	finished = 0
	jobs = len(results)
	while finished < jobs:
		time.sleep(1)
		finished = len([r for r in results if r.ready()])
		update_progress(finished, jobs, sys.stdout)
	sys.stdout.write('\n\n')


def regroup(oldgs):
	newgs = []
	for og in oldgs:
		if len(oldgs[og]) <= 500:
			newgs.append((og, oldgs[og]))
		else:
			for ng in chunker(oldgs[og], 500):
				newgs.append((og, ng))
	return newgs


def remove_padding(collection):
	c = db[collection]
	print_remove_padding()
	c.update({}, {'$unset': {'padding': ''}}, multi=True)




# ================================================
#
#                  FIGURES
#
# ================================================



def make_figure(standard_id, scores, collection):
	print_fig_info()
	fig_file = os.path.join(args.output_dir, '{0}_{1}_{2}.pdf'.format(args.db, collection, standard_id))
	x = list(scores['germ_divergence'].values)
	y = list(scores['identity'].values)
	xy_vals = zip(x, y)
	trunc_xy_vals = [v for v in xy_vals if v[0] <= args.x_max and v[1] >= args.y_min]
	x = [v[0] for v in trunc_xy_vals]
	y = [v[1] for v in trunc_xy_vals]

	# To make sure the gridsize is correct (since it's based on the actual values)
	# I need to add a single value near the max and min of each axis.
	# They're added just outside the visible plot, so there's no effect on the plot.
	x.extend([args.x_min - 1, args.x_max + 1])
	y.extend([args.y_min - 1, args.y_max + 1])

	# plot params
	plt.subplots_adjust(hspace=0.95)
	plt.subplot(111)
	plt.hexbin(x, y, bins='log', cmap=mpl.cm.jet, mincnt=3, gridsize=set_gridsize())
	plt.title(standard_id, fontsize=18)

	# set and label axes
	plt.axis([args.x_min, args.x_max, args.y_min, args.y_max])
	plt.xlabel('Germline divergence')
	plt.ylabel('{0} identity'.format(standard_id))

	# make and label the colorbar
	cb = plt.colorbar()
	cb.set_label('Sequence count (log10)', labelpad=10)

	# save figure and close
	plt.savefig(fig_file)
	plt.close()
	print_done()


def set_gridsize():
	if args.gridsize:
		return args.gridsize
	elif args.is_aa:
		return 36
	return 50




# ================================================
#
#                  PRINTING
#
# ================================================



def print_standards_info(standards):
	print('')
	print('')
	print('Found {} standard sequence(s):'.format(len(standards)))
	print(', '.join([s.id for s in standards]))


def print_collections_info(collections):
	print('')
	print('Found {} collection(s):'.format(len(collections)))
	print(', '.join(collections))


def print_single_standard(standard):
	standard_id_string = 'Standard ID: {}'.format(standard.id)
	print('-' * len(standard_id_string))
	print(standard_id_string)
	print('-' * len(standard_id_string))
	print('')


def print_single_collection(collection):
	collection_string = '      Collection: {}      '.format(collection)
	print('')
	print('')
	print('=' * len(collection_string))
	print(collection_string)
	print('=' * len(collection_string))
	print('')


def print_query_info():
	print('Querying for comparison sequences...', end='')
	sys.stdout.flush()


def print_remove_padding():
	print('Removing MongoDB padding...')


def print_query_done():
	print('Done.')
	print('\n')


def print_fig_info():
	print('Making the identity/divergence figure...', end='')
	sys.stdout.flush()


def print_update_info():
	print('')
	print('Updating the MongoDB database with identity scores:')


def print_done():
	print('Done.')
	print('')




# ================================================
#
#                IDENTITY JOBS
#
# ================================================



def run_jobs(files, standard, log):
	log.write('Running AbCompare...\n')
	if args.cluster:
		return _run_jobs_via_celery(files, standard, log)
	else:
		return _run_jobs_via_multiprocessing(files, standard, log)


def _run_jobs_via_multiprocessing(files, standard, log):
	from utils.identity import identity
	results = []
	if args.debug:
		for f in files:
			results.extend(identity(f, standard, args.is_aa, args.debug))
	else:
		p = mp.Pool()
		async_results = []
		for f in files:
			async_results.append(p.apply_async(identity, (f, standard, args.is_aa)))
		monitor_mp_jobs(async_results, log)
		for a in async_results:
			results.extend(a.get())
		p.close()
		p.join()
	ids = [r[0] for r in results]
	identities = pd.Series([r[1] for r in results], index=ids)
	divergences = pd.Series([100. - r[2] for r in results], index=ids)
	d = {'identity': identities, 'germ_divergence': divergences}
	df = pd.DataFrame(d)
	return df


def monitor_mp_jobs(results, log):
	finished = 0
	jobs = len(results)
	while finished < jobs:
		time.sleep(1)
		ready = [ar for ar in results if ar.ready()]
		finished = len(ready)
		update_progress(finished, jobs, log)
	log.write('\n\n')


def _run_jobs_via_celery(files, standard, log):
	from utils.vdj import run as run_vdj
	async_results = []
	for f in files:
		async_results.append(run_vdj.delay(f, standard, args.is_aa))
	succeeded, failed = monitor_celery_jobs(async_results, log)
	# retry any failed jobs
	if failed:
		retry_results = []
		log.write('{} jobs failed and will be retried:\n'.format(len(failed)))
		files_to_retry = [f for i, f in enumerate(files) if async_results[i].failed()]
		for f in files_to_retry:
			retry_results.append(run_vdj.delay(f, standard))
		retry_succeeded, retry_failed = monitor_celery_jobs(retry_results, log)
		succeeded.extend(retry_succeeded)
	scores = []
	for s in succeeded:
		scores.extend(s.get())
	ids = [r[0] for r in scores]
	identities = pd.Series([r[1] for r in scores], index=ids)
	divergences = pd.Series([r[2] for r in scores], index=ids)
	d = {'identity': identities, 'germ_divergence': divergences}
	df = pd.DataFrame(d)
	return df


def monitor_celery_jobs(results, log):
	finished = 0
	jobs = len(results)
	while finished < jobs:
		time.sleep(1)
		succeeded = [ar for ar in results if ar.successful()]
		failed = [ar for ar in results if ar.failed()]
		finished = len(succeeded) + len(failed)
		update_progress(finished, jobs, log, failed=len(failed))
	log.write('\n\n')
	return succeeded, failed


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





def main():
	log = args.log if args.log else sys.stdout
	make_directories()
	standards = get_standards()
	print_standards_info(standards)
	collections = get_collections()
	print_collections_info(collections)
	for collection in collections:
		indexed = False
		print_single_collection(collection)
		remove_padding(collection)
		seq_files = get_sequences(collection, args.temp_dir, log)
		for standard in standards:
			print_single_standard(standard)
			scores = run_jobs(seq_files, standard, log)
			if args.output_dir:
				make_figure(standard.id, scores, collection)
			if args.update:
				if not indexed:
					ensure_index(collection, 'seq_id')
					indexed = True
				update_db(standard.id, scores, collection)
		clean_up(seq_files)


if __name__ == '__main__':
	conn = MongoClient(args.ip, args.port, max_pool_size=2000)
	db = conn[args.db]
	main()
