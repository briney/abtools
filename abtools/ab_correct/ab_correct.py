#!/usr/bin/env python
# filename: ab_correct.py



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
import math
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

from pymongo import MongoClient

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo


parser = argparse.ArgumentParser("Clusters sequences using either an identity threshold or unique antibody identifiers (UAIDs). \
								  Calculates either the centroid or consensus sequence for each identity/UAID cluster passing a cluster size threshold.")
parser.add_argument('-d', '--database', dest='db', required=True,
					help="Name of the MongoDB database to query. Required")
parser.add_argument('-c', '--collection', dest='collection', default=None,
					help="Name of the MongoDB collection to query. If not provided, all collections in the given database will be processed iteratively.")
parser.add_argument('-o', '--output', dest='output', required=True,
					help="Output directory for the FASTA files. Required")
parser.add_argument('-t', '--temp', dest='temp_dir', required=True,
					help="Directory for temporary files, including the SQLite database. Required")
parser.add_argument('-i', '--ip', dest='ip', default='localhost',
					help="IP address for the MongoDB server.  Defaults to 'localhost'.")
parser.add_argument('-u', '--user', dest='user', default=None,
					help="Username for the MongoDB server. Not used if not provided.")
parser.add_argument('-p', '--password', dest='password', default=None,
					help="Password for the MongoDB server. Not used if not provided.")
parser.add_argument('-m', '--min', dest='min_seqs', default=1, type=int,
					help="Minimum number of sequences for finding centroids from a UAID group.  Defaults to 1.")
parser.add_argument('-U', '--no_uaid', dest='uaid', action='store_false', default=True,
					help="Clusters sequences by identity rather than using universal antibody IDs (UAIDs).")
parser.add_argument('-P', '--parse_uaids', dest='parse_uaids', type=int, default=0,
					help="Length of the UAID to parse, if the UAID was not parsed during AbStar processing. \
					If the '--no_uaid' flag is also used, this option will be ignored. \
					For a UAID of length 20, option should be passed as '--parse_uaids 20'. \
					Default is to not parse UAIDs.")
parser.add_argument('-C', '--consensus', dest='consensus', action='store_true', default=False,
					help="Generates consensus sequences for each UAID or homology cluster. Default is to identify cluster centroids.")
parser.add_argument('-I', '--identity', dest='identity_threshold', default=0.975, type=float,
					help="Identity threshold for sequence clustering. Not used for UAID-based correction. Default is 0.975.")
parser.add_argument('--only-largest-cluster', default=False, action='store_true',
					help="If set while calculating centroids using UAIDs, only the largest centroid for each UAID cluster is retained.")
parser.add_argument('-g', '--germs', dest='germs', default=None,
					help="Path to a FASTA-formatted file of germline V gene sequences. Required if building consensus sequences, not required for centroids.")
parser.add_argument('-D', '--debug', dest='debug', action='store_true', default=False,
					help="If set, will run in debug mode.")
parser.add_argument('-s', '--sleep', dest='sleep', type=int, default=0,
					help="Delay, in minutes, until the script starts executing. Default is 0.")
args = parser.parse_args()





# =========================================
#
#        SEQUENCES AND DATABASES
#
# =========================================



def get_collections():
	if args.collection:
		return [args.collection, ]
	conn = MongoClient(args.ip, 27017)
	db = conn[args.db]
	collections = db.collection_names(include_system_collections=False)
	return sorted(collections)


def get_seqs(collection):
	seqs = query(collection)
	return build_seq_db(seqs)


def query(collection):
	print('Getting sequences from MongoDB...', end='')
	sys.stdout.flush()
	if args.user and args.password:
		password = urllib.quote_plus(password)
		uri = 'mongodb://{}:{}@{}'.format(args.user, password, args.ip)
		conn = MongoClient(uri)
	else:
		conn = MongoClient(args.ip, 27017)
	db = conn[args.db]
	coll = db[collection]
	results = coll.find({'prod': 'yes'}, {'_id': 0, 'v_gene.full': 1, 'seq_id': 1, 'uaid': 1, 'vdj_nt': 1, 'raw_query': 1})
	if args.uaid:
		seqs = []
		for r in results:
			if 'uaid' in r:
				seqs.append((r['seq_id'], r['uaid'], r['vdj_nt'], r['v_gene']['full']))
			elif args.parse_uaids:
				seqs.append((r['seq_id'], r['raw_query'][:args.parse_uaids], r['vdj_nt'], r['v_gene']['full']))
			else:
				err = 'UAID field was not found. '
				err += 'Ensure that UAIDs are being parsed by AbStar, \
				use the -parse_uaids option to parse the UAIDs from the raw query input, \
				or use the -u option for identity-based clustering.'
				raise ValueError(err)
	else:
		seqs = [(r['seq_id'], r['vdj_nt'], r['v_gene']['full']) for r in results]
	print('Done.\nFound {} sequences\n'.format(len(seqs)))
	return seqs


def build_seq_db(seqs):
	print('Building a SQLite database of sequences...', end='')
	sys.stdout.flush()
	db_path = os.path.join(args.temp_dir, 'seq_db')
	conn = sqlite3.connect(db_path)
	c = conn.cursor()
	create_cmd = get_seq_db_creation_cmd()
	insert_cmd = get_seq_db_insert_cmd()
	c.execute('DROP TABLE IF EXISTS seqs')
	c.execute(create_cmd)
	c.executemany(insert_cmd, seqs)
	print('Done')
	print('Indexing the SQLite database...', end='')
	sys.stdout.flush()
	start = time.time()
	c.execute('CREATE INDEX seq_index ON seqs (seq_id)')
	print('Done')
	print('Indexing took {} seconds\n'.format(round(time.time() - start, 2)))
	return c


def get_seq_db_creation_cmd():
	if args.uaid:
		return '''CREATE TABLE seqs (seq_id text, uaid text, vdj_nt text, v_gene text)'''
	return '''CREATE TABLE seqs (seq_id text, vdj_nt text, v_gene text)'''


def get_seq_db_insert_cmd():
	if args.uaid:
		return 'INSERT INTO seqs VALUES (?,?,?,?)'
	return 'INSERT INTO seqs VALUES (?,?,?)'


def remove_sqlite_db():
	db_path = os.path.join(args.temp_dir, 'seq_db')
	os.unlink(db_path)


def parse_germs():
	germ_handle = open(args.germs, 'r')
	germs = {}
	for seq in SeqIO.parse(germ_handle, 'fasta'):
		germs[seq.id] = str(seq.seq).upper()
	return germs





# =========================================
#
#           CD-HIT CLUSTERING
#
# =========================================



def cdhit_clustering(seq_db, uaid=True, centroid=False):
	infile = make_cdhit_input(seq_db, uaid=uaid)
	outfile = os.path.join(args.temp_dir, 'clust')
	logfile = open(os.path.join(args.temp_dir, 'log'), 'a')
	threshold = 1.0 if uaid else args.identity_threshold
	do_cdhit(infile.name, outfile, logfile, threshold)
	if centroid:
		cent_handle = open(outfile, 'r')
		if not uaid:
			clust_handle = open('{}.clstr'.format(outfile), 'r')
			sizes = parse_cluster_sizes(clust_handle)
			seqs = parse_centroids(cent_handle, sizes=sizes)
		else:
			seqs = parse_centroids(cent_handle)
	else:
		clust_handle = open('{}.clstr'.format(outfile), 'r')
		seqs, sizes = parse_clusters(clust_handle, seq_db)
	if not args.debug:
		os.unlink(infile.name)
		os.unlink(os.path.join(args.temp_dir, 'log'))
		os.unlink(outfile)
		os.unlink(outfile + '.clstr')
	if uaid:
		return seqs
	return seqs, sizes


def make_cdhit_input(seq_db, uaid=True):
	infile = tempfile.NamedTemporaryFile(dir=args.temp_dir, delete=False)
	fastas = []
	if uaid:
		seqs = seq_db.execute('''SELECT seqs.seq_id, seqs.uaid FROM seqs''')
	else:
		seqs = seq_db.execute('''SELECT seqs.seq_id, seqs.vdj_nt FROM seqs''')
	for s in seqs:
		fastas.append('>{}\n{}'.format(s[0], s[1]))
	infile.write('\n'.join(fastas))
	infile.close()
	return infile


def do_cdhit(fasta, clust, log, threshold):
	seq_type = 'UAIDs' if args.uaid else 'sequences'
	print('Clustering {} with CD-HIT...'.format(seq_type), end='')
	sys.stdout.flush()
	start_time = time.time()
	cdhit_cmd = 'cd-hit -i {} -o {} -c {} -n 5 -d 0 -T 0 -M 35000'.format(fasta, clust, threshold)
	cluster = sp.Popen(cdhit_cmd, shell=True, stdout=log)
	cluster.communicate()
	print('Done.\nClustering took {} seconds\n'.format(round(time.time() - start_time, 2)))


def parse_centroids(centroid_handle, sizes=None):
	counter = 0
	centroids = []
	for seq in SeqIO.parse(centroid_handle, 'fasta'):
		if sizes:
			size = sizes[counter]
			centroids.append('>{}_{}\n{}'.format(seq.id, size, str(seq.seq)))
		else:
			centroids.append('>{}\n{}'.format(seq.id, str(seq.seq)))
		counter += 1
	return centroids


def parse_clusters(cluster_handle, seq_db):
	print('Parsing CD-HIT cluster file...', end='')
	sys.stdout.flush()
	clusters = [c.split('\n') for c in cluster_handle.read().split('\n>')]
	print('Done.\n{} total clusters identified.\n'.format(len(clusters)))
	print('Retrieving cluster sequences...', end='')
	sys.stdout.flush()
	cluster_seqs = []
	start = time.time()
	lengths = []
	for cluster in clusters:
		lengths.append(len(cluster) - 1)
		if len(cluster) >= args.min_seqs + 1:
			cluster_ids = get_cluster_ids(cluster)
			cluster_seqs.append(get_cluster_seqs(cluster_ids, seq_db))
	print('Done.\n{} clusters meet the minimum size cutoff ({} sequences)'.format(len(cluster_seqs), args.min_seqs))
	print('The average cluster contains {} sequences; the largest contains {} sequences.'.format(round(1. * sum(lengths) / len(lengths), 2), max(lengths)))
	print('Cluster parsing and sequence retrieval took {} seconds\n'.format(round(time.time() - start, 2)))
	return cluster_seqs, lengths


def parse_cluster_sizes(cluster_handle):
	clusters = [c.split('\n') for c in cluster_handle.read().split('\n>')]
	lengths = []
	for cluster in clusters:
		lengths.append(len(cluster) - 1)
	return lengths


def get_cluster_ids(cluster):
	ids = []
	for c in cluster[1:]:
		if c:
			ids.append(c.split()[2][1:-3])
	return ids


def chunker(l, size=900):
	return (l[pos:pos + size] for pos in xrange(0, len(l), size))


def get_cluster_seqs(seq_ids, seq_db):
	seqs = []
	for chunk in chunker(seq_ids):
		seq_chunk = seq_db.execute('''SELECT seqs.seq_id, seqs.vdj_nt
								   FROM seqs
								   WHERE seqs.seq_id IN ({})'''.format(','.join('?' * len(chunk))), chunk)
		seqs.extend(seq_chunk)
	return ['>{}\n{}'.format(s[0], s[1]) for s in seqs]





# =========================================
#
#       UAID CENTROID CLUSTERING
#
# =========================================



def get_uaid_centroids(uaid_clusters):
	print('Calculating centroid sequences with USEARCH:')
	sys.stdout.flush()
	start_time = time.time()
	centroids = []
	singletons = [c[0] for c in uaid_clusters if len(c) == 1]
	for s in singletons:
		seq_id = s.split('\n')[0]
		seq = s.split('\n')[1]
		centroids.append('>{}_{}\n{}'.format(seq_id, 1, seq))
	sizes = [1] * len(centroids)
	clusters = [c for c in uaid_clusters if len(c) > 1]
	if args.debug:
		for cluster in clusters:
			centroid, size = do_usearch_centroid(cluster)
			centroids.extend(centroid)
			sizes.extend(size)
	else:
		p = mp.Pool(maxtasksperchild=100)
		async_results = []
		for c in clusters:
			async_results.append(p.apply_async(do_usearch_centroid, (c, )))
		monitor_mp_jobs(async_results)
		for a in async_results:
			centroid, size = a.get()
			centroids.extend(centroid)
			sizes.extend(size)
		p.close()
		p.join()
		print('Centroids were calculated in {} seconds.'.format(round(time.time() - start_time), 2))
	return centroids, sizes


def do_usearch_centroid(uaid_group_seqs):
	'''
	Clusters sequences at 90% identity using USEARCH.

	Inputs
	uaid_group_seqs: a list of fasta strings corresponding to sequences from a single UAID group.

	Outputs
	A list of fasta strings, containing centroid sequences for each cluster.
	'''
	fasta = tempfile.NamedTemporaryFile(dir=args.temp_dir, prefix='cluster_input_', delete=False)
	results = tempfile.NamedTemporaryFile(dir=args.temp_dir, prefix='results_', delete=False)
	centroids = tempfile.NamedTemporaryFile(dir=args.temp_dir, prefix='centroids_', delete=False)
	fasta.write('\n'.join(uaid_group_seqs))
	fasta.close()
	usearch = ['usearch',
			   '-cluster_fast',
			   fasta.name,
			   '-maxaccepts', '0',
			   '-maxrejects', '0',
			   '-id', '0.9',
			   '-sizeout',
			   '-uc', results.name,
			   '-centroids', centroids.name]
	p = sp.Popen(usearch, stdout=sp.PIPE, stderr=sp.PIPE)
	stdout, stderr = p.communicate()
	centroid_seqs = []
	sizes = []
	for cent in SeqIO.parse(open(centroids.name, 'r'), 'fasta'):
		seq_id = cent.id.split(';')[0]
		size = int(cent.id.split(';')[1].replace('size=', ''))
		centroid_seqs.append('>{}_{}\n{}'.format(seq_id, size, str(cent.seq)))
		sizes.append(size)
	if args.only_largest_cluster:
		cents_plus_sizes = sorted(zip(centroid_seqs, sizes), key=lambda x: x[1], reverse=True)
		centroid_seqs = [cents_plus_sizes[0][0], ]
		sizes = [cents_plus_sizes[0][1], ]
	os.unlink(fasta.name)
	os.unlink(results.name)
	os.unlink(centroids.name)
	return centroid_seqs, sizes





# =========================================
#
#          CONSENSUS SEQUENCES
#
# =========================================



def get_consensus(clusters):
	print('Building consensus sequences...')
	sys.stdout.flush()
	if args.debug:
		consensus_seqs = [calculate_consensus(cluster) for cluster in clusters]
	else:
		p = mp.Pool()
		async_results = [p.apply_async(calculate_consensus, (cluster, )) for cluster in clusters]
		monitor_mp_jobs(async_results)
		results = [a.get() for a in async_results]
		p.close()
		p.join()
	fastas = []
	for r in results:
		seq = r[0]
		size = r[1]
		fastas.append('>{}_{}\n{}'.format(uuid.uuid4(), size, seq))
	return fastas, [r[1] for r in results]


def calculate_consensus(cluster):
	if len(cluster) == 1:
		return (cluster[0].split('\n')[1].upper(), 1)
	fasta_string = consensus_alignment_input(cluster)
	if len(cluster) < 100:
		muscle_cline = 'muscle -clwstrict'
	elif len(cluster) < 1000:
		muscle_cline = 'muscle -clwstrict -maxiters 2'
	else:
		muscle_cline = 'muscle -clwstrict -maxiters 1 -diags'
	muscle = sp.Popen(str(muscle_cline),
					  stdin=sp.PIPE,
					  stdout=sp.PIPE,
					  stderr=sp.PIPE,
					  universal_newlines=True,
					  shell=True)
	alignment = muscle.communicate(input=fasta_string)[0]
	alignment_string = AlignIO.read(StringIO(alignment), 'clustal')
	summary_align = AlignInfo.SummaryInfo(alignment_string)
	consensus = summary_align.gap_consensus(threshold=0.51, ambiguous='n')
	consensus_string = str(consensus).replace('-', '')
	return (consensus_string.upper(), len(cluster))


def consensus_alignment_input(cluster):
	try:
		v_gene = vgene_lookup(cluster)
		germ = '>{}\n{}'.format(v_gene, germs[v_gene])
	except KeyError:
		germ = ''
	cluster.append(germ)
	return '\n'.join(cluster)


def vgene_lookup(cluster):
	v_genes = []
	seq_ids = [seq.split('\n')[0].replace('>', '') for seq in cluster]
	db_path = os.path.join(args.temp_dir, 'seq_db')
	conn = sqlite3.connect(db_path)
	seq_db = conn.cursor()
	for chunk in chunker(seq_ids):
		v_chunk = seq_db.execute('''SELECT seqs.v_gene
								    FROM seqs
								    WHERE seqs.seq_id IN ({})'''.format(','.join('?' * len(chunk))), chunk)
		v_genes.extend(v_chunk)
	v_gene = Counter(v_genes).most_common()[0][0]
	conn.close()
	return v_gene





# =========================================
#
#         PROGRESS AND PRINTING
#
# =========================================



def monitor_mp_jobs(results):
	finished = 0
	jobs = len(results)
	while finished < jobs:
		time.sleep(1)
		ready = [ar for ar in results if ar.ready()]
		finished = len(ready)
		update_progress(finished, jobs, sys.stdout)
	sys.stdout.write('\n')


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


def write_output(collection, fastas, sizes, collection_start_time):
	seq_type = 'consensus' if args.consensus else 'centroid'
	print('\nWriting {} sequences to output file...'.format(seq_type), end='')
	sys.stdout.flush()
	write_fasta_output(collection, fastas)
	write_stats_output(collection, sizes)
	print('Done. \n{} {} sequences were identified.'.format(len(fastas), seq_type))
	print('{} was processed in {} seconds.\n'.format(collection, round(time.time() - collection_start_time, 2)))


def write_fasta_output(collection, fastas):
	seq_type = 'consensus' if args.consensus else 'centroids'
	outfile = os.path.join(args.output, '{}_{}.fasta'.format(collection, seq_type))
	out_handle = open(outfile, 'w')
	out_handle.write('\n'.join(fastas))
	out_handle.close()


def write_stats_output(collection, sizes):
	sizes = [int(s) for s in sizes]
	bin_counts = np.bincount(sizes)[1:]
	bins = range(1, len(bin_counts))
	binned_data = zip(bins, bin_counts)
	bin_string = 'Cluster Size\tCount\n'
	bin_string += '\n'.join(['{}\t{}'.format(b[0], b[1]) for b in binned_data])
	outfile = os.path.join(args.output, '{}_cluster_sizes.txt'.format(collection))
	out_handle = open(outfile, 'w')
	out_handle.write(bin_string)
	out_handle.close()



def print_collection_info(collection):
	print('')
	print('')
	print('-' * 25)
	print(collection)
	print('-' * 25)
	print('')


def countdown():
	start = time.time()
	h = int(args.sleep / 60)
	m = int(args.sleep % 60)
	hz = '0' if len(str(h)) == 1 else ''
	mz = '0' if len(str(m)) == 1 else ''
	print('Countdown: {}{}:{}{}:{}'.format(hz, h, mz, m, '00'), end='')
	sys.stdout.flush()
	done = False
	while not done:
		time.sleep(0.25)
		elapsed = int(time.time() - start)
		if elapsed < args.sleep * 60:
			h = int((args.sleep - elapsed / 60.) / 60)
			m = int((args.sleep - elapsed / 60.) % 60)
			s = int((60 * args.sleep - elapsed) % 60)
			hz = '0' if len(str(h)) == 1 else ''
			mz = '0' if len(str(m)) == 1 else ''
			sz = '0' if len(str(s)) == 1 else ''
			print('\rCountdown: {}{}:{}{}:{}{}'.format(hz, h, mz, m, sz, s), end='')
			sys.stdout.flush()
		else:
			print('\rCountdown: 00:00:00')
			done = True


def main():
	for collection in get_collections():
		collection_start = time.time()
		print_collection_info(collection)
		seq_db = get_seqs(collection)
		if args.uaid:
			uaid_clusters = cdhit_clustering(seq_db)
			if args.consensus:
				sequences, sizes = get_consensus(uaid_clusters)
			else:
				sequences, sizes = get_uaid_centroids(uaid_clusters)
		else:
			if args.consensus:
				seq_clusters, sizes = cdhit_clustering(seq_db, uaid=False)
				sequences, sizes = get_consensus(seq_clusters)
			else:
				filtered_seqs = []
				filtered_sizes = []
				sequences, sizes = cdhit_clustering(seq_db, uaid=False, centroid=True)
				for seq, size in zip(sequences, sizes):
					if size >= args.min_seqs:
						filtered_seqs.append(seq)
						filtered_sizes.append(size)
				sequences = filtered_seqs
				sizes = filtered_sizes
		write_output(collection, sequences, sizes, collection_start)
		remove_sqlite_db()


if __name__ == '__main__':
	if args.sleep:
		countdown()
	if args.consensus:
		germs = parse_germs()
	main()
