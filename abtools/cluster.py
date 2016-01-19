#!/usr/bin/env python
# filename: cluster.py


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

import logging
import os
import sqlite3
import subprocess as sp
import sys
import tempfile
import time

from Bio import SeqIO
from Bio.Align import AlignInfo

from abtools import log
from abtools.alignment import mafft
from abtools.sequence import Sequence


class Cluster(object):
	"""
	Contains a single CD-HIT cluster of sequences.

	Input is a raw cluster (a list of lines from the CD-HIT output
	that contain the information from a single cluster), as well as
	a connection to a SQLite3 sequence database (created by default
	when running cluster.cdhit())

	cluster.sequences returns a list of abtools Sequence objects for
	each sequence in the cluster

	both cluster.centroid and cluster.consensus return abtools Sequence
	objects for the cluster centroid or consensus, respectively.
	"""
	def __init__(self, raw_cluster, seq_db):
		super(Cluster, self).__init__()
		self.raw_cluster = raw_cluster
		self._ids = None
		self._size = None
		self._sequences = None
		self._consensus = None
		self._centroid = None


	@property
	def ids(self):
		if self._ids is None:
			self._ids = self._get_ids()
		return self._ids

	@property
	def size(self):
		if self._size is None:
			self._size = len(self.ids)
		return self._size

	@property
	def sequences(self, db=None):
		if self._sequences is None:
			self._sequences = self._get_sequences
		return self._sequences

	@property
	def consensus(self):
		if self._consensus is None:
			self._consensus = self._make_consensus()
		return self._consensus

	@property
	def centroid(self):
		if self._centroid is None:
			self._centroid = self._get_centroid()
		return self._centroid


	def _get_ids(self):
		ids = []
		for c in self.raw_cluster[1:]:
			if c:
				ids.append(c.split()[2][1:-3])
		return ids

	def _get_sequences(self):
		seqs = []
		for chunk in self._chunker(seq_ids):
			sql_cmd = '''SELECT seqs.id, seqs.sequence
						 FROM seqs
						 WHERE seqs.id IN ({})'''.format(','.join('?' * len(chunk)))
			seq_chunk = seq_db.execute(sql_cmd, chunk)
			seqs.extend(seq_chunk)
		return [Sequence(s) for s in seqs]

	def _get_centroid(self):
		for line in self.raw_cluster:
			if '*' in line:
				centroid_id = line.split()[2][1:-3]
				break
		centroid = [s for s in self.sequences if s.id == centroid_id][0]
		return centroid

	def _make_consensus(self):
		if len(self.sequences) == 1:
			return self.sequences[0]
		aln = mafft(self.sequences)
		summary_align = AlignInfo.SummaryInfo(aln)
		consensus = summary_align.gap_consensus(threshold=0.51, ambiguous='n')
		consensus_string = str(consensus).replace('-', '')
		return Sequence(consensus_string.upper())

	@staticmethod
	def _chunker(l, size=900):
		return (l[pos:pos + size] for pos in xrange(0, len(l), size))



def cdhit(seqs, out_file=None, temp_dir=None, threshold=0.975, make_db=True):
	'''
	Perform CD-HIT clustering on a set of sequences.

	Inputs are an iterable of sequences, which can be in any format that abtools.sequence.Sequence
	can handle.

	Returns the centroid file name and cluster file name (from CD-HIT).
	If ::make_db:: is True (default), a SQLite3 connection and database path are also returned.
	'''
	logger = log.get_logger('cluster')
	start_time = time.time()
	seqs = [Sequence(s) for s in seqs]
	logger.info('CD-HIT: clustering {} seqeunces'.format(len(seqs)))
	if not out_file:
		out_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
		ofile = out_file.name
	else:
		ofile = os.path.expanduser(out_file)
	ifile = _make_cdhit_input(seqs, temp_dir)
	cdhit_cmd = 'cd-hit -i {} -o {} -c {} -n 5 -d 0 -T 0 -M 35000'.format(ifile,
																		  ofile,
																		  threshold)
	cluster = sp.Popen(cdhit_cmd,
					   shell=True,
					   stdout=sp.PIPE,
					   stderr=sp.PIPE)
	stdout, stderr = cluster.communicate()
	os.unlink(ifile)
	logger.info('CD-HIT: clustering took {:.2f} seconds'.format(time.time() - start_time))
	cfile = ofile + '.clstr'
	if make_db:
		logger.info('CD-HIT: building a SQLite3 database')
		seq_db, db_path = _build_seq_db(seqs, direc=temp_dir)
		return ofile, cfile, seq_db, db_path
	return ofile, cfile


def parse_clusters(clust_file, seq_db):
	'''
	Parses clustered sequences.

	Inputs are a CD-HIT cluster file (ends with '.clstr') and a connection to a
	SQLite3 database of sequence IDs and sequences.

	Returns a list of Cluster objects (one per cluster).
	'''
	raw_clusters = [c.split('\n') for c in open(clust_file, 'r').read().split('\n>')]
	return [Cluster(rc, seq_db) for rc in raw_clusters]


def _make_cdhit_input(seqs, temp_dir):
	ifile = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
	fastas = [s.fasta for s in seqs]
	ifile.write('\n'.join(fastas))
	ifile.close()
	return ifile.name


def _build_seq_db(seqs, direc=None):
	'''
	Builds a SQLite3 database of sequences.

	Inputs are a list of Sequence objects and an optional directory to store the database.
	If ::direc:: is not provided, '/tmp' will be used.

	Returns a SQLite3 connection object and the database path.
	'''
	direc = direc if direc else '/tmp'
	db_path = os.path.join(direc, 'seq_db')
	conn = sqlite3.connect(db_path)
	c = conn.cursor()
	create_cmd = '''CREATE TABLE seqs (id text, sequence text)'''
	insert_cmd = 'INSERT INTO seqs VALUES (?,?)'
	c.execute('DROP TABLE IF EXISTS seqs')
	c.execute(create_cmd)
	c.executemany(insert_cmd, [(str(s.id), str(s.sequence)) for s in seqs])
	c.execute('CREATE INDEX seq_index ON seqs (id)')
	return c, db_path
