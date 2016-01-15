#!/usr/bin/env python
# filename: alignment.py


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
from StringIO import StringIO
import subprocess as sp
import tempfile

from skbio.alignment import StripedSmithWaterman

import nwalign as nw

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

from abtools.pipeline import list_files
from abtools.sequence import Sequence



# -------------------------------------
#
#     MULTIPLE SEQUENCE ALIGNMENT
#
# -------------------------------------



def mafft(sequences=None, alignment_file=None, fasta=None, fmt='fasta', threads=-1, as_file=False):
	'''
	Performs multiple sequence alignment with MAFFT

	MAFFT must be installed for this to work

	Input: sequences to be aligned, or a FASTA file of sequences to be aligned
	Returns: A biopython AlignIO object, or path to the alignment file (if ::as_file::)

	::sequences:: can be one of four different things:
		1) a FASTA-formatted string of sequences
		2) a list of Biopython SeqRecord objects
		3) a list of AbTools Sequence objects
		4) a list of lists/tuples, of the format (seq_id, sequence)

	::fasta:: can be used to provide a FASTA-formatted file of sequences
	instead of providing ::sequences::

	::threads:: is the number of threads that MAFFT should use.
	Default is -1, which uses all available cores.

	If a name for the alignment file is not provided (via ::alignment_file::),
	a NamedTemporaryFile will be used

	Options for alignment output format (::fmt::) are "fasta" and "clustal".
	'''
	if sequences:
		fasta_string = _get_fasta_string(sequences)
		fasta_file = tempfile.NamedTemporaryFile(delete=False)
		fasta_file.write(fasta_string)
		ffile = fasta_file.name
	elif fasta:
		ffile = fasta
	if alignment_file is None:
		alignment_file = tempfile.NamedTemporaryFile(delete=False).name
	aln_format = ''
	if fmt == 'clustal':
		aln_format = '--clustalout '
	mafft_cline = 'mafft --thread {} {}{} > {}'.format(threads, aln_format, ffile, alignment_file)
	mafft = sp.Popen(str(mafft_cline),
					 stdout=sp.PIPE,
					 stderr=sp.PIPE,
					 universal_newlines=True,
					 shell=True)
	stdout, stderr = mafft.communicate()
	os.unlink(ffile)
	if as_file:
		return alignment_file
	if os.stat(alignment_file).st_size == 0:
		return None
	aln = AlignIO.read(open(alignment_file), fmt)
	os.unlink(alignment_file)
	return aln


def muscle(sequences=None, alignment_file=None, fasta=None,
	fmt='fasta', as_file=False, maxiters=None, diags=False):
	'''
	Performs multiple sequence alignment with MUSCLE

	MUSCLE must be installed for this to work

	Input: sequences to be aligned, or a FASTA file of sequences to be aligned
	Returns: A biopython AlignIO object, or path to the alignment file (if ::as_file::)

	::sequences:: can be one of four different things:
		1) a FASTA-formatted string of sequences
		2) a list of Biopython SeqRecord objects
		3) a list of AbTools Sequence objects
		4) a list of lists/tuples, of the format (seq_id, sequence)

	::fasta:: can be used to provide a FASTA-formatted file of sequences
	instead of providing ::sequences::

	If a name for the alignment file is not provided (via ::alignment_file::),
	a NamedTemporaryFile will be used

	Options for alignment output format (::fmt::) are "fasta" and "clustal".
	'''
	if sequences:
		fasta_string = _get_fasta_string(sequences)
	elif fasta:
		fasta_string = open(fasta, 'r').read()
	aln_format = ''
	if fmt == 'clustal':
		aln_format = ' -clwstrict'
	muscle_cline = 'muscle{} '.format(aln_format)
	if maxiters is not None:
		muscle_cline += ' -maxiters {}'.format(maxiters)
	if diags:
		muscle_cline += ' -diags'
	muscle = sp.Popen(str(muscle_cline),
					  stdin=sp.PIPE,
					  stdout=sp.PIPE,
					  stderr=sp.PIPE,
					  universal_newlines=True,
					  shell=True)
	alignment = muscle.communicate(input=fasta_string)[0]
	aln = AlignIO.read(StringIO(alignment), fmt)
	if as_file:
		if not alignment_file:
			alignment_file = tempfile.NamedTemporaryFile().name
		AlignIO.write(aln, alignment_file, fmt)
		return alignment_file
	return aln


def consensus(aln, name=None, threshold=0.51, ambiguous='N'):
	summary_align = AlignInfo.SummaryInfo(aln)
	consensus = summary_align.gap_consensus(threshold=threshold, ambiguous=ambiguous)
	if name is None:
		name = uuid.uuid4()
	consensus_string = str(consensus).replace('-', '')
	return (name, consensus_string.upper())


def _get_fasta_string(sequences):
	if type(sequences) == str:
		return sequences
	elif all([type(s) == Sequence for s in sequences]):
		return '\n'.join([s.fasta for s in sequences])
	else:
		return '\n'.join([Sequence(s).fasta for s in sequences])
	# elif type(sequences[0]) == SeqRecord:
	# 	return '\n'.join(['>{}\n{}'.format(seq.id, str(seq.seq).upper()) for seq in sequences])
	# # elif type(sequences[0]) == Sequence:
	# # 	return '\n'.join(['>{}\n{}'.format(seq.id, seq.seq) for seq in sequences])
	# elif type(sequences[0]) in [list, tuple]:
	# 	return '\n'.join(['>{}\n{}'.format(seq[0], seq[1]) for seq in sequences])



# ----------------------------
#
#     PAIRWISE ALIGNMENT
#
# ----------------------------



def local_alignment(query, target=None, targets=None, match=3, mismatch=-2, matrix=None,
	gap_open_penalty=5, gap_extend_penalty=2, aa=False):
	'''
	Wrapper for SSWAlignment, which performs fast Striped Smith-Waterman local pairwise alignment

	Input: query and target sequences
	Returns: a single SSWAlignment object or a list of multiple SSWAlignment objects

	Sequences can be one of four things:
		1) a nucleotide or amino acid sequence, as a string
		2) a Biopython SeqRecord object
		3) a AbTools Sequence object
		4) an iterable of the format (seq_id, sequence)

	::query:: a single sequence

	::target:: can be one of two things:
		1) a single sequence, as a string
		2) an iterable containing one or more sequences as strings

	default scoring parameters:
		match = 3
		mismatch = -2
		gap_open = 5
		gap_extend = 2

	For protein sequences, set ::aa:: to True and optionally provide a scoring matrix.
	::matrix:: can be one of two things:
		1) the name of a built-in matrix (current options are 'blosum62' and 'pam250')
		2) a 2D dict containing match scores for each residue pair (either aa or nt)
	'''
	if aa and not matrix:
		err = 'ERROR: You must supply a scoring matrix for amino acid alignments'
		raise RuntimeError(err)
	if not target and not targets:
		err = 'ERROR: You must supply a target sequence (or sequences).'
		raise RuntimeError(err)
	if target:
		targets = [target, ]
	alignments = []
	for t in targets:
		alignment = SSWAlignment(query=query,
								 target=t,
								 match=match,
								 mismatch=mismatch,
								 matrix=matrix,
								 gap_open=gap_open_penalty,
								 gap_extend=gap_extend_penalty,
								 aa=aa)
		alignments.append(alignment)
	if len(alignments) == 1:
		return alignments[0]
	return alignments


def global_alignment(query, target=None, targets=None, match=3, mismatch=-2, gap_open=-5, gap_extend=-2,
		score_match=None, score_mismatch=None, score_gap_open=None,
		score_gap_extend=None, matrix=None, aa=False):
	'''

	'''
	if not target and not targets:
		err = 'ERROR: You must supply a target sequence (or sequences).'
		raise RuntimeError(err)
	if target:
		targets = [target, ]
	if type(targets) not in (list, tuple):
		err = 'ERROR: ::targets:: requires an iterable (list or tuple).'
		err += 'For a single sequence, use ::target::'
		raise RuntimeError(err)
	alignments = []
	for t in targets:
		alignment = NWAlignment(query=query,
								target=t,
								match=match,
								mismatch=mismatch,
								gap_open=gap_open,
								gap_extend=gap_extend,
								score_match=score_match,
								score_mismatch=score_mismatch,
								score_gap_open=score_gap_open,
								score_gap_extend=score_gap_extend,
								matrix=matrix,
								aa=aa)
		alignments.append(alignment)
	if len(alignments) == 1:
		return alignments[0]
	return alignments


class BaseAlignment(object):
	"""docstring for BaseAlignment"""
	def __init__(self, query, target, matrix,
		match, mismatch, gap_open, gap_extend, aa):
		super(BaseAlignment, self).__init__()
		self.query = self._process_sequence(query, aa=aa)
		self.target = self._process_sequence(target, aa=aa)
		self.raw_query = query
		self.raw_target = target
		self._matrix = matrix
		self._match = int(match)
		self._mismatch = int(mismatch)
		self._gap_open = int(gap_open)
		self._gap_extend = int(gap_extend)
		self._aa = bool(aa)

	def __repr__(self):
		if len(self.aligned_query) > 20:
			qstring = '{}...{}'.format(self.aligned_query[:10], self.aligned_query[-10:])
			mstring = '{}...{}'.format(self.alignment_midline[:10], self.alignment_midline[-10:])
			tstring = '{}...{}'.format(self.aligned_target[:10], self.aligned_target[-10:])
		else:
			qstring = self.aligned_query
			mstring = self.alignment_midline
			tstring = self.aligned_target
		return_string = '\n\n'
		return_string += 'Pairwise Alignment\n'
		return_string += '------------------\n\n'
		return_string += 'query:  {}\n'.format(qstring)
		return_string += '        {}\n'.format(mstring)
		return_string += 'target: {}\n\n'.format(tstring)
		return_string += 'score: {}\n'.format(str(self.score))
		return_string += 'type: {}\n'.format(self.alignment_type)
		return_string += 'length: {}'.format(str(len(self.aligned_query)))
		print(return_string)
		return ''

	def __str__(self):
		return_string = ''
		return_string += '{}\n'.format(self.aligned_query)
		return_string += '{}\n'.format(self.alignment_midline)
		return_string += '{}\n'.format(self.aligned_target)
		return return_string

	def __len__(self):
		return len(self.aligned_query)

	def __eq__(self, other):
		if not hasattr(other, 'score'):
			if type(other) in [int, float]:
				return self.score == other
			return False
		return self.score == other.score

	def __lt__(self, other):
		if not hasattr(other, 'score'):
			if type(other) in [int, float]:
				return self.score == other
			return False
		return self.score < other.score

	def __le__(self, other):
		if not hasattr(other, 'score'):
			if type(other) in [int, float]:
				return self.score == other
			return False
		return self.score <= other.score

	def __gt__(self, other):
		if not hasattr(other, 'score'):
			if type(other) in [int, float]:
				return self.score == other
			return False
		return self.score > other.score

	def __ge__(self, other):
		if not hasattr(other, 'score'):
			if type(other) in [int, float]:
				return self.score == other
			return False
		return self.score >= other.score


	@property
	def target_id(self):
		return self._target_id

	@target_id.setter
	def target_id(self, target_id):
		self._target_id = target_id


	@staticmethod
	def _process_sequence(sequence, aa):
		if type(sequence) == Sequence:
			return sequence
		return Sequence(sequence, aa=aa)

	def _alignment_midline(self):
		midline = ''
		for q, t in zip(self.aligned_query, self.aligned_target):
			if q == t:
				midline += '|'
			else:
				midline += ' '
		return midline


class SSWAlignment(BaseAlignment):
	"""
	Stucture for performing and analyzing a Smith-Waterman alignment
	using the skbio StripedSmithWaterman method.

	Inputs are fairly straight-forward, with a few notable things:
		1) Sequences (::query:: and ::target::) can be one of several things:
			-- just the raw sequence, as a string
			-- an iterable, of the format (sequence ID, sequence)
			-- a Biopython SeqRecord object
			-- a VaxTools Sequence object
		2) ::match:: should be a POSITIVE integer, ::mismatch:: should be a NEGATIVE integer
		3) both gap penalties (gap_open and gap_extend) should be POSITIVE integers
		4) more complex scoring matrices can be specified by name with ::matrix::
			built-in matrices are 'blosum62' and 'pam250'
		5) if you'd like to align using one set of scoring parameters and score using
			a different set, provide all of the 'score_*' parameters

	Exposed properties and convenience methods are the same as NWAlignment objects, so local
	and global alignments can be handled the same way. In fact, since comparisons are made
	based on score, local and global alignments can be directly compared with constructions like
	local_aln == global_aln and local_aln > global_aln.
	"""
	def __init__(self, query, target, match=3, mismatch=-2, matrix=None,
		gap_open=5, gap_extend=2, aa=False):
		super(SSWAlignment, self).__init__(query, target, matrix,
			match, mismatch, gap_open, gap_extend, aa)

		self.alignment_type = 'local'
		self._alignment = self._align()
		self.aligned_query = self._alignment.aligned_query_sequence
		self.aligned_target = self._alignment.aligned_target_sequence
		self.alignment_midline = self._alignment_midline()
		self.score = self._alignment.optimal_alignment_score
		self.query_begin = self._alignment.query_begin
		self.query_end = self._alignment.query_end
		self.target_begin = self._alignment.target_begin
		self.target_end = self._alignment.target_end_optimal

	def _align(self):
		aligner = StripedSmithWaterman(self.query.sequence,
									   match_score=self._match,
									   mismatch_score=self._mismatch,
									   gap_open_penalty=self._gap_open,
									   gap_extend_penalty=self._gap_extend,
									   substitution_matrix=self._matrix,
									   protein=self._aa)
		return aligner(self.target.sequence)


class NWAlignment(BaseAlignment):
	"""
	Stucture for performing and analyzing a Needleman-Wunch alignment
	using the nwalign package.

	Inputs are fairly straight-forward, with a few notable things:
		1) Sequences (::query:: and ::target::) can be one of several things:
			-- just the raw sequence, as a string
			-- an iterable, of the format (sequence ID, sequence)
			-- a Biopython SeqRecord object
			-- a VaxTools Sequence object
		2) ::match:: (and ::score_match::) should be POSITIVE integers. ::mismatch:: (and
			::score_mismatch::) should be NEGATIVE integers or 0.
		2) Gap penalties (gap_open, gap_extend) should be NEGATIVE integers or 0.
		3) more complex scoring matrices can be specified by name with ::matrix::
			built-in matrices are 'blosum62' and 'pam250'
		4) if you'd like to align using one set of scoring parameters and score using
			a different set, provide all of the 'score_*' parameters

	Exposed properties and convenience methods are the same as SWAlignment objects, so local
	and global alignments can be handled the same way. In fact, since comparisons are made
	based on score, local and global alignments can be directly compared with constructions like
	local_aln == global_aln and local_aln > global_aln.
	"""

	def __init__(self, query, target, match=3, mismatch=-2,
		gap_open=-5, gap_extend=-2,
		score_match=None, score_mismatch=None,
		score_gap_open=None, score_gap_extend=None,
		matrix=None, aa=False):
		super(NWAlignment, self).__init__(query, target, matrix,
			match, mismatch, gap_open, gap_extend, aa)
		self.alignment_type = 'global'
		self._score_match = int(score_match) if score_match is not None else None
		self._score_mismatch = int(score_mismatch) if score_mismatch is not None else None
		self._score_gap_open = int(score_gap_open) if score_gap_open is not None else None
		self._score_gap_extend = int(score_gap_extend) if score_gap_extend is not None else None
		self._matrix = matrix
		self._alignment = self._align()
		self.aligned_query = self._alignment[0]
		self.aligned_target = self._alignment[1]
		self.alignment_midline = self._alignment_midline()
		self.score = self._score_alignment()


	def _get_matrix_file(self, match=None, mismatch=None, matrix=None):
		matrix_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utils/matrices')
		builtins = ['blosum62', 'match3mismatch2', 'match1mismatch0']
		if self._matrix is not None:
			matrix_name = self._matrix
		else:
			matrix_name = 'match{}mismatch{}'.format(abs(match), abs(mismatch))
		if matrix_name.lower() in builtins:
			return os.path.join(matrix_dir, matrix_name)
		builtin_names = [os.path.basename(f) for f in list_files(matrix_dir)]
		if self._matrix is not None:
			if self._matrix.lower() in builtin_namess:
				return os.path.join(matrix_dir, self._matrix.lower())
			else:
				err = 'The supplied matrix name ({}) does not exist.'.format(matrix)
				err += 'Built-in matrices are: '.format(', '.join(builtins))
				raise RuntimeError(err)
		else:
			self._build_matrix_from_params(match, mismatch)

	def _align(self):
		# matrix = self._matrix
		# if self._matrix_file is not None:
		# 	matrix = self._matrix_file
		# if self._matrix is None:
		# 	matrix = self._build_matrix(match=self._match,
		# 								mismatch=self._mismatch)
		# elif type(self._matrix) in [str, unicode]:
		# 	self._delete_score_matrix = False
		# 	matrix = self._get_builtin_matrix(self._matrix)
		# elif type(self._matrix) == dict:
		# 	matrix = self._build_matrix(matrix=self._matrix)
		matrix = self._get_matrix_file(match=self._match,
										mismatch=self._mismatch,
										matrix=self._matrix)
		aln = nw.global_align(self.query.sequence,
							  self.target.sequence,
							  gap_open=self._gap_open,
							  gap_extend=self._gap_extend,
							  matrix=matrix)
		# if self._delete_score_matrix:
		# 	os.unlink(matrix)
		return aln

	def _score_alignment(self):
		# matrix = self._matrix
		# if all([self._score_match is not None, self._score_mismatch is not None]):
		# 	matrix = self._build_matrix(match=self._score_match,
		# 								mismatch=self._score_mismatch)
		# elif type(self._matrix) in [str, unicode]:
		# 	self._delete_score_matrix = False
		# 	matrix = self._get_builtin_matrix(self._matrix)
		# elif self._matrix is None:
		# 	matrix = self._build_matrix(match=self._match,
		# 								mismatch=self._mismatch)
		if all([self._score_match is not None, self._score_mismatch is not None]):
			matrix = self._get_matrix_file(match=self._score_match,
										   mismatch=self._score_mismatch)
		elif self._matrix is not None:
			matrix = self._get_matrix_file(matrix=self._matrix)
		else:
			matrix = self._get_matrix_file(match=self._match,
										   mismatch=self._mismatch)
		gap_open = self._score_gap_open if self._score_gap_open is not None else self._gap_open
		gap_extend = self._score_gap_extend if self._score_gap_extend is not None else self._gap_extend
		aln = nw.score_alignment(self.aligned_query,
								self.aligned_target,
								gap_open=gap_open,
								gap_extend=gap_extend,
								matrix=matrix)
		return aln

	# def _build_matrix(self, match=None, mismatch=None, matrix=None):
	# 	if matrix:
	# 		return self._build_matrix_from_dict(matrix)
	# 	return self._build_matrix_from_params(match, mismatch)

	# @staticmethod
	# def _build_matrix_from_dict(matrix):
	# 	matrix_file = tempfile.NamedTemporaryFile(delete=False)
	# 	residues = sorted(matrix.keys())
	# 	header = '   ' + '  '.join(residues)
	# 	matlist = [header, ]
	# 	for r1 in residues:
	# 		resline = [r1, ]
	# 		for r2 in residues:
	# 			s = str(matrix[r1][r2])
	# 			sstring = ' {}'.format(s) if len(s) == 1 else s
	# 		matlist.append(' '.join(resline))
	# 	matrix_file.write('\n'.join(matlist))
	# 	return matrix_file.name

	@staticmethod
	def _build_matrix_from_params(match, mismatch, matrix_file):
		mstring = ' {}'.format(match) if len(str(match)) == 1 else str(match)
		mmstring = ' {}'.format(mismatch) if len(str(mismatch)) == 1 else str(mismatch)
		# matrix_file = tempfile.NamedTemporaryFile(delete=False)
		residues = ['A', 'C', 'D', 'E', 'F',
					'G', 'H', 'I', 'K', 'L',
					'M', 'N', 'P', 'Q', 'R',
					'S', 'T', 'V', 'W', 'Y', '*']
		header = '   ' + '  '.join(residues)
		matlist = [header, ]
		for r1 in residues:
			resline = [r1, ]
			for r2 in residues:
				resline.append(mstring if r1 == r2 else mmstring)
			matlist.append(' '.join(resline))
		open(matrix_file, 'w').write('\n'.join(matlist))
		return matrix_file


	@staticmethod
	def _get_builtin_matrix(matrix_name):
		matrix_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'matrices')
		matrices = [os.path.basename(f) for f in list_files(matrix_dir)]
		if matrix_name.lower() not in matrices:
			err = 'The maxtrix name you provided ({}) is not built-in.'.format(matrix_name)
			err += 'Built in matrices are: {}'.format(', '.join(matrices))
			raise RuntimeError()
		return os.path.join(matrix_dir, matrix_name.lower())
		# matrix_file = tempfile.NamedTemporaryFile(delete=False)
		# matrix_file.write(matrices[matrix_name])
		# return matrix_file.name


# BLOSUM62 = '''#  Matrix made by matblas from blosum62.iij
# #  * column uses minimum score
# #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
# #  Blocks Database = /data/blocks_5.0/blocks.dat
# #  Cluster Percentage: >= 62
# #  Entropy =   0.6979, Expected =  -0.5209
#    A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  * 
# A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
# R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
# N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
# D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
# C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
# Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
# E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
# G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
# H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
# I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
# L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
# K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
# M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
# F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
# P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
# S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
# T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
# W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
# Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
# V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
# B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
# Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
# X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
# * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 '''

# PAM250 = ''
