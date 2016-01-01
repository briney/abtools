#!/usr/bin/python
# filename: identity.py


###########################################################################
#
# Copyright (c) 2015 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################


import sys

import nwalign as nw
from Bio import SeqIO

from queue.celery import celery


# @celery.task
def identity(f, standard, is_aa, debug=False):
	seqs = parse_seqs(f)
	scores = []
	for s in seqs:
		seq_id = '_'.join(s.id.split('_')[:-1])
		germ_identity = float(s.id.split('_')[-1])
		seq = str(s.seq).upper()
		align_matrix = 'nw/blosum62' if is_aa else 'nw/match3mismatch2'
		alignment = nw.global_align(seq,
									str(standard.seq.upper()),
									matrix=align_matrix,
									gap_open=-15,
									gap_extend=-2)
		if debug:
			try:
				score = nw.score_alignment(alignment[0],
								   alignment[1],
								   matrix='nw/match1mismatch0',
								   gap_open=0,
								   gap_extend=0)
			except IndexError:
				print(alignment[0])
				print(alignment[1])
				"Unexpected error:", sys.exc_info()[0]
    			raise
		score = nw.score_alignment(alignment[0],
								   alignment[1],
								   matrix='nw/match1mismatch0',
								   gap_open=0,
								   gap_extend=0)
		seq_trunc_align_length = len(alignment[0].rstrip('-').lstrip('-'))
		standard_trunc_align_length = len(alignment[1].rstrip('-').lstrip('-'))
		align_length = min(seq_trunc_align_length, standard_trunc_align_length)
		norm_score = 100. * score / align_length
		norm_score = round(norm_score, 1)
		scores.append((seq_id, norm_score, germ_identity))
	return scores


def parse_seqs(f):
	return [s for s in SeqIO.parse(open(f, 'r'), 'fasta')]
