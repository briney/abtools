#!/usr/bin/python
# filename: tasks.py


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


from abtools.queue.celery import celery
from abtools.utils.alignment import global_alignment


@celery.task
def identity(f, standard, is_aa, debug=False):
	seqs = parse_seqs(f)
	scores = []
	for s in seqs:
		seq_id = '_'.join(s.id.split('_')[:-1])
		germ_identity = float(s.id.split('_')[-1])
		seq = str(s.seq).upper()
		matrix = 'blosum62' if is_aa else None
		aln = nw.global_align(seq,
							  str(standard.seq.upper()),
							  matrix=matrix,
							  gap_open=-15,
							  gap_extend=-2,
							  score_match=1,
							  score_mismatch=0,
							  score_gap_open=0,
							  score_gap_extend=0,
							  aa=is_aa)
		score = aln.score
		# align_matrix = 'nw/blosum62' if is_aa else 'nw/match3mismatch2'
		# alignment = nw.global_align(seq,
		# 							str(standard.seq.upper()),
		# 							matrix=align_matrix,
		# 							gap_open=-15,
		# 							gap_extend=-2)
		# if debug:
		# 	try:
		# 		score = nw.score_alignment(alignment[0],
		# 						   alignment[1],
		# 						   matrix='nw/match1mismatch0',
		# 						   gap_open=0,
		# 						   gap_extend=0)
		# 	except IndexError:
		# 		print(alignment[0])
		# 		print(alignment[1])
		# 		"Unexpected error:", sys.exc_info()[0]
  		# 		raise
		# score = nw.score_alignment(alignment[0],
		# 						   alignment[1],
		# 						   matrix='nw/match1mismatch0',
		# 						   gap_open=0,
		# 						   gap_extend=0)
		aln_length = min(len(aln.aligned_query.strip('-')), len(aln.aligned_target.strip('-')))
		# seq_trunc_align_length = len(alignment[0].rstrip('-').lstrip('-'))
		# standard_trunc_align_length = len(alignment[1].rstrip('-').lstrip('-'))
		# align_length = min(seq_trunc_align_length, standard_trunc_align_length)
		norm_score = 100. * score / aln_length
		norm_score = round(norm_score, 1)
		scores.append((seq_id, norm_score, germ_identity))
	return scores
