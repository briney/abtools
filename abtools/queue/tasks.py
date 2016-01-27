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


from Bio import SeqIO

from abtools.queue.celery import celery
from abtools.alignment import global_alignment
from abtools.sequence import Sequence


@celery.task
def identity(f, standard, is_aa, debug=False):
    seqs = [s for s in SeqIO.parse(open(f, 'r'), 'fasta')]
    scores = []
    for s in seqs:
        seq_id = '_'.join(s.id.split('_')[:-1])
        germ_identity = float(s.id.split('_')[-1])
        seq = str(s.seq).upper()
        matrix = 'blosum62' if is_aa else None
        aln = global_alignment(seq,
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
        aln_length = min(len(aln.aligned_query.strip('-')), len(aln.aligned_target.strip('-')))
        norm_score = 100. * score / aln_length
        norm_score = round(norm_score, 1)
        scores.append((seq_id, norm_score, germ_identity))
    return scores
