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
import traceback

from skbio.alignment import StripedSmithWaterman

import nwalign as nw

from Bio import AlignIO, pairwise2
from Bio.SeqRecord import SeqRecord

from abtools import log
from abtools.pipeline import list_files
from abtools.sequence import Sequence



# -------------------------------------
#
#     MULTIPLE SEQUENCE ALIGNMENT
#
# -------------------------------------



def mafft(sequences=None, alignment_file=None, fasta=None, fmt='fasta', threads=-1, as_file=False,
        print_stdout=False, print_stderr=False):
    '''
    Performs multiple sequence alignment with MAFFT.

    MAFFT is a required dependency.

    Args:

        sequences (list): Sequences to be aligned. ``sequences`` can be one of four things:

            1. a FASTA-formatted string

            2. a list of BioPython ``SeqRecord`` objects

            3. a list of AbTools ``Sequence`` objects

            4. a list of lists/tuples, of the format ``[sequence_id, sequence]``

        alignment_file (str): Path for the output alignment file. If not supplied,
            a name will be generated using ``tempfile.NamedTemporaryFile()``.

        fasta (str): Path to a FASTA-formatted file of sequences. Used as an
            alternative to ``sequences`` when suppling a FASTA file.

        fmt (str): Format of the alignment. Options are 'fasta' and 'clustal'. Default
            is 'fasta'.

        threads (int): Number of threads for MAFFT to use. Default is ``-1``, which
            results in MAFFT using ``multiprocessing.cpu_count()`` threads.

        as_file (bool): If ``True``, returns a path to the alignment file. If ``False``,
            returns a BioPython ``MultipleSeqAlignment`` object (obtained by calling
            ``Bio.AlignIO.read()`` on the alignment file).

    Returns:

        Returns a BioPython ``MultipleSeqAlignment`` object, unless ``as_file`` is ``True``,
            in which case the path to the alignment file is returned.

    '''
    if sequences:
        fasta_string = _get_fasta_string(sequences)
        fasta_file = tempfile.NamedTemporaryFile(delete=False)
        fasta_file.write(fasta_string)
        ffile = fasta_file.name
        fasta_file.close()
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
    if print_stdout:
        print(mafft_cline)
        print(stdout)
    if print_stderr:
        print(stderr)
    os.unlink(ffile)
    if as_file:
        return alignment_file
    if os.stat(alignment_file).st_size == 0:
        return None
    aln = AlignIO.read(open(alignment_file), fmt)
    os.unlink(alignment_file)
    return aln


def muscle(sequences=None, alignment_file=None, fasta=None,
    fmt='fasta', as_file=False, maxiters=None, diags=False,
    gap_open=None, gap_extend=None):
    '''
    Performs multiple sequence alignment with MUSCLE.

    MUSCLE is a required dependency.

    Args:

        sequences (list): Sequences to be aligned. ``sequences`` can be one of four things:

            1. a FASTA-formatted string

            2. a list of BioPython ``SeqRecord`` objects

            3. a list of AbTools ``Sequence`` objects

            4. a list of lists/tuples, of the format ``[sequence_id, sequence]``

        alignment_file (str): Path for the output alignment file. If not supplied,
            a name will be generated using ``tempfile.NamedTemporaryFile()``.

        fasta (str): Path to a FASTA-formatted file of sequences. Used as an
            alternative to ``sequences`` when suppling a FASTA file.

        fmt (str): Format of the alignment. Options are 'fasta' and 'clustal'. Default
            is 'fasta'.

        threads (int): Number of threads for MAFFT to use. Default is ``-1``, which
            results in MAFFT using ``multiprocessing.cpu_count()`` threads.

        as_file (bool): If ``True``, returns a path to the alignment file. If ``False``,
            returns a BioPython ``MultipleSeqAlignment`` object (obtained by calling
            ``Bio.AlignIO.read()`` on the alignment file).

        maxiters (int): Passed directly to MUSCLE using the ``-maxiters`` flag.

        diags (int): Passed directly to MUSCLE using the ``-diags`` flag.

        gap_open (float): Passed directly to MUSCLE using the ``-gapopen`` flag. Ignored
            if ``gap_extend`` is not also provided.

        gap_extend (float): Passed directly to MUSCLE using the ``-gapextend`` flag. Ignored
            if ``gap_open`` is not also provided.

    Returns:

        Returns a BioPython ``MultipleSeqAlignment`` object, unless ``as_file`` is ``True``,
            in which case the path to the alignment file is returned.
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
    if all([gap_open is not None, gap_extend is not None]):
        muscle_cline += ' -gapopen {} -gapextend {}'.format(gap_open, gap_extend)
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
    #     return '\n'.join(['>{}\n{}'.format(seq.id, str(seq.seq).upper()) for seq in sequences])
    # # elif type(sequences[0]) == Sequence:
    # #     return '\n'.join(['>{}\n{}'.format(seq.id, seq.seq) for seq in sequences])
    # elif type(sequences[0]) in [list, tuple]:
    #     return '\n'.join(['>{}\n{}'.format(seq[0], seq[1]) for seq in sequences])



# ----------------------------
#
#     PAIRWISE ALIGNMENT
#
# ----------------------------



def local_alignment(query, target=None, targets=None, match=3, mismatch=-2,
        gap_open=-5, gap_extend=-2, matrix=None, aa=False, gap_open_penalty=None, gap_extend_penalty=None):
    '''
    Striped Smith-Waterman local pairwise alignment.

    Args:

        query: Query sequence. ``query`` can be one of four things:

            1. a nucleotide or amino acid sequence, as a string

            2. a Biopython ``SeqRecord`` object

            3. an AbTools ``Sequence`` object

            4. a list/tuple of the format ``[seq_id, sequence]``

        target: A single target sequence. ``target`` can be anything that
            ``query`` accepts.

        targets (list): A list of target sequences, to be proccssed iteratively.
            Each element in the ``targets`` list can be anything accepted by
            ``query``.

        match (int): Match score. Should be a positive integer. Default is 3.

        mismatch (int): Mismatch score. Should be a negative integer. Default is -2.

        gap_open (int): Penalty for opening gaps. Should be a negative integer.
            Default is -5.

        gap_extend (int): Penalty for extending gaps. Should be a negative
            integer. Default is -2.

        matrix (str, dict): Alignment scoring matrix. Two options for passing the
            alignment matrix:

            - The name of a built-in matrix. Current options are ``blosum62`` and ``pam250``.

            - A nested dictionary, giving an alignment score for each residue pair. Should be formatted
              such that retrieving the alignment score for A and G is accomplished by::

                matrix['A']['G']

        aa (bool): Must be set to ``True`` if aligning amino acid sequences. Default
            is ``False``.

    Returns:

        If a single target sequence is provided (via ``target``), a single ``SSWAlignment``
        object will be returned. If multiple target sequences are supplied (via ``targets``),
        a list of ``SSWAlignment`` objects will be returned.
    '''
    if aa and not matrix:
        err = 'ERROR: You must supply a scoring matrix for amino acid alignments'
        raise RuntimeError(err)
    if not target and not targets:
        err = 'ERROR: You must supply a target sequence (or sequences).'
        raise RuntimeError(err)
    if target:
        targets = [target, ]
    # to maintain backward compatibility with earlier AbTools API
    if gap_open_penalty is not None:
        gap_open = -1 * gap_open_penalty
    if gap_extend_penalty is not None:
        gap_extend = -1 * gap_extend_penalty
    alignments = []
    for t in targets:
        try:
            alignment = SSWAlignment(query=query,
                                     target=t,
                                     match=match,
                                     mismatch=mismatch,
                                     matrix=matrix,
                                     gap_open=-1 * gap_open,
                                     gap_extend=-1 * gap_extend,
                                     aa=aa)
            alignments.append(alignment)
        except IndexError:
            continue
    if len(alignments) == 1:
        return alignments[0]
    return alignments


def local_alignment_biopython(query, target=None, targets=None, match=3, mismatch=-2, matrix=None,
        gap_open=-5, gap_extend=-2, aa=False):
    if not target and not targets:
        err = 'ERROR: You must supply a target sequence (or sequences).'
        raise RuntimeError(err)
    if target:
        targets = [target, ]
    alignments = []
    for t in targets:
        try:
            alignment = alignment = BiopythonAlignment(query=query,
                                                       target=t,
                                                       match=match,
                                                       mismatch=mismatch,
                                                       matrix=matrix,
                                                       gap_open=gap_open,
                                                       gap_extend=gap_extend,
                                                       aa=aa)
            alignments.append(alignment)
        except IndexError:
            continue
    if len(alignments) == 1:
        return alignments[0]
    return alignments


def global_alignment(query, target=None, targets=None, match=3, mismatch=-2, gap_open=-5, gap_extend=-2,
        score_match=None, score_mismatch=None, score_gap_open=None,
        score_gap_extend=None, matrix=None, aa=False):
    '''
    Needleman-Wunch global pairwise alignment.

    With ``global_alignment``, you can score an alignment using different
    paramaters than were used to compute the alignment. This allows you to
    compute pure identity scores (match=1, mismatch=0) on pairs of sequences
    for which those alignment parameters would be unsuitable. For example::

        seq1 = 'ATGCAGC'
        seq2 = 'ATCAAGC'

    using identity scoring params (match=1, all penalties are 0) for both alignment
    and scoring produces the following alignment::

        ATGCA-GC
        || || ||
        AT-CAAGC

    with an alignment score of 6 and an alignment length of 8 (identity = 75%). But
    what if we want to calculate the identity of a gapless alignment? Using::

        global_alignment(seq1, seq2,
                         gap_open=20,
                         score_match=1,
                         score_mismatch=0,
                         score_gap_open=10,
                         score_gap_extend=1)

    we get the following alignment::

        ATGCAGC
        ||  |||
        ATCAAGC

    which has an score of 5 and an alignment length of 7 (identity = 71%). Obviously,
    this is an overly simple example (it would be much easier to force gapless alignment
    by just iterating over each sequence and counting the matches), but there are several
    real-life cases in which different alignment and scoring paramaters are desirable.

    Args:

        query: Query sequence. ``query`` can be one of four things:

            1. a nucleotide or amino acid sequence, as a string

            2. a Biopython ``SeqRecord`` object

            3. an AbTools ``Sequence`` object

            4. a list/tuple of the format ``[seq_id, sequence]``

        target: A single target sequence. ``target`` can be anything that
            ``query`` accepts.

        targets (list): A list of target sequences, to be proccssed iteratively.
            Each element in the ``targets`` list can be anything accepted by
            ``query``.

        match (int): Match score for alignment. Should be a positive integer. Default is 3.

        mismatch (int): Mismatch score for alignment. Should be a negative integer. Default is -2.

        gap_open (int): Penalty for opening gaps in alignment. Should be a negative integer.
            Default is -5.

        gap_extend (int): Penalty for extending gaps in alignment. Should be a negative
            integer. Default is -2.

        score_match (int): Match score for scoring the alignment. Should be a positive integer.
            Default is to use the score from ``match`` or ``matrix``, whichever is provided.

        score_mismatch (int): Mismatch score for scoring the alignment. Should be a negative
            integer. Default is to use the score from ``mismatch`` or ``matrix``, whichever
            is provided.

        score_gap_open (int): Gap open penalty for scoring the alignment. Should be a negative
            integer. Default is to use ``gap_open``.

        score_gap_extend (int): Gap extend penalty for scoring the alignment. Should be a negative
            integer. Default is to use ``gap_extend``.

        matrix (str, dict): Alignment scoring matrix. Two options for passing the alignment matrix:

            - The name of a built-in matrix. Current options are ``blosum62`` and ``pam250``.

            - A nested dictionary, giving an alignment score for each residue pair. Should be
              formatted such that retrieving the alignment score for A and G is accomplished by::

                matrix['A']['G']

        aa (bool): Must be set to ``True`` if aligning amino acid sequences. Default
            is ``False``.

    Returns:

        If a single target sequence is provided (via ``target``), a single ``NWAlignment``
        object will be returned. If multiple target sequences are supplied (via ``targets``),
        a list of ``NWAlignment`` objects will be returned.
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
    if target is not None:
        return alignments[0]
    return alignments


class BaseAlignment(object):
    """
    Base class for local and global pairwise alignments.

    .. note::

        All comparisons between ``BaseAlignments``
        are done on the ``score`` attribute (which must be implemented
        by any classes that subclass ``BaseAlignment``). This was done
        so that sorting alignments like so::

            alignments = [list of alignments]
            alignments.sort(reverse=True)

        results in a sorted list of alignments from the highest alignment
        score to the lowest.

    Attributes:

        query (Sequence): The input query sequence, as an AbTools
            ``Sequence`` object.

        target (Sequence): The input target sequence, as an AbTools
            ``Sequence`` object.

        target_id (str): ID of the target sequence.

        raw_query: The raw query, before conversion to a ``Sequence``.

        raw_target: The raw target, before conversion to a ``Sequence``.

    """
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
        return Sequence(sequence)

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
    Structure for performing and analyzing a Smith-Waterman local alignment.

    .. note:

        Exposed attributes and methods are the same as ``NWAlignment``, so
        local and global alignmnts can be handled in the same way. In fact,
        since comparisons are made based on score, local and global alignments
        can be directly compared with constructions like::

            local_aln == global_aln
            local_aln > global_aln
            alignments = sorted([global_aln, local_aln])

    Attributes:

        alignment_type (str): Is 'local' for all ``SSWAlignment`` objects.

        aligned_query (str): The aligned query sequence (including gaps).

        aligned_target (str): The aligned target sequence (including gaps).

        alignment_midline (str): Midline for the aligned sequences, with ``|`` indicating
          matches and a gap indicating mismatches::

            print(aln.aligned_query)
            print(aln.alignment_midline)
            print(aln.aligned_target)

            # ATGC
            # || |
            # ATCC

        score (int): Alignment score.

        query_begin (int): Position in the raw query sequence at which
            the optimal alignment begins.

        query_end (int): Position in the raw query sequence at which the
            optimal alignment ends.

        target_begin (int): Position in the raw target sequence at which
            the optimal alignment begins.

        target_end (int): Position in the raw target sequence at which the
            optimal alignment ends.
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


class BiopythonAlignment(BaseAlignment):
    def __init__(self, query, target, match=3, mismatch=-2, matrix=None,
            gap_open=5, gap_extend=2, aa=False):
        super(BiopythonAlignment, self).__init__(query, target, matrix,
            match, mismatch, gap_open, gap_extend, aa)

        self.alignment_type = 'local'
        self._aln = self._align()
        aln_query, aln_target, score, begin, end = self._aln
        self.aligned_query = aln_query[begin:end]
        self.aligned_target = aln_target[begin:end]
        self.alignment_midline = self._alignment_midline()
        self.score = score
        self.query_begin = self._get_begin_pos(aln_query, begin)
        self.query_end = self._get_end_pos(aln_query, end)
        self.target_begin = self._get_begin_pos(aln_target, begin)
        self.target_end = self._get_end_pos(aln_target, end)

    def _align(self):
        aln = pairwise2.align.localms(self.query.sequence,
                                      self.target.sequence,
                                      self._match,
                                      self._mismatch,
                                      self._gap_open,
                                      self._gap_extend)
        return aln[0]


    def _get_begin_pos(self, seq, begin):
        dashes = seq.count('-', 0, begin)
        return begin - dashes


    def _get_end_pos(self, seq, end):
        return len(seq[:end].replace('-', ''))


class NWAlignment(BaseAlignment):
    """
    Structure for performing and analyzing a Needleman-Wunch global alignment.

    .. note:

        Exposed attributes and methods are the same as ``SSWAlignment``, so
        local and global alignmnts can be handled in the same way. In fact,
        since comparisons are made based on score, local and global alignments
        can be directly compared with constructions like::

            local_aln == global_aln
            local_aln > global_aln
            alignments = sorted([global_aln, local_aln])

    Attributes:

        alignment_type (str): Is 'global' for all ``NWAlignment`` objects.

        aligned_query (str): The aligned query sequence (including gaps).

        aligned_target (str): The aligned target sequence (including gaps).

        alignment_midline (str): Midline for the aligned sequences, with
          ``|`` indicating matches and a gap indicating mismatches::

            print(aln.aligned_query)
            print(aln.alignment_midline)
            print(aln.aligned_target)

            # ATGC
            # || |
            # ATCC

        score (int): Alignment score.

        query_begin (int): Position in the raw query sequence at which
            the optimal alignment begins.

        query_end (int): Position in the raw query sequence at which the
            optimal alignment ends.

        target_begin (int): Position in the raw target sequence at which
            the optimal alignment begins.

        target_end (int): Position in the raw target sequence at which the
            optimal alignment ends.
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
            if self._matrix.lower() in builtin_names:
                return os.path.join(matrix_dir, self._matrix.lower())
            else:
                err = 'The supplied matrix name ({}) does not exist. '.format(matrix)
                err += 'Built-in matrices are: {}'.format(', '.join(builtins))
                raise RuntimeError(err)
        else:
            self._build_matrix_from_params(match, mismatch, os.path.join(matrix_dir, matrix_name))
            return os.path.join(matrix_dir, matrix_name)

    def _align(self):
        matrix = self._get_matrix_file(match=self._match,
                                       mismatch=self._mismatch,
                                       matrix=self._matrix)
        aln = nw.global_align(self.query.sequence,
                              self.target.sequence,
                              gap_open=self._gap_open,
                              gap_extend=self._gap_extend,
                              matrix=matrix)
        return aln

    def _score_alignment(self):
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

    @staticmethod
    def _build_matrix_from_params(match, mismatch, matrix_file):
        mstring = ' {}'.format(match) if len(str(match)) == 1 else str(match)
        mmstring = ' {}'.format(mismatch) if len(str(mismatch)) == 1 else str(mismatch)
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
