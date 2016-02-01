#!/usr/bin/env python
# filename: sequence.py


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


import uuid

from Bio.SeqRecord import SeqRecord


class Sequence(object):
    """
    Container for biological (RNA and DNA) sequences.

    <seq> can be one of several things:
        1) a raw sequence, as a string
        2) an iterable, formatted as (seq_id, sequence)
        3) a dict, containing at least the ID (default key = 'seq_id') and a
            sequence (default key = 'vdj_nt'). Alternate <id_key> and <seq_key>
            can be provided.
        4) a Biopython SeqRecord object
        5) an AbTools Sequence object

    If <seq> is provided as a string, the sequence ID can optionally be
    provided via <id>.  If <seq> is a string and <id> is not provided,
    a random sequence ID will be generated with uuid.uuid4().

    Quality scores can be supplied with <qual> or as part of the SeqRecord object.
    If providing a SeqRecord object with quality scores and quality scores via <qual>,
    the <qual> scores will override the SeqRecord quality scores.

    If <seq> is a dictionary, typically the result of a MongoDB query, the dictionary
    can be accessed directly from the Sequence instance. To retrive the value
    for 'junc_aa' in the instantiating dictionary, you would simply:

        s = Sequence(dict)
        junc = s['junc_aa']

    A note on comparing Sequence objects:
    Sequence objects are equal only if their sequences and IDs are identical.
    This means that two Sequence objects with identical sequences but without
    user-supplied IDs won't be equal, because their IDs will have been randomly
    generated.
    """
    def __init__(self, seq, id=None, qual=None, id_key='seq_id', seq_key='vdj_nt'):
        super(Sequence, self).__init__()
        self._input_sequence = None
        self._input_id = id
        self._input_qual = qual
        self.id = None
        self.sequence = None
        self.qual = None
        self._mongo = None
        self._fasta = None
        self._fastq = None
        self._reverse_complement = None
        self._strand = None
        self.id_key = id_key
        self.seq_key = seq_key

        self._process_input(seq, id, qual)


    def __len__(self):
        return len(self.sequence)

    def __iter__(self):
        return iter(self.sequence)

    def __reversed__(self):
        return ''.join(reversed(self.sequence))

    def __contains__(self, item):
        return item in self.sequence

    def __getitem__(self, key):
        return self.mongo.get(key, None)

    def __setitem__(self, key, val):
        self.mongo[key] = val

    def __eq__(self, other):
        if hasattr(other, 'sequence'):
            return self.sequence == other.sequence
        return False


    @property
    def fasta(self):
        if not self._fasta:
            self._fasta = '>{}\n{}'.format(self.id, self.sequence)
        return self._fasta

    @property
    def fastq(self):
        if self.qual is None:
            self._fastq = None
        else:
            if self._fastq is None:
                self._fastq = '@{}\n{}\n+\n{}'.format(self.id, self.sequence, self.qual)
        return self._fastq

    @property
    def reverse_complement(self):
        if self._reverse_complement is None:
            self._reverse_complement = self._get_reverse_complement()
        return self._reverse_complement

    @property
    def mongo(self):
        if self._mongo is None:
            self._mongo = {}
        return self._mongo

    @mongo.setter
    def mongo(self, val):
        self._mongo = val

    @property
    def strand(self):
        if self._strand is None:
            self._strand = 'plus'
        return self._strand

    @strand.setter
    def strand(self, strand):
        self._strand = strand


    def region(self, start=0, end=None):
        if end is None:
            end = len(self.sequence)
        return '>{}\n{}'.format(self.id, self.sequence[start:end])


    def _get_reverse_complement(self):
        rc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
              'Y': 'R', 'R': 'Y', 'S': 'S', 'W': 'W',
              'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
              'H': 'D', 'V': 'B', 'N': 'N'}
        return ''.join([rc.get(res, res) for res in self.sequence[::-1]])


    def _process_input(self, seq, id, qual):
        if type(seq) in [str, unicode]:
            self.sequence = str(seq).upper()
            self._input_sequence = self.sequence
            if id is None:
                id = uuid.uuid4()
            self.id = id
            self.qual = qual
        elif type(seq) == Sequence:
            self.id = seq.id
            self.sequence = seq.sequence
            self._input_sequence = self.sequence
            self.qual = seq.qual
        elif type(seq) in [list, tuple]:
            self.id = str(seq[0])
            self.sequence = str(seq[1]).upper()
            self._input_sequence = self.sequence
            self.qual = qual
        elif type(seq) == SeqRecord:
            self.id = str(seq.id)
            self.sequence = str(seq.seq).upper()
            self._input_sequence = self.sequence
            if qual:
                self.qual = qual
            elif 'phred_quality' in seq.letter_annotations:
                self.qual = seq.letter_annotations['phred_quality']
            elif 'solexa_quality' in seq.letter_annotations:
                self.qual = seq.letter_annotations['solexa_quality']
            else:
                self.qual = None
        elif type(seq) == dict:
            self.id = seq[self.id_key]
            self.sequence = seq[self.seq_key]
            self._input_sequence = self.sequence
            self.qual = qual
            self._mongo = seq
