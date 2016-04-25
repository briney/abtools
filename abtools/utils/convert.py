#!/usr/bin/env python
# filename: convert.py


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


from itertools import chain
import os
from zipfile import ZipFile

from Bio import SeqIO

from abtools.pipeline import list_files, make_dir


def abi_to_fasta(input, output):
    '''
    Converts ABI or AB1 files to FASTA format.

    Args:

         input (str): Path to a file or directory containing abi/ab1 files or
            zip archives of abi/ab1 files

        output (str): Path to a directory for the output FASTA files
    '''
    direcs = [input, ]
    # unzip any zip archives
    zip_files = list_files(input, ['zip'])
    if zip_files:
        direcs.extend(_process_zip_files(zip_files))
    # convert files
    for d in direcs:
        files = list_files(d, ['ab1', 'abi'])
        seqs = [SeqIO.read(open(f, 'rb'), 'abi') for f in files]
        # seqs = list(chain.from_iterable(seqs))
        fastas = ['>{}\n{}'.format(s.id, str(s.seq)) for s in seqs]
        ofile = os.path.basename(os.path.normpath(d)) + '.fasta'
        opath = os.path.join(output, ofile)
        open(opath, 'w').write('\n'.join(fastas))


def _process_zip_files(zip_files):
    out_dirs = []
    for z in zip_files:
        out_path = '.'.join(z.split('.')[:-1])
        make_dir(out_path)
        zhandle = ZipFile(z)
        zhandle.extractall(out_path)
        out_dirs.append(out_path)
    return out_dirs
