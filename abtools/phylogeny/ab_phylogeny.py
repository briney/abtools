#!/usr/bin/python
# filename: ab_phylogeny.py


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


import os
import subprocess as sp
import sys

import matplotlib as mpl
import seaborn as sns

from Bio import SeqIO

from abtools.ab_phylogeny.utils import msa, tree
from abtools.ab_phylogeny.utils.timepoint import Timepoint



def parse_args():
	import argparse
	parser = argparse.ArgumentParser("")
	parser.add_argument('-i', '--input', dest='input', required=True,
						help="Input file, containing sequences in FASTA format with the timepoint appended to the sequence ID. \
						Required.")
	parser.add_argument('-o', '--output', dest='output', required=True,
						help="Output directory for the FASTA, alignment, tree and figure files. \
						Required")
	parser.add_argument('-r', '--root', dest='root', required=True,
						help="FASTA file containing a single sequence that will be used to root the phylogenetic tree. \
						Required.")
	parser.add_argument('-m', '--mabs', dest='mabs', default=None,
						help="FASTA file containing sequences of monoclonal antibodies of the appropriate chain.")
	parser.add_argument('-a', '--alignment', dest='alignment', default=None,
						help="Provide a pre-aligned dataset. \
						Ensure that the root sequence is labeled as 'root' and that all mAb sequences are pre-pended with 'mab_'. \
						If provided, -i, -r and -m are ignored.")
	parser.add_argument('-n', '--newick', dest='newick', default=None,
						help="Provide a pre-generated Newick tree file. \
						Ensure that the root sequence is labeled as 'root' and that all mAb sequences are pre-pended with 'mab_'. \
						If provided -i, -r, -m and -a are ignored.")
	parser.add_argument('-s', '--sample_id', dest='sample_id', default=None,
						help="Sample ID. If not provided, the sample ID will be inferred from the input file.")
	parser.add_argument('-t', '--timepoints', dest='timepoints', default=None,
						help="Tab-delimited file, of the following format (one per line): \
						TimepointName	TimepointOrder	TimepointColor\
						TimepointName is the name appended to the sequences in the input file.\
						TimepointOrder is an integer that indicates the order in which the timepoints should be sorted.\
						TimepointColor is an RGB or hex value that will be used to color the phylogenetic tree.\
						If mAb sequences are provided, the 'mab' entry will be used to sort/color the mAb sequences.\
						If not provided, colors will be automatically selected and timepoints will be determined by a simple \
						sort of the raw timepoint values parsed from the input file.")
	parser.add_argument('-A', '--aa', dest='is_aa', default=False, action='store_true',
						help="If used, all input files contain amino acid sequences. \
						Default is nucleotide sequences.")
	parser.add_argument('-D', '--delimiter', dest='delimiter', default='_',
						help="Delimiter that separates the timepoint and sequence ID. \
						Cannot use ':', ';' or '=', since these can screw up the tree file. \
						Default is '_'.")
	parser.add_argument('-S', '--scale', dest='scale', default=None,
						help="Scale for the resulting ete2 tree. \
						If not provided, the ete2 default value will be used.")
	parser.add_argument('-b', '--branch_vertical_margin', dest='branch_vertical_margin', default=None,
						help="Branch vertical margin for the resulting ete2 tree. \
						If not provided, the ete2 default value will be used..")
	return parser.parse_args()


class Args(object):
	"""docstring for Args"""
	def __init__(self, input=None, output=None,
				 root=None, mabs=None,
				 alignment=None, newick=None,
				 sample_id=None, is_aa=False,
				 timepoints=None, timepoint_delimiter='_',
				 tree_scale=None, branch_vertical_margin=None):
		super(Args, self).__init__()
		if not all([input, output]):
			print('Input and Output files must be provided.')
			sys.exit(1)
		self.input = input
		self.output = output
		self.root = root
		self.mabs = mabs
		self.alignment = alignment
		self.newick = newick
		self.sample_id = sample_id
		self.timepoints = timepoints
		self.is_aa = is_aa
		self.timepoint_delimiter = timepoint_delimiter



# ================================================
#
#                 INPUT PARSING
#
# ================================================



def parse_seqs(args):
	seqs = parse_input_file(args.input)
	if args.root:
		seqs += parse_root(args.root)
	if args.mabs:
		seqs += parse_mabs(args.mabs, args.delimiter)
	timepoints = list(set([s.id.split(args.delimiter)[0] for s in seqs]))
	return seqs, timepoints


def parse_input_file(input_file):
	seqs = []
	for seq in SeqIO.parse(open(input_file, 'r'), 'fasta'):
		seqs.append(seq)
	return seqs


def parse_root(root_file):
	root = SeqIO.read(open(root_file, 'r'), 'fasta')
	root.id = 'root'
	return [root, ]


def parse_mabs(mabs_file, delimiter):
	seqs = []
	for seq in SeqIO.parse(open(mabs_file, 'r'), 'fasta'):
		if seq.id.split(delimiter)[0] != 'mab':
			seq.id = 'mab{}{}'.format(delimiter, seq.id)
		seqs.append(seq)
	return seqs


def parse_timepoints(tps, args):
	timepoints = []
	if args.timepoints:
		with open(args.timepoints, 'r') as f:
			for line in f:
				name, order, color = line.strip().split('\t')
				if name in tps:
					timepoints.append(Timepoint(name, order, color))
	else:
		colors = sns.hls_palette(len(tps), l=0.5, s=0.9)
		timepoints.append(Timepoint('root', 0, colors[0]))
		for i, tp in enumerate(sorted([t for t in tps if t != 'root'])):
			timepoints.append(Timepoint(tp, i + 1, colors[i + 1]))
	return timepoints


def make_msa(seqs, args):
	if args.sample_id:
		sample = args.sample_id
	else:
		sample = '.'.join(os.path.basename(args.input).split('.')[:-1])
	fasta_file = os.path.join(args.output, '{}.fasta'.format(sample))
	alignment = msa.align(seqs, fasta_file)
	return alignment


def make_tree(alignment, timepoints, args):
	return tree.make_tree(alignment,
						  timepoints,
						  args.delimiter,
						  args.is_aa,
						  args.scale,
						  args.branch_vertical_margin)


def run(**kwargs):
	args = Args(**kwargs)
	main(args)


def main(args):
	seqs, tps = parse_seqs(args)
	timepoints = parse_timepoints(tps, args)
	alignment = make_msa(seqs, args)
	make_tree(alignment, timepoints, args)


if __name__ == '__main__':
	args = parse_args()
	main(args)
