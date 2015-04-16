#!/usr/bin/python
# filename: ab_phylogeny.py



###########################################################################
#
# Copyright (c) 2014 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################



import os
import sys
import argparse
import subprocess as sp

import matplotlib as mpl
import seaborn as sns

from Bio import SeqIO

from timepoint import Timepoint
from msa import align



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
parser.add_argument('-s', '--sample_id', dest='sample_id', default=None,
					help="Sample ID. If not provided, the sample ID will be inferred from the input file.")
parser.add_argument('-t', '--timepoints', dest='timepoints', default=None,
					help="Tab-delimited file, of the following format (one per line): \
					TimepointName	TimepointOrder	TimepointColor\
					TimepointName is the name appended to the sequences in the input file.\
					TimepointOrder is an integer that indicates the order in which the timepoints should be sorted.\
					TimepointColor is an RGB or hex value that will be used to color the phylogenetic tree.\
					If mab sequences are provided, the 'mabs' entry will be used to sort/color the mab sequences.\
					If not provided, colors will be automatically selected and timepoints will be determined by a simple \
					sort of the raw timepoint values parsed from the input file.")
parser.add_argument('-a', '--aa', dest='is_aa', default=False, action='store_true',
					help="If used, all input files contain amino acid sequences. \
					Default is nucleotide sequences.")
parser.add_argument('-D', '--delimiter', dest='delimiter', default='_',
					help="Delimiter that separates the timepoint and sequence ID. \
					Cannot use ':' or ';', since these can screw up the tree file. \
					Default is '_'.")
args = parser.parse_args()





# ================================================
#
#                 INPUT PARSING
#
# ================================================



def parse_seqs():
	seqs = parse_input_file()
	seqs += parse_root()
	if args.mabs:
		seqs += parse_mabs()
	timepoints = list(set([s.id.split(args.delimiter)[0] for s in seqs]))
	return seqs, timepoints


def parse_input_file():
	seqs = []
	for seq in SeqIO.parse(open(args.input, 'r'), 'fasta'):
		seqs.append(seq)
	return seqs


def parse_root():
	root = SeqIO.read(open(args.root, 'r'), 'fasta')
	root.id = 'root'
	return [root, ]


def parse_mabs():
	seqs = []
	for seq in SeqIO.parse(open(args.mabs, 'r'), 'fasta'):
		if seq.id.split(args.delimiter)[0] != 'mab':
			seq.id = 'mab{}{}'.format(args.delimiter, seq.id)
		seqs.append(seq)
	return seqs


def parse_timepoints(tps):
	if args.timepoints:
		timepoints = []
		with open(args.timepoints, 'r') as f:
			for line in f:
				name, order, color = line.strip().split('\t')
				if name in tps:
					timepoints.append(Timepoint(name, order, color))
	else:
		colors = sns.hls_palette(len(tps), l=0.5, s=0.9)
		for i, tp in enumerate(sorted(tps)):
			timepoints.append(Timepoint(tp, i + 1, colors[i]))
	return timepoints




def make_msa(seqs):
	import msa
	if args.sample_id:
		sample = args.sample_id
	else:
		sample = os.path.basename(args.input).replace('.fasta', '')
	fasta_file = os.path.join(args.output, '{}.fasta'.format(sample))
	alignment = msa.align(seqs, fasta_file)
	return alignment



def make_tree_file(alignment):
	import tree
	return tree.make_tree(alignment, args.is_aa)






def main():
	seqs, tps = parse_seqs()
	timepoints = parse_timepoints(tps)
	alignment = make_msa(seqs)
	tree_file = make_tree_file(alignment)
	make_figure(tree_file, timepoints)



if __name__ == '__main__':
	main()








