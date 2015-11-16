#!/usr/bin/env python
# filename: build_germline_db.py



###########################################################################
#
# Copyright (c) 2015 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################



from __future__ import print_function

import shelve
import argparse

from Bio import SeqIO


parser = argparse.ArgumentParser("")
parser.add_argument('-i', '--in', dest='in_file', required=True,
					help="Input FASTA file containing germline nucleotide sequences. Required.")
parser.add_argument('-d', '--database', dest='db', required=True,
					help="The name of the shelve DB to be built. Required.")
args = parser.parse_args()


def build_db(db, seqs):
	for seq in seqs:
		db[seq.id] = str(seq.seq)
	db.close()


def main():
	db = shelve.open(args.db, flag='c')
	seqs = SeqIO.parse(open(args.in_file, 'r'), 'fasta')
	build_db(db, seqs)


if __name__ == '__main__':
	main()
