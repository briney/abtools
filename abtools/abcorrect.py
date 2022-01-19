#!/usr/bin/env python
# filename: _correct.py



###########################################################################
#
# Copyright (c) 2015 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################



# from __future__ import print_function

# from collections import Counter
from datetime import datetime
# from itertools import cycle
# import json
# import math
import multiprocessing as mp
import os
# import paramiko
# import shelve
import sqlite3
# from io import StringIO
import subprocess as sp
import sys
import tempfile
import time
# import urllib.request, urllib.parse, urllib.error
# import uuid

import numpy as np

# from pymongo import MongoClient

# from Bio import SeqIO
# from Bio import AlignIO
# from Bio.Align import AlignInfo

# from abstar.utils.queue.celery import celery

from abutils import Sequence
from abutils.io import read_sequences
from abutils.utils import log
from abutils.utils import mongodb
# from abutils.utils.alignment import muscle
from abutils.utils.cluster import cluster
from abutils.utils.pipeline import make_dir, list_files
from abutils.utils.progbar import progress_bar


def parse_args():
    import argparse
    parser = argparse.ArgumentParser("Clusters sequences using either an identity threshold or unique \
                                      molecular identifiers (UMIs). \
                                      Calculates either the centroid or consensus sequence for each \
                                      identity/UMI cluster that passes a cluster size threshold.")
    parser.add_argument('-d', '--database', dest='db', default=None,
                        help="Name of the MongoDB database to query.")
    parser.add_argument('-c', '--collection', dest='collection', default=None,
                        help="Name of the MongoDB collection to query. \
                        If not provided, all collections in the given database will be processed iteratively.")
    parser.add_argument('--json', dest='json', default=None,
                        help="Input JSON file or directory of JSON files.")
    parser.add_argument('--tabular', dest='tabular', default=None,
                        help="Input TABULAR file or directory of TABULAR files.")
    parser.add_argument('--sep', dest='sep', default='\t',
                        help="Character used to separate fields in TABULAR-formated inputs. Default is '\t'.")
    parser.add_argument('-o', '--output', dest='output', required=True,
                        help="Output directory for the FASTA files. Required")
    parser.add_argument('-l', '--log', dest='log',
                        help="Location for the log file. \
                        If not provided, log will be written to <output>/abcorrect.log")
    parser.add_argument('-t', '--temp', dest='temp_dir', default='/tmp',
                        help="Directory for temporary files, including the SQLite database. \
                        Default is /tmp.")
    parser.add_argument('-i', '--ip', dest='ip', default='localhost',
                        help="IP address for the MongoDB server.  Defaults to 'localhost'.")
    parser.add_argument('--port', dest='port', default=27017, type=int,
                        help="IP address for the MongoDB server.  Defaults to 'localhost'.")
    parser.add_argument('-u', '--user', dest='user', default=None,
                        help="Username for the MongoDB server. Not used if not provided.")
    parser.add_argument('-p', '--password', dest='password', default=None,
                        help="Password for the MongoDB server. Not used if not provided.")
    parser.add_argument('-m', '--min', dest='min_seqs', default=1, type=int,
                        help="Minimum number of sequences for computing a centroid/consensus. \
                        UMI bins or identity-based clusters with fewer sequences will be discarded. Defaults to 1, which \
                        results in centroid/consensus sequences being generated for all UMI bins or identity-based clusters.")
    parser.add_argument('--cluster-by-identity', dest='umi', action='store_false', default=True,
                        help="Clusters sequences by identity rather than using universal antibody IDs (UMIs).")
    parser.add_argument('--parse-umis', dest='parse_umis', action='append', default=[],
                        help="Parameters for UMI parsing. \
                        If UMIs are present in multiple locations (for example, short UMIs at both the start \
                        and end of a sequencing read), multiple parameters may be supplied by providing ``'--parse-umis'`` \
                        multiple times at runtime. \
                        Supplying a single positive integer will result in that number of nucleotides being parsed \
                        from the start (5' end) of the raw sequence. A single negative integer will result in that number of nucleotides \
                        being parsed from the end (3' end) of the raw sequence. \
                        For locations other than the start/end of the sequence, two positions may be supplied (0-based indexing) \
                        separated by a comma (no whitespace) which indicate the start and end of the index position. \
                        If the ``--no-umi`` option is also used, this option will be ignored. \
                        For a UMI of length 20 at the start of the sequence, option should be passed as ``--parse-umis 20``. \
                        For a UMI of length 20 at the end of the sequence, option should be passed as ``--parse-umis -20``. \
                        For a UMI of length 20 that starts at the 5th position in the sequence, option should be passed as \
                        ``--parse-umis 4,24``. \
                        For two UMIs, both length 20, one that starts at the 5th position and one at the 3' end of the sequence, \
                        options should be passed as ``--parse-umis 4,24 --parse-umis -20``. \
                        Default is to not parse UMIs, which results in using pre-parsed UMI data (found in the ``--umi-key`` field \
                        of the input data file).")
    parser.add_argument('--centroid', dest='consensus', action='store_false', default=True,
                        help="Generates centroid sequences for each UMI or identity-based cluster. \
                        Default is to compute consensus sequences for each cluster.")
    parser.add_argument('--identity-threshold', dest='identity_threshold', default=0.975, type=float,
                        help="If performing identity-based clustering, this threshold is used for the clustering. \
                        If performing UMI-based correction, this threshold is used when clustering sequences assigned \
                        to the same UMI bin. Default is 0.975.")
    parser.add_argument('--only-largest-cluster', default=False, action='store_true',
                        help="If set while calculating centroids using UMIs, \
                        only the largest centroid for each UMI cluster is retained.")
    parser.add_argument('--non-redundant', dest='non_redundant', action='store_true', default=False,
                        help='If set, will make a non-redundant FASTA file using Unix sort. \
                        This is much faster than clustering using mmseqs, but can only remove exact duplicates. \
                        Note that this is not exactly equivalent to clustering with a threshold of 1.0, \
                        because sequences of different length that align perfectly for the entirety \
                        of the shorter sequence will be collapsed when using mmseqs but will not be \
                        collapsed by Unix sort.')
    parser.add_argument('--umi_key', dest='umi_key', default='umi',
                        help="The field to be used for the pre-parsed UMI sequence. Default is 'umi'.")                    
    parser.add_argument('--clustering_key', dest='clustering_key', default='sequence',
                        help="The field to be used for clustering. Default is 'sequence'.")
    parser.add_argument('--id_key', dest='id_key', default='sequence_id',
                        help="The field containing the sequence ID. Default is 'sequence_id'.")
    parser.add_argument('--raw_key', dest='raw_key', default='raw_input',
                        help="The field containing the raw input sequence, used for parsing UMIs. \
                        Default is 'raw_input'.")
    parser.add_argument('--output_key', dest='output_key', default='raw_input',
                        help="The field to be used for creating consensus/centroid sequences. \
                        Only required when parsing UMIs. Default is 'raw_input'.")
    parser.add_argument('--no-cluster-sizes', dest='include_cluster_size', action='store_false', default=True,
                        help="If set, the sequence name of the consensus/centroid will not include the cluster size. \
                        Default is to include the cluster size in the consensus/centroid name.")
    parser.add_argument('-D', '--debug', dest='debug', action='store_true', default=False,
                        help="If set, will run in debug mode.")
    return parser


class Args(object):
    def __init__(self, db=None, collection=None, json=None, tabular=None, sep='\t',
                 output=None, log=None, temp_dir=None, include_cluster_size=True,
                 ip='localhost', port=27017, user=None, password=None,
                 min_seqs=1, identity_threshold=0.975,
                 umi=True, parse_umis=None, non_redundant=False,
                 consensus=True, only_largest_cluster=False,
                 id_key='sequence_id', umi_key='umi', clustering_key='sequence',
                 raw_key='raw_input', output_key='raw_input', debug=False):
        super(Args, self).__init__()
        if all([db is None, json is None, tabular is None]):
            print("\nERROR: One of 'db', 'json' or 'tabular' must be supplied.\n")
            sys.exit(1)
        self.db = db
        self.collection = collection
        self.json = json
        self.tabular = tabular
        self.sep = sep
        self.output = output
        self.log = log
        self.temp_dir = temp_dir if temp_dir is not None else '/tmp'
        self.ip = ip
        self.port = int(port)
        self.user = user
        self.password = password
        self.min_seqs = int(min_seqs)
        self.umi = False if non_redundant else umi
        self.parse_umis = parse_umis if parse_umis is not None else []
        self.consensus = consensus
        self.identity_threshold = float(identity_threshold)
        self.include_cluster_size = include_cluster_size
        self.only_largest_cluster = only_largest_cluster
        self.non_redundant = non_redundant
        self.id_key = id_key
        self.umi_key = umi_key
        self.clustering_key = clustering_key
        self.raw_key = raw_key
        self.output_key = output_key
        self.debug = debug



# =========================================
#
#        SEQUENCES AND DATABASES
#
# =========================================


def get_seqs(input, args):
    # read JSON data
    if args.json is not None:
        logger.info('reading input file...')
        seqs = read_sequences(file=input, format='json',
                              id_key=args.id_key,
                              sequence_key=args.clustering_key)
    # read tabular data
    elif args.tabular is not None:
        logger.info('reading input file...')
        seqs = read_sequences(file=input, format='tabular', sep=args.sep,
                              id_key=args.id_key,
                              sequence_key=args.clustering_key)
    # load MongoDB data
    elif args.db is not None:
        logger.info('querying MongoDB...')
        mongodb_kwargs = {'ip': args.ip,
                          'port': args.port,
                          'user': args.user,
                          'password': args.password}
        seqs = read_sequences(format='mongodb',
                              mongodb_kwargs=mongodb_kwargs,
                              db=args.db, collection=input,
                              id_key=args.id_key,
                              sequence_key=args.clustering_key)
    else:
        seqs = []
    # make sure all of the necessary keys are present in the annotated data
    required_keys = [args.id_key, args.clustering_key, args.output_key]
    if args.parse_umis:
        required_keys.append(args.raw_key)
    check_sequence_keys(seqs, required_keys)
    logger.info(f'found {len(seqs)} input sequences')
    # parse UMIs
    if args.parse_umis:
        logger.info('parsing UMIs...')
        for s in seqs:
            s[args.umi_key] = parse_umi(s[args.raw_key], args)
    return seqs


def check_sequence_keys(seqs, keys):
    for k in keys:
        missing = 0
        for s in seqs:
            if k not in s.annotations:
                missing += 1
        missing_frac = missing / len(seqs)
        if missing_frac > 0.1:
            error = f'ERROR: The "{k}" field was missing in {round(missing_frac * 100, 2)}% of the input sequences. '
            error += 'Please ensure that the id, clustering, raw and output keys are present in the input data.'
            print(f'\n{error}\n')
            sys.exit()


def parse_umi(raw_seq, args):
    umi = ''
    for umi_slice in args.parse_umis:
        if len(umi_slice.split(',')) == 2:
            pos = [int(s) for s in umi_slice.strip().split(',')]
            start = min(pos)
            end = max(pos)
            umi += raw_seq[start:end]
        else:
            pos = int(umi_slice)
            if pos >= 0:
                umi += raw_seq[:pos]
            else:
                umi += raw_seq[pos:]
    return umi


# def query(db, collection, args):
#     # parse MINIMAL file...
#     if args.tabular:
#         logger.info('Reading MINIMAL file...')
#         results = []
#         with open(collection) as f:
#             for i, line in enumerate(f):
#                 # parse header file
#                 if i == 0:
#                     header = line.strip().split(',')
#                     seq_id_index = header.index('seq_id')
#                     uid_index = header.index('uid')
#                     v_gene_index = header.index('v_gene')
#                     clustering_index = header.index(args.clustering_key)
#                     output_index = header.index(args.output_field)
#                     raw_index = header.index('raw_input')
#                 # parse sequence data
#                 else:
#                     try:
#                         l = line.strip().split(',')
#                         # identify the UID (either by parsing or by field)
#                         if args.parse_uaids:
#                             uid = parse_uid(l[raw_index], args)
#                         elif not args.uaid:
#                             uid = ''
#                         else:
#                             uid = l[uid_index]
#                         d = {'seq_id': l[seq_id_index],
#                              'uaid': uid,
#                              args.clustering_key: l[clustering_index],
#                              args.output_field: l[output_index],
#                              'v_gene': {'full': l[v_gene_index]},
#                              'raw_query': l[raw_index]}
#                         results.append(d)
#                     except IndexError:
#                         continue
#     # parse JSON file...
#     elif db is None:
#         logger.info('Reading JSON file...')
#         results = []
#         _results = []
#         with open(collection) as f:
#             for line in f:
#                 if line.strip():
#                     _results.append(json.loads(line.strip()))
#         for r in _results:
#             raw_field = 'raw_query' if 'raw_query' in r else 'raw_input'
#             try:
#                 d = {'seq_id': r['seq_id'],
#                      args.clustering_key: r[args.clustering_key],
#                      args.output_field: r[args.output_field],
#                      'raw_query': r[raw_field],
#                      'v_gene': {'full': r['v_gene']['full']}}
#                 if args.uaid and not args.non_redundant:
#                     # identify the UID (either by parsing or by field)
#                     uid_field = 'uid' if 'uid' in r else 'uaid'
#                     if args.parse_uaids:
#                         d['uaid'] = parse_uid(r[raw_field], args)
#                     elif uid_field in r:
#                         d['uaid'] = r[uid_field]
#                     else:
#                         err = 'ERROR: UAID field was not found. '
#                         err += 'Ensure that UAIDs were parsed by AbStar, '
#                         err += 'use the -parse_uaids option to parse the UAIDs from the raw query input, '
#                         err += 'or use the -u option for identity-based clustering.'
#                         raise ValueError(err)
#                 results.append(d)
#             except KeyError:
#                 continue
#     # ...or get sequences from a MongoDB database
#     else:
#         logger.info('Getting sequences from MongoDB...')
#         coll = db[collection]
#         match = {'prod': 'yes'}
#         project = {'_id': 0, 'seq_id': 1, args.clustering_key: 1, args.output_field: 1, 'raw_query': 1, 'raw_input': 1, 'v_gene.full': 1}
#         if args.uaid is not None:
#             project['uaid'] = 1
#             project['uid'] = 1
#         results = coll.find(match, project)
#     if args.non_redundant:
#         return results
#     elif args.uaid:
#         seqs = []
#         for r in results:
#             # raw_field and uid_field are necessary to maintain compatibility with legacy AbStar versions
#             raw_field = 'raw_query' if 'raw_query' in r else 'raw_input'
#             uid_field = 'uid' if 'uid' in r else 'uaid'
#             if args.parse_uaids:
#                 uid = parse_uid(r[raw_field], args)
#             else:
#                 try:
#                     uid = r[uid_field, '']
#                 except:
#                     err = 'ERROR: UAID field was not found. '
#                     err += 'Ensure that UAIDs were parsed by AbStar, '
#                     err += 'use the -parse_uaids option to parse the UAIDs from the raw query input, '
#                     err += 'or use the -u option for identity-based clustering.'
#                     raise ValueError(err)
#             seqs.append((r['seq_id'], uid, r[args.clustering_key], r[args.output_field], r[raw_field], r['v_gene']['full']))
#     else:
#         seqs = [(r['seq_id'], r[seq_field], r[args.clustering_key], r[args.output_field], r[raw_field], r['v_gene']['full']) for r in results]
#     logger.info('Found {} sequences\n'.format(len(seqs)))
#     return seqs


def build_seq_db(seqs, args):
    logger.info('building a SQLite database of sequences...')
    # get DB data
    keys = list(set([args.id_key, args.umi_key, args.clustering_key, args.output_key, args.raw_key]))
    db_data = get_sequence_db_data(seqs, keys)
    db_path = os.path.join(args.temp_dir, 'seq_db')
    # connect to the DB
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    # create a SQLite table
    c.execute('DROP TABLE IF EXISTS seqs')
    key_string = ', '.join([f'{k} text' for k in keys])
    create_cmd = f'CREATE TABLE seqs ({key_string})'
    c.execute(create_cmd)
    # insert data
    insert_cmd = f"INSERT INTO seqs VALUES ({','.join(['?'] * len(keys))})"
    c.executemany(insert_cmd, db_data)
    # index on the id key
    logger.info('indexing the SQLite database...')
    index_cmd = f'CREATE INDEX seq_index ON seqs ({args.id_key})'
    c.execute(index_cmd)
    return c


def get_sequence_db_data(seqs, keys):
    seq_data = []
    for s in seqs:
        d = [s.get(k, '') for k in keys]
        seq_data.append(d)
    return seq_data


def remove_sqlite_db(args):
    db_path = os.path.join(args.temp_dir, 'seq_db')
    os.unlink(db_path)


# def get_seq_db_creation_cmd(args):
#     if args.uaid:
#         return '''CREATE TABLE seqs (seq_id text, uaid text, clustering_seq text, output_seq text, raw text)'''
#     return '''CREATE TABLE seqs (seq_id text, clustering_seq text, output_seq text, raw text, v_gene text)'''


# def get_seq_db_insert_cmd(args):
#     if args.uaid:
#         return 'INSERT INTO seqs VALUES (?,?,?,?,?,?)'
#     return 'INSERT INTO seqs VALUES (?,?,?,?,?)'


def chunker(l, size=900):
    return (l[pos:pos + size] for pos in range(0, len(l), size))


def get_clustering_seqs(seq_db, args):
    query = f'''SELECT seqs.{args.id_key}, seqs.{args.clustering_key} FROM seqs'''
    seqs = seq_db.execute(query).fetchall()


def get_clustering_seqs_by_id(seq_ids, seq_db, args):
    seqs = []
    for chunk in chunker(seq_ids):
        query = f'''SELECT seqs.{args.id_key}, seqs.{args.clustering_key}
                    FROM seqs
                    WHERE seqs.{args.id_key} in ({','.join('?' * len(chunk))})'''
        seq_chunk = seq_db.execute(query, chunk).fetchall()
        seqs.extend(seq_chunk)
    return [Sequence(s[1], id=s[0]) for s in seqs]


def get_output_seqs_by_id(seq_ids, seq_db, args):
    seqs = []
    for chunk in chunker(seq_ids):
        query = f'''SELECT seqs.{args.id_key}, seqs.{args.output_key}
                    FROM seqs
                    WHERE seqs.{args.id_key} in ({','.join('?' * len(chunk))})'''
        seq_chunk = seq_db.execute(query, chunk).fetchall()
        seqs.extend(seq_chunk)
    return [Sequence(s[1], id=s[0]) for s in seqs]


# def parse_germs(germ_file):
#     germ_handle = open(germ_file, 'r')
#     germs = {}
#     for seq in SeqIO.parse(germ_handle, 'fasta'):
#         germs[seq.id] = str(seq.seq).upper()
#     return germs





# =========================================
#
#             NON-REDUNDANT
#
# =========================================


def sort_unique(seqs, args):
    sort_file, input_count = make_sort_unique_input(seqs, args)
    # logger.info(f'Found {input_count} sequences')
    unique_file = tempfile.NamedTemporaryFile(dir=args.temp_dir, delete=False)
    sort_unique_cmd = f'sort -k2,2 -u -o {unique_file.name} {sort_file}'
    logger.info('running sort/uniq...')
    p = sp.Popen(sort_unique_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    p.communicate()
    os.unlink(sort_file)
    return unique_file.name


def make_sort_unique_input(seqs, args):
    sort_file = tempfile.NamedTemporaryFile(dir=args.temp_dir, delete=False)
    sort_file.close()
    count = 0
    with open(sort_file.name, 'a') as f:
        for s in seqs:
            f.write('{} {} {}\n'.format(s[args.id_key], s[args.clustering_key], s[args.raw_key]))
            count += 1
    return sort_file.name, count


# def parse_unique_file(unique_file):
#     with open(unique_file) as f:
#         for line in f:
#             if len(line.strip().split()) < 3:
#                 continue
#             seq_id, vdj, raw = line.strip().split()
#             fastas.append('>{}\n{}'.format(seq_id, raw))
#     return fastas




# =========================================
#
#              CLUSTERING
#
# =========================================



# def cdhit_clustering(seq_db, args, uaid=True, centroid=False):
#     infile = make_cdhit_input(seq_db, args, uaid=uaid)
#     outfile = os.path.join(args.temp_dir, 'clust')
#     logfile = open(os.path.join(args.temp_dir, 'log'), 'a')
#     threshold = 1.0 if uaid else args.identity_threshold
#     do_cdhit(infile.name, outfile, logfile, threshold, args)
#     if centroid:
#         cent_handle = open(outfile, 'r')
#         # if not uaid:
#         clust_handle = open('{}.clstr'.format(outfile), 'r')
#         sizes = parse_cluster_sizes(clust_handle)
#         seqs = parse_centroids(cent_handle, seq_db, args, sizes=sizes)
#         # else:
#         #     seqs = parse_centroids(cent_handle, seq_db, args)
#     else:
#         clust_handle = open('{}.clstr'.format(outfile), 'r')
#         seqs, sizes = parse_clusters(clust_handle, seq_db, args)
#     if not args.debug:
#         os.unlink(infile.name)
#         os.unlink(os.path.join(args.temp_dir, 'log'))
#         os.unlink(outfile)
#         os.unlink(outfile + '.clstr')
#     if uaid:
#         return seqs
#     return seqs, sizes


def cluster_umis(seq_db, args):
    # sort sequences into UMI bins
    sort_input, input_count = make_umi_sort_input(seq_db, args)
    # logger.info('Found {} sequences'.format(input_count))
    start_time = time.time()
    sorted_file = tempfile.NamedTemporaryFile(dir=args.temp_dir, delete=False)
    sort_cmd = 'sort -k2,2 -o {} {}'.format(sorted_file.name, sort_input)
    logger.info('\nsorting sequences by UMI...')
    p = sp.Popen(sort_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if args.debug:
        logger.info(f'STDOUT\n{stdout}')
        logger.info(f'STDERR\n{stderr}\n')
    logger.info(f'sorting took {round(time.time() - start_time, 2)} seconds')
    umi_bins = process_sorted_umi_file(sorted_file.name)
    # cluster sequences in each UMI bin
    logger.info('\nclustering UMI bins...')
    start_time = datetime.now()
    progress_bar(0, len(umi_bins), start_time=start_time)
    clusters = []
    for i, umi_bin in enumerate(umi_bins, 1):
        seqs = get_clustering_seqs_by_id(umi_bin, seq_db, args)
        bin_clusters = cluster(seqs, threshold=args.identity_threshold)
        if args.only_largest_cluster:
            clusters.append(bin_clusters.largest_cluster)
        else:
            clusters.extend(bin_clusters.clusters)
        progress_bar(i, len(umi_bins), start_time=start_time)
    if not args.debug:
        os.unlink(sort_input)
        os.unlink(sorted_file.name)
    return clusters


def make_umi_sort_input(seq_db, args):
    query = f'SELECT seqs.{args.id_key}, seqs.{args.umi_key} FROM seqs'
    seqs = seq_db.execute(query).fetchall()
    seq_strings = [' '.join(s) for s in seqs]
    sort_input = tempfile.NamedTemporaryFile(dir=args.temp_dir, delete=False)
    sort_input.close()
    with open(sort_input.name, 'w') as f:
        f.write('\n'.join(seq_strings))
    return sort_input.name, len(seq_strings)


def process_sorted_umi_file(sorted_file):
    '''
    Processes a sorted UMI file and returns a list of sequence IDs for each UMI bin.
    '''
    # logger.info('')
    logger.info('processing the sorted UMI file...')
    start = time.time()
    umi_bins = []
    bin_ids = []
    prev = ''
    with open(sorted_file) as f:
        for line in f:
            if not line.strip():
                continue
            seq_id, umi = line.strip().split()
            if umi != prev:
                if bin_ids:
                    umi_bins.append(bin_ids)
                prev = umi
                bin_ids = [seq_id, ]
            else:
                bin_ids.append(seq_id)
        # process the last cluster
        if bin_ids:
            umi_bins.append(bin_ids)
    bin_sizes = [len(u) for u in umi_bins]
    logger.info(f'{len(umi_bins)} UMI bins were identified.\n')
    logger.info(f'the average UMI bin contains {round(np.mean(bin_sizes), 2)} sequences.') 
    logger.info(f'the largest UMI bin contains {max(bin_sizes)} sequences.')
    return umi_bins


def cluster_sequences(seq_db, args):
    # get sequences for clustering
    logger.info('retrieving sequences for clustering...')
    query = f'''SELECT seqs.{args.id_key}, seqs.{args.clustering_key} FROM seqs'''
    result = seq_db.execute(query).fetchall()
    seqs = [Sequence(r[1], id=r[0]) for r in result]
    logger.info(f'Found {len(seqs)} sequences')
    # do the clustering
    logger.info(f'clustering sequences at {args.identity_threshold * 100}% identity...')
    start_time = time.time()
    clusters = cluster(seqs,
                       threshold=args.identity_threshold,
                       temp_dir=args.temp_dir)
    sizes = [c.size for c in clusters]
    logger.info(f'clustering completed in {round(time.time() - start_time, 2)} seconds.')
    logger.info(f'{len(clusters)} clusters were identified.')
    logger.info(f'the average cluster contains {np.mean(sizes)} sequences.')
    return clusters.clusters




# def run_distributed_cdhit(collections, args):
#     seq_dbs = {}
#     fasta_files = {}
#     cluster_files = {}
#     cdhit_output = []
#     logger.info('Building CD-HIT input files...')
#     for i, collection in enumerate(collections):
#         seq_db_name = os.path.basename(collection).replace('.txt', '')
#         update_progress(i, len(collections), sys.stdout, extra_info=seq_db_name)
#         seq_db = get_seqs(None, collection, args, make_seq_db=True, seq_db_name=seq_db_name)
#         seq_dbs[collection] = seq_db
#         fasta_file = make_cdhit_input(seq_db, args)
#         fasta_files[collection] = fasta_file.name
#     update_progress(len(collections), len(collections), sys.stdout, extra_info=' ' * len(seq_db_name))
#     print('\n')
#     node_ips = get_running_celery_node_ips()
#     zip_list = zip(node_ips, cycle(collections)) if len(node_ips) > len(collections) else zip(cycle(node_ips), collections)
#     async_results = []
#     p = mp.Pool()
#     for chunk in chunker(zip_list, size=len(node_ips)):
#         for node_ip, collection in chunk:
#             fasta = fasta_files[collection]
#             clust_filename = seq_db_name = os.path.basename(collection).replace('.txt', '') + '_clust'
#             clust = os.path.join(args.temp_dir, clust_filename)
#             cdhit_output.append(clust)
#             cluster_files[collection] = clust + '.clstr'
#             cdhit_cmd = 'cd-hit -i {} -o {} -c {} -n 5 -d 0 -T 0 -M 35000'.format(fasta, clust, args.identity_threshold)
#             async_results.append(p.apply_async(run_remote, (node_ip, cdhit_cmd)))
#         monitor_mp_jobs(async_results, len(zip_list), trailing_newline=False)
#     print('')
#     for ff in fasta_files.values():
#         os.unlink(ff)
#     for cho in cdhit_output:
#         os.unlink(cho)
#     return cluster_files, seq_dbs


# def run_remote(ip, cmd):
#     with paramiko.SSHClient() as ssh:
#         ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#         ssh.connect(ip)
#         _stdin, stdout, stderr = ssh.exec_command(cmd)
#         while not stdout.channel.exit_status_ready():
#             time.sleep(1)
#     return stdout.read(), stderr.read()





# def make_cdhit_input(seq_db, args, uaid=True):
#     infile = tempfile.NamedTemporaryFile(dir=args.temp_dir, delete=False, mode='w')
#     fastas = []
#     if uaid:
#         seqs = seq_db.execute('''SELECT seqs.seq_id, seqs.uaid FROM seqs''')
#     else:
#         seqs = seq_db.execute('''SELECT seqs.seq_id, seqs.clustering_seq FROM seqs''')
#     for s in seqs:
#         fastas.append('>{}\n{}'.format(s[0], s[1]))
#     infile.write('\n'.join(fastas))
#     infile.close()
#     return infile


# def do_cdhit(fasta, clust, log, threshold, args):
#     seq_type = 'UAIDs' if args.uaid else 'sequences'
#     logger.info('Clustering {} with CD-HIT...'.format(seq_type))
#     start_time = time.time()
#     cdhit_cmd = 'cd-hit -i {} -o {} -c {} -n 5 -d 0 -T 0 -M 35000'.format(fasta, clust, threshold)
#     cluster = sp.Popen(cdhit_cmd, shell=True, stdout=log)
#     cluster.communicate()
#     logger.info('Clustering took {} seconds'.format(round(time.time() - start_time, 2)))
#     logger.info('')


# def parse_centroids(centroid_handle, seq_db, args, sizes=None):
#     cent_ids = []
#     counter = 0
#     for seq in SeqIO.parse(centroid_handle, 'fasta'):
#         # seq_id = '{}_{}'.format(seq.id, sizes[counter]) if sizes is not None else seq.id
#         seq_id = seq.id
#         cent_ids.append(seq_id)
#     return get_clustering_seqs_by_id(cent_ids, seq_db, sizes, args)
    # counter = 0
    # centroids = []
    # for seq in SeqIO.parse(centroid_handle, 'fasta'):
    #     if sizes:
    #         size = sizes[counter]
    #         centroids.append('>{}_{}\n{}'.format(seq.id, size, str(seq.seq)))
    #     else:
    #         centroids.append('>{}\n{}'.format(seq.id, str(seq.seq)))
    #     counter += 1
    # return centroids


# def parse_clusters(cluster_handle, seq_db, args):
#     logger.info('Parsing CD-HIT cluster file...')
#     clusters = [c.split('\n') for c in cluster_handle.read().split('\n>')]
#     logger.info('{} total clusters identified.\n'.format(len(clusters)))
#     logger.info('Retrieving cluster sequences...')
#     cluster_seqs = []
#     start = time.time()
#     lengths = []
#     for cluster in clusters:
#         size = len(cluster) - 1
#         lengths.append(size)
#         if size >= args.min_seqs:
#             cluster_ids = get_cluster_ids(cluster)
#             cluster_seqs.append(get_clustering_seqs_by_id(cluster_ids, seq_db, None, args))
#     logger.info('{} clusters meet the minimum size cutoff ({} sequences)'.format(len(cluster_seqs), args.min_seqs))
#     logger.info('The average cluster contains {} sequences; the largest contains {} sequences.'.format(round(1. * sum(lengths) / len(lengths), 2), max(lengths)))
#     logger.info('Cluster parsing and sequence retrieval took {} seconds\n'.format(round(time.time() - start, 2)))
#     return cluster_seqs, lengths


# def parse_cluster_sizes(cluster_handle):
#     clusters = [c.split('\n') for c in cluster_handle.read().split('\n>')]
#     lengths = []
#     for cluster in clusters:
#         lengths.append(len(cluster) - 1)
#     return lengths


# def get_cluster_ids(cluster):
#     ids = []
#     for c in cluster[1:]:
#         if c:
#             ids.append(c.split()[2][1:-3])
#     return ids


# def get_running_celery_node_ips():
#     celery_info_cmd = '/home/ubuntu/anaconda2/bin/celery -A abstar.utils.queue.celery status'
#     ips = []
#     c = sp.Popen(celery_info_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
#     stdout, stderr = c.communicate()
#     for inst in stdout.split('\n'):
#         if all(['celery@' in inst, 'OK' in inst]):
#             ip = inst.split(':')[0][10:].replace('-', '.')
#             ips.append(ip)
#     return ips





# =========================================
#
#       UAID CENTROID CLUSTERING
#
# =========================================



# def get_uaid_centroids(uaid_clusters, args):
#     logger.info('Calculating centroid sequences with USEARCH:')
#     start_time = time.time()
#     centroids = []
#     singletons = [c[0] for c in uaid_clusters if len(c) == 1]
#     for s in singletons:
#         seq_id = s.split('\n')[0].replace('>', '')
#         seq = s.split('\n')[1]
#         centroids.append('>{}\n{}'.format(seq_id, seq))
#     sizes = [1] * len(centroids)
#     clusters = [c for c in uaid_clusters if len(c) > 1]
#     if args.debug:
#         for cluster in clusters:
#             centroid, size = do_usearch_centroid(cluster, vars(args))
#             centroids.extend(centroid)
#             sizes.extend(size)
#     elif args.cluster:
#         async_results = [do_usearch_centroid.delay(clusters[cs:cs + args.chunksize], vars(args)) for cs in range(0, len(clusters), args.chunksize)]
#         succeeded, failed = monitor_celery_jobs(async_results, len(clusters), args.chunksize)
#         results = []
#         for ar in async_results:
#             results.extend(ar.get())
#     else:
#         p = mp.Pool(maxtasksperchild=100)
#         async_results = []
#         for c in clusters:
#             async_results.append(p.apply_async(do_usearch_centroid, (c, vars(args))))
#         monitor_mp_jobs(async_results, len(clusters), args.chunksize)
#         for a in async_results:
#             centroid, size = a.get()
#             centroids.extend(centroid)
#             sizes.append(size)
#         p.close()
#         p.join()
#         logger.info('Centroids were calculated in {} seconds.'.format(round(time.time() - start_time), 2))
#     return centroids, sizes


# @celery.task
# def do_usearch_centroid(uaid_groups, arg_dict):
#     # '''
#     # Clusters sequences at 90% identity using USEARCH.

#     # Inputs
#     # uaid_group_seqs: a list of fasta strings corresponding to sequences from a single UAID group.

#     # Outputs
#     # A list of fasta strings, containing centroid sequences for each cluster.
#     # '''
#     args = Args(**arg_dict)
#     all_centroid_seqs = []
#     all_sizes = []
#     for uaid_group_seqs in uaid_groups:
#         fasta = tempfile.NamedTemporaryFile(dir=args.temp_dir, prefix='cluster_input_', delete=False, mode='w')
#         results = tempfile.NamedTemporaryFile(dir=args.temp_dir, prefix='results_', delete=False, mode='w')
#         centroids = tempfile.NamedTemporaryFile(dir=args.temp_dir, prefix='centroids_', delete=False, mode='w')
#         fasta.write('\n'.join(uaid_group_seqs))
#         fasta.close()
#         usearch = ['usearch',
#                    '-cluster_fast',
#                    fasta.name,
#                    '-maxaccepts', '0',
#                    '-maxrejects', '0',
#                    '-id', '0.9',
#                    '-sizeout',
#                    '-uc', results.name,
#                    '-centroids', centroids.name]
#         p = sp.Popen(usearch, stdout=sp.PIPE, stderr=sp.PIPE)
#         stdout, stderr = p.communicate()
#         centroid_seqs = []
#         sizes = []
#         for cent in SeqIO.parse(open(centroids.name, 'r'), 'fasta'):
#             seq_id = cent.id.split(';')[0]
#             size = int(cent.id.split(';')[1].replace('size=', ''))
#             centroid_seqs.append('>{}\n{}'.format(seq_id, str(cent.seq)))
#             sizes.append(size)
#         if args.only_largest_cluster:
#             cents_plus_sizes = sorted(zip(centroid_seqs, sizes), key=lambda x: x[1], reverse=True)
#             centroid_seqs = [cents_plus_sizes[0][0], ]
#             sizes = [cents_plus_sizes[0][1], ]
#         os.unlink(fasta.name)
#         os.unlink(results.name)
#         os.unlink(centroids.name)
#         all_centroid_seqs += centroid_seqs
#         all_sizes += sizes
#     return all_centroid_seqs, all_sizes





# =========================================
#
#         CONSENTROID SEQUENCES
#
# =========================================


def calculate_consentroids(clusters, seq_db, args):
    seq_type = 'consensus' if args.consensus else 'centroid'
    logger.info('')
    logger.info(f'computing {seq_type} sequences...')
    start_time = datetime.now()
    extra_info = 'passed cluster size threshold: {} / {}'
    progress_bar(0, len(clusters), start_time=start_time, extra_info=extra_info.format(0, 0))
    sequences = []
    passed = 0
    for i, c in enumerate(clusters):
        if c.size < args.min_seqs:
            progress_bar(i, len(clusters), start_time=start_time, extra_info=extra_info.format(passed, i))
            continue
        elif c.size == 1:
            output_seqs = get_output_seqs_by_id(c.seq_ids, seq_db, args)
            consentroid = output_seqs[0]
            output_cluster_size = 1
        else:
            output_seqs = get_output_seqs_by_id(c.seq_ids, seq_db, args)
            output_clusters = cluster(output_seqs,
                                    threshold=args.identity_threshold * 0.8,
                                    temp_dir=args.temp_dir)
            output_cluster = output_clusters.largest_cluster
            consentroid = output_cluster.consensus if args.consensus else output_cluster.centroid
            output_cluster_size = output_cluster.size
        if args.include_cluster_size:
            consentroid.id = f'{consentroid.id}_{output_cluster_size}'
        sequences.append(consentroid)
        passed += 1
        progress_bar(i, len(clusters), start_time=start_time, extra_info=extra_info.format(passed, i))
    logger.info('')
    return sequences














# def get_consensus(clusters, germs, args):
#     fastas = []
#     all_sizes = []
#     if args.min_seqs == 1:
#         singletons = [c[0] for c in clusters if len(c) == 1]
#         clusters = [c for c in clusters if len(c) > 1]
#         singleton_count = len(singletons)
#         logger.info('')
#         logger.info('{} clusters contained only a single sequence. Processing singletons...'.format(len(singletons)))
#         progress_bar(0, singleton_count)
#         for i, s in enumerate(singletons):
#             seq_id = str(uuid.uuid4())
#             if args.include_cluster_size:
#                 seq_id += '_1'
#             seq = s.split('\n')[-1]
#             fastas.append('>{}\n{}'.format(seq_id, seq))
#             all_sizes.append(1)
#             progress_bar(i + 1, singleton_count)
#         logger.info('')
#     logger.info('')
#     logger.info('Building consensus sequences...')
#     if args.debug:
#         results = []
#         cluster_count = len(clusters)
#         progress_bar(0, cluster_count)
#         for i, cluster in enumerate(clusters):
#             results.append(do_usearch_consensus(cluster, germs, vars(args)))
#             progress_bar(i + 1, cluster_count)
#         logger.info('')
#         logger.info('')
#         # results = [do_usearch_consensus(cluster, germs, args) for cluster in clusters]
#     elif args.cluster:
#         async_results = [do_usearch_consensus.delay(clusters[cs:cs + args.chunksize], germs, vars(args)) for cs in range(0, len(clusters), args.chunksize)]
#         succeeded, failed = monitor_celery_jobs(async_results, len(clusters), args.chunksize)
#         results = []
#         for ar in async_results:
#             results.append(ar.get())
#     else:
#         p = mp.Pool(maxtasksperchild=50)
#         # async_results = [p.apply_async(calculate_consensus, (cluster, germs, args)) for cluster in clusters]
#         async_results = [p.apply_async(do_usearch_consensus, (clusters[cs:cs + args.chunksize], germs, vars(args))) for cs in range(0, len(clusters), args.chunksize)]
#         monitor_mp_jobs(async_results, len(clusters), args.chunksize)
#         results = []
#         for ar in async_results:
#             results.append(ar.get())
#         p.close()
#         p.join()
#     for r in results:
#         seqs = r[0]
#         sizes = r[1]
#         all_sizes.extend(sizes)
#         fastas.extend(seqs)
#     return fastas, all_sizes


# @celery.task
# def do_usearch_consensus(clusters, germs, arg_dict):
#     # '''
#     # Clusters sequences at using USEARCH and calculates consensus sequences for each cluster.

#     # Inputs
#     # cluster, germs: a list of fasta strings corresponding to sequences from a single UAID group.

#     # Outputs
#     # A list of fasta strings, containing consensus sequences for each cluster.
#     # '''
#     args = Args(**arg_dict)
#     all_consensus_seqs = []
#     all_sizes = []
#     for cluster in clusters:
#         fasta = tempfile.NamedTemporaryFile(dir=args.temp_dir, prefix='cluster_input_', delete=False, mode='w')
#         results = tempfile.NamedTemporaryFile(dir=args.temp_dir, prefix='results_', delete=False, mode='w')
#         consensus = tempfile.NamedTemporaryFile(dir=args.temp_dir, prefix='consensus_', delete=False, mode='w')
#         fasta.write('\n'.join(cluster))
#         fasta.close()
#         usearch = ['usearch',
#                    '-cluster_fast',
#                    fasta.name,
#                    '-maxaccepts', '0',
#                    '-maxrejects', '0',
#                    '-id', str(args.identity_threshold),
#                    '-sizeout',
#                    '-uc', results.name,
#                    '-consout', consensus.name,
#                    ]
#         p = sp.Popen(usearch, stdout=sp.PIPE, stderr=sp.PIPE)
#         stdout, stderr = p.communicate()
#         consensus_seqs = []
#         sizes = []
#         with open(consensus.name, 'r') as f:
#             for cons in SeqIO.parse(f, 'fasta'):
#                 seq_id = str(uuid.uuid4())
#                 size = int(cons.id.split(';')[1].replace('size=', ''))
#                 if args.include_cluster_size:
#                     seq_id += '_{}'.format(size)
#                 consensus_seqs.append('>{}\n{}'.format(seq_id, str(cons.seq)))
#                 sizes.append(size)
#             if args.only_largest_cluster:
#                 cons_plus_sizes = sorted(zip(consensus_seqs, sizes), key=lambda x: x[1], reverse=True)
#                 consensus_seqs = [cons_plus_sizes[0][0], ]
#                 sizes = [cons_plus_sizes[0][1], ]
#         if not args.debug:
#             os.unlink(fasta.name)
#             os.unlink(results.name)
#             os.unlink(consensus.name)
#         all_consensus_seqs += consensus_seqs
#         all_sizes += sizes
#     return all_consensus_seqs, all_sizes


# def calculate_consensus(cluster, germs, args):
#     if len(cluster) == 1:
#         return (cluster[0].split('\n')[1].upper(), 1)
#     fasta_string = consensus_alignment_input(cluster, germs, args)

#     if len(cluster) < 100:
#         alignment = muscle(fasta_string)
#     elif len(cluster) < 1000:
#         alignment = muscle(fasta_string, maxiters=2)
#     else:
#         alignment = muscle(fasta_string, maxiters=1, diags=True)
#     ambig = 'X' if args.aa else 'N'
#     summary_align = AlignInfo.SummaryInfo(alignment)
#     consensus = summary_align.gap_consensus(threshold=0.51, ambiguous=ambig)
#     consensus = str(consensus).replace('-', '')
#     return (consensus.upper(), len(cluster))


# def consensus_alignment_input(cluster, germs, args):
#     if germs is not None:
#         try:
#             v_gene = vgene_lookup(cluster, args)
#             germ = '>{}\n{}'.format(v_gene, germs[v_gene])
#         except KeyError:
#             germ = ''
#     else:
#         germ = ''
#     cluster.append(germ)
#     return '\n'.join(cluster)


# def vgene_lookup(cluster, args):
#     v_genes = []
#     seq_ids = [seq.split('\n')[0].replace('>', '') for seq in cluster]
#     db_path = os.path.join(args.temp_dir, 'seq_db')
#     conn = sqlite3.connect(db_path)
#     seq_db = conn.cursor()
#     for chunk in chunker(seq_ids):
#         v_chunk = seq_db.execute('''SELECT seqs.v_gene
#                                     FROM seqs
#                                     WHERE seqs.seq_id IN ({})'''.format(','.join('?' * len(chunk))), chunk)
#         v_genes.extend(v_chunk)
#     v_gene = Counter(v_genes).most_common()[0][0]
#     conn.close()
#     return v_gene





# =========================================
#
#         PROGRESS AND PRINTING
#
# =========================================



def print_start_info(args):
    logger.info('')
    logger.info('')
    logger.info('')
    logger.info('-' * 25)
    logger.info('ERROR CORRECTION')
    logger.info('-' * 25)
    log_params(args)


def log_params(args):
    logger.info('')
    logger.info('INPUT DB: {}'.format(args.db))
    logger.info('OUTPUT DIR: {}'.format(args.output))
    logger.info('TEMP DIR: {}'.format(args.temp_dir))
    logger.info('UMIs: {}'.format(args.umi))
    logger.info('CLUSTERING TYPE: {}'.format('sort/uniq' if args.non_redundant else 'CD-HIT'))
    if args.non_redundant:
        logger.info('IDENTITY THRESHOLD: 1.0')
        logger.info('MIN SEQS: 1')
        return
    elif args.umi == 0:
        logger.info('IDENTITY THRESHOLD: {}'.format(args.identity_threshold))
        logger.info('GERMLINES: {}'.format(args.germs))
    else:
        logger.info('PARSE UMIs: {}'.format(args.parse_umis))
    logger.info('MIN SEQS: {}'.format(args.min_seqs))
    logger.info('LARGEST CLUSTER ONLY: {}'.format(args.only_largest_cluster))



# def monitor_mp_jobs(results, total, chunksize=1, trailing_newline=True):
#     finished = 0
#     # jobs = len(results)
#     jobs = total
#     while finished < jobs:
#         time.sleep(1)
#         ready = [ar for ar in results if ar.ready()]
#         finished = min(len(ready) * chunksize, jobs)
#         update_progress(finished, jobs, sys.stdout)
#     if trailing_newline:
#         sys.stdout.write('\n')


# def monitor_celery_jobs(results, total, chunksize):
#     finished = 0
#     # jobs = len(results)
#     jobs = total
#     while finished < jobs:
#         time.sleep(1)
#         succeeded = [ar for ar in results if ar.successful()]
#         failed = [ar for ar in results if ar.failed()]
#         finished = min((len(succeeded) + len(failed)) * chunksize, jobs)
#         update_progress(finished, jobs, sys.stdout)
#     sys.stdout.write('\n\n')
#     return succeeded, failed


# def update_progress(finished, jobs, log, failed=None, extra_info=None):
#     pct = int(100. * finished / jobs)
#     ticks = int(pct / 2)
#     spaces = int(50 - ticks)
#     if failed is not None:
#         prog_bar = '\r({}/{}) |{}{}|  {}% ({}, {})'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct, finished - failed, failed)
#     elif extra_info is not None:
#         prog_bar = '\r({}/{}) |{}{}|  {}% ({})'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct, extra_info)
#     else:
#         prog_bar = '\r({}/{}) |{}{}|  {}%'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct)
#     sys.stdout.write(prog_bar)


def write_output(group, sequences, sizes, group_start_time, args):
    seq_type = 'consensus' if args.consensus else 'centroid'
    logger.info('writing {} sequences to output file...'.format(seq_type))
    write_fasta_output(group, sequences, args)
    if sizes is not None:
        write_stats_output(group, sizes, args)
    logger.info('{} {} sequences were identified.'.format(len(sequences), seq_type))
    logger.info('{} was processed in {} seconds.'.format(group,
        round(time.time() - group_start_time, 2)))
    logger.info('')


def write_nr_output(group, unique_file, group_start_time, args):
    logger.info('writing non-redundant sequences to output file...')
    oname = group
    if os.path.isfile(group):
        oname = '.'.join(os.path.basename(group).split('.')[:-1])
    outfile = os.path.join(args.output, '{}.fasta'.format(oname))
    open(outfile, 'w').write('')
    ohandle = open(outfile, 'a')
    count = 0
    with open(unique_file) as f:
        for line in f:
            if len(line.strip().split()) != 3:
                continue
            seq_id, vdj, raw = line.strip().split()
            ohandle.write('>{}\n{}\n'.format(seq_id, raw))
            count += 1
    logger.info('{} non-redundant sequences were identified.'.format(count))
    logger.info('{} was processed in {} seconds.'.format(group,
        round(time.time() - group_start_time, 2)))
    logger.info('')


def write_fasta_output(group, sequences, args):
    seq_type = 'consensus' if args.consensus else 'centroids'
    oname = group
    if os.path.isfile(group):
        oname = '.'.join(os.path.basename(group).split('.')[:-1])
    out_dir = os.path.join(args.output, seq_type)
    make_dir(out_dir)
    outfile = os.path.join(out_dir, '{}.fasta'.format(oname))
    with open(outfile, 'w') as f:
        f.write('\n'.join([s.fasta for s in sequences]))


def write_stats_output(group, sizes, args):
    sizes = [int(s) for s in sizes]
    bin_counts = np.bincount(sizes)[1:]
    bins = list(range(1, len(bin_counts)))
    bin_string = 'Cluster Size,Count\n'
    bin_string += '\n'.join([f'{size},{count}' for size, count in zip(bins, bin_counts)])
    oname = group
    if os.path.isfile(group):
        oname = '.'.join(os.path.basename(group).split('.')[:-1])
    out_dir = os.path.join(args.output, 'cluster_sizes')
    make_dir(out_dir)
    outfile = os.path.join(out_dir, '{}.csv'.format(oname))
    with open(outfile, 'w') as f:
        f.write(bin_string)


def print_group_info(group):
    logger.info('')
    logger.info('')
    logger.info(group)
    logger.info('-' * len(group))


# def countdown(args):
#     start = time.time()
#     h = int(args.sleep / 60)
#     m = int(args.sleep % 60)
#     hz = '0' if len(str(h)) == 1 else ''
#     mz = '0' if len(str(m)) == 1 else ''
#     print('Countdown: {}{}:{}{}:{}'.format(hz, h, mz, m, '00'), end='')
#     sys.stdout.flush()
#     done = False
#     while not done:
#         time.sleep(0.25)
#         elapsed = int(time.time() - start)
#         if elapsed < args.sleep * 60:
#             h = int((args.sleep - elapsed / 60.) / 60)
#             m = int((args.sleep - elapsed / 60.) % 60)
#             s = int((60 * args.sleep - elapsed) % 60)
#             hz = '0' if len(str(h)) == 1 else ''
#             mz = '0' if len(str(m)) == 1 else ''
#             sz = '0' if len(str(s)) == 1 else ''
#             print('\rCountdown: {}{}:{}{}:{}{}'.format(hz, h, mz, m, sz, s), end='')
#             sys.stdout.flush()
#         else:
#             print('\rCountdown: 00:00:00')
#             done = True


def run(**kwargs):
    '''
    Corrects antibody reads using UMIs (unique molecular identifiers) or identity-based clustering.

    Either ``json``, ``tabular`` or ``db`` is required.


    Args:

        db (str): Name of a MongoDB database to query.

        collection (str): Name of a MongoDB collection. If not provided, all
            collections in ``db`` will be iteratively processed.

        json (str): Can be one of two things:

            1. Path to a JSON file, containing sequence data annotated by AbStar.

            2. Path to a directory, containing one or more JSON files of
                AbStar-annotated data.

        tabular (str): Can be one of two things:

            1. Path to a tabular file, containing annotated sequence data.

            2. Path to a directory, containing one or more tabular files of
                AbStar-annotated data.

        output (str): Path to the output directory, into which corrected FASTA
            files will be deposited. If it does not exist, it will be created.

        log (str): Path to the log file. If the parent directory doesn't exist,
            it will be created.

        ip (str): IP address of the MongoDB server. Default is ``localhost``.

        port (int): Port of the MongoDB server. Default is ``27017``.

        sep (str): Character used to separate fields in tabular input files. Default is ``'\t'``.

        user (str): Username with which to connect to the MongoDB database. If either
            of ``user`` or ``password`` is not provided, the connection to the MongoDB
            database will be attempted without authentication.

        password (str): Password with which to connect to the MongoDB database. If either
            of ``user`` or ``password`` is not provided, the connection to the MongoDB
            database will be attempted without authentication.

        min_seqs (int): Minimum number of sequences for a centroid/consensus sequence to be
            calculated. After clustering (either by identity or using UMIs), clusters with
            at least ``min_seqs`` sequences will be retained for consensus/centroid calculation.
            Default is ``1``.

        umi (bool): If ``True``, use unique molecular identifiers (UMIs) for error correction.
            Sequences will be binned by UMI and the sequences in each bin will be used to
            compute a centroid or consensus sequence. If ``False``, sequences will be clustered
            by identity and each cluster will be used for consensus/centroid determination.

        parse_umis (int): If UMIs haven't been pre-parsed, indicate the length of the
            UMI sequence (in nucleotides) and the UMIs will be parsed from the 5' end of each 
            sequence prior to clustering and correction. If ``parse_umis`` is negative, the 
            UMI will be parsed from the 3' end of the sequence. Default is ``0``, which does not 
            parse a UMI sequence.

        consensus (bool): If ``True``, consensus sequences are calculated. If ``False``, centroid
            sequences are calculated. Default is ``True``.

        identity_threshold (float): Identity threshold, if clustering by identity (not UMIs).
            Must be a floating point number between 0 and 1. Default is 0.975.

        only_largest_cluster (bool): When clustering using UMIs, there is a some probability that
            different sequences get labeled with the same UMI or that PCR crossover can produce hybrid sequences
            within the same UMI bin that may confound consensus computation. To limit incorrect consensus/centroid
            calculation, sequences in each UMI bin are clustered using ``identity_threshold`` before
            calculating consens/centroid sequences. By default, all UMI clusters that meet the
            ``min_seqs`` size threshold are used to generate consensus/centroid sequences. If that
            behavior is not desired, setting ``only_largest_cluster`` to ``True`` results in only
            the largest cluster in each UMI bin being used to generate centroid/consensus sequences.

        nr (bool): If ``True``, a non-redundant sequence dataset will be generated using ``sort | uniq``.
            This is much faster than normal sequence clustering with mmseqs, but can only be performed at an
            identity threshold of 100%.

            .. note:

                Using ``nr`` may produce different results than clustering sequences with ``identity_threshold``
                set to ``1.0``. This is because sequences of different lengths that are otherwise identical
                will not be collapsed when using ``nr`` but will be collapsed using normal clustering.

        debug (bool): If ``True``, logging is more verbose.
    '''
    args = Args(**kwargs)
    global logger
    logger = log.get_logger('abcorrect')
    main(args)


def run_standalone(args):
    logfile = args.log if args.log else os.path.join(args.output, 'abcorrect.log')
    log.setup_logging(logfile)
    global logger
    logger = log.get_logger('abcorrect')
    main(args)


def main(args):
    print_start_info(args)
    for d in [args.output, args.temp_dir]:
        make_dir(d)
    # check whether JSON files have been passed
    if args.json is not None:
        if os.path.isfile(args.json):
            groups = [args.json, ]
        elif os.path.isdir(args.json):
            groups = list_files(args.json)
        else:
            err = f'ERROR: The supplied JSON path ({args.json}) is neither a file nor a directory. '
            err += 'Please double-check your JSON path.'
            print('\n')
            print(err)
            print('\n')
            sys.exit()
    # check whether TABULAR files have been passed
    elif args.tabular is not None:
        if os.path.isfile(args.tabular):
            groups = [args.tabular, ]
        elif os.path.isdir(args.tabular):
            groups = list_files(args.tabular)
        else:
            err = f'ERROR: The supplied TABULAR path ({args.tabular}) is neither a file nor a directory. '
            err += 'Please double-check your TABULAR path.'
            print('\n')
            print(err)
            print('\n')
            sys.exit()
    # otherwise, get sequences from MongoDB
    else:
        if args.collection is not None:
            groups = [args.collection, ]
        else:
            db = mongodb.get_db(args.db, args.ip, args.port, args.user, args.password)
            groups = mongodb.get_collections(db)

    for group in groups:
        start_time = time.time()
        print_group_info(group)
        if args.non_redundant:
            seqs = get_seqs(group, args)
            unique_file = sort_unique(seqs, args)
            write_nr_output(group, unique_file, start_time, args)
        else:
            seqs = get_seqs(group, args)
            seq_db = build_seq_db(seqs, args)
            if args.umi:
                clusters = cluster_umis(seq_db, args)
            else:
                clusters = cluster_sequences(seq_db, args)
            sizes = [c.size for c in clusters]
            consentroids = calculate_consentroids(clusters, seq_db, args)
            write_output(group, consentroids, sizes, start_time, args)
            remove_sqlite_db(args)


if __name__ == '__main__':
    parser = parse_args()
    args = parser.parse_args()
    logfile = args.log if args.log else os.path.join(args.output, 'abcorrect.log')
    log.setup_logging(logfile)
    logger = log.get_logger('abcorrect')
    main(args)
