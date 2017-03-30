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



from __future__ import print_function

from collections import Counter
import itertools
import json
import math
import multiprocessing as mp
import os
import shelve
import sqlite3
from StringIO import StringIO
import subprocess as sp
import sys
import tempfile
import time
import urllib
import uuid

import numpy as np

from pymongo import MongoClient

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo

from abtools import log, mongodb
from abtools.alignment import mafft, muscle
from abtools.cluster import cdhit, cluster
from abtools.pipeline import make_dir, list_files
from abtools.sequence import Sequence


def parse_args():
    import argparse
    parser = argparse.ArgumentParser("Clusters sequences using either an identity threshold or unique \
                                      antibody identifiers (UAIDs). \
                                      Calculates either the centroid or consensus sequence for each \
                                      identity/UAID cluster passing a cluster size threshold.")
    parser.add_argument('-d', '--database', dest='db', default=None,
                        help="Name of the MongoDB database to query. Required")
    parser.add_argument('-c', '--collection', dest='collection', default=None,
                        help="Name of the MongoDB collection to query. \
                        If not provided, all collections in the given database will be processed iteratively.")
    parser.add_argument('-j', '--json', dest='json', default=None,
                        help="Input JSON file or directory of JSON files.")
    parser.add_argument('-o', '--output', dest='output', required=True,
                        help="Output directory for the FASTA files. Required")
    parser.add_argument('-l', '--log', dest='log',
                        help="Location for the log file. \
                        If not provided, log will be written to <output>/abcorrect.log")
    parser.add_argument('-t', '--temp', dest='temp_dir', required=True,
                        help="Directory for temporary files, including the SQLite database. Required")
    parser.add_argument('-i', '--ip', dest='ip', default='localhost',
                        help="IP address for the MongoDB server.  Defaults to 'localhost'.")
    parser.add_argument('--port', dest='port', default=27017, type=int,
                        help="IP address for the MongoDB server.  Defaults to 'localhost'.")
    parser.add_argument('-u', '--user', dest='user', default=None,
                        help="Username for the MongoDB server. Not used if not provided.")
    parser.add_argument('-p', '--password', dest='password', default=None,
                        help="Password for the MongoDB server. Not used if not provided.")
    parser.add_argument('-m', '--min', dest='min_seqs', default=1, type=int,
                        help="Minimum number of sequences for finding centroids from a UAID group. \
                        Defaults to 1.")
    parser.add_argument('-U', '--no-uaid', dest='uaid', action='store_false', default=True,
                        help="Clusters sequences by identity rather than using universal antibody IDs (UAIDs).")
    parser.add_argument('-P', '--parse-uaids', dest='parse_uaids', type=int, default=0,
                        help="Length of the UAID to parse, if the UAID was not parsed during AbStar processing. \
                        If the '--no-uaid' flag is also used, this option will be ignored. \
                        For a UAID of length 20, option should be passed as '--parse-uaids 20'. \
                        Default is to not parse UAIDs.")
    parser.add_argument('-C', '--centroid', dest='consensus', action='store_false', default=True,
                        help="Generates consensus sequences for each UAID or homology cluster. \
                        Default is to identify cluster centroids.")
    parser.add_argument('-I', '--identity_threshold', dest='identity_threshold', default=0.975, type=float,
                        help="If performing identity-based clustering, this threshold is used for the clustering. \
                        If performing UID-based correction, this threshold is used when clustering sequences assigned \
                        to the same UID bin. Default is 0.975.")
    parser.add_argument('--only-largest-cluster', default=False, action='store_true',
                        help="If set while calculating centroids using UAIDs, \
                        only the largest centroid for each UAID cluster is retained.")
    parser.add_argument('--non-redundant', dest='non_redundant', action='store_true', default=False,
                        help='If set, will make a non-redundant FASTA file using Unix sort. \
                        This is much faster than CD-HIT, but can only remove exact duplicates. \
                        Note that this is not exactly equivalent to clustering with a threshold of 1.0, \
                        because sequences of different length that align perfectly for the entirety \
                        of the shorter sequence will be collapsed when using CD-HIT but will not be \
                        collapsed by Unix sort.')
    parser.add_argument('-g', '--germs', dest='germs', default=None,
                        help="Path to a FASTA-formatted file of germline V gene sequences. \
                        Required if building consensus sequences, not required for centroids.")
    parser.add_argument('--aa', dest='aa', action='store_true', default=False,
                        help="If set, use amino acid sequences for identity-based clustering.")
    parser.add_argument('--clustering_field', dest='clustering_field', default='vdj_nt',
                        help="The MongoDB field to be used for clustering. \
                        Default is 'vdj_nt'.")
    parser.add_argument('--output_field', dest='output_field', default='oriented_input',
                        help="The MongoDB field to be used for creating consensus/centroid sequences. \
                        Default is 'oriented_input'.")
    parser.add_argument('-D', '--debug', dest='debug', action='store_true', default=False,
                        help="If set, will run in debug mode.")
    parser.add_argument('-s', '--sleep', dest='sleep', type=int, default=0,
                        help="Delay, in minutes, until the script starts executing. Default is 0.")
    return parser


class Args(object):
    def __init__(self, db=None, collection=None, json=None,
                 output=None, log=None, temp=None,
                 ip='localhost', port=27017, user=None, password=None,
                 min_seqs=1, identity_threshold=0.975,
                 uaid=True, parse_uaids=0, non_redundant=False,
                 consensus=True, only_largest_cluster=False, germs=None,
                 aa=False, field='vdj_nt', debug=False, sleep=0):
        super(Args, self).__init__()
        if not all([db, output, temp]):
            print('\nERROR: Output and temp directories must be provided, \
                as well as a MongoDB database name.\n')
            sys.exit(1)
        self.db = db
        self.collection = collection
        self.json = json
        self.output = output
        self.log = log
        self.temp_dir = temp
        self.ip = ip
        self.port = int(port)
        self.user = user
        self.password = password
        self.min_seqs = int(min_seqs)
        self.uaid = False if non_redundant else uaid
        self.parse_uaids = int(parse_uaids)
        self.consensus = consensus
        self.identity_threshold = float(identity_threshold)
        self.only_largest_cluster = only_largest_cluster
        self.non_redundant = non_redundant
        self.germs = germs
        self.aa = aa
        self.field = field
        self.debug = debug
        self.sleep = int(sleep)



# =========================================
#
#        SEQUENCES AND DATABASES
#
# =========================================



def get_seqs(db, collection, args, make_seq_db=True):
    seqs = query(db, collection, args)
    if not make_seq_db:
        return seqs
    return build_seq_db(seqs, args)


def query(db, collection, args):
    # parse JSON file...
    if db is None:
        logger.info('Reading JSON file...')
        results = []
        _results = []
        with open(collection) as f:
            for line in f:
                if line.strip():
                    _results.append(json.loads(line.strip()))
        for r in _results:
            raw_field = 'raw_query' if 'raw_query' in r else 'raw_input'
            try:
                d = {'seq_id': r['seq_id'],
                     args.clustering_field: r[args.clustering_field],
                     args.output_field: r[args.output_field],
                     'raw_query': r[raw_field],
                     'v_gene': {'full': r['v_gene']['full']}}
                if args.uaid and not args.non_redundant:
                    uid_field = 'uid' if 'uid' in r else 'uaid'
                    if uid_field in r:
                        d['uaid'] = r[uid_field]
                    elif args.parse_uaids:
                        if args.parse_uaids > 0:
                            d['uaid'] = r[raw_field][:args.parse_uaids]
                        else:
                            d['uaid'] = r[raw_field][args.parse_uaids:]
                    else:
                        err = 'ERROR: UAID field was not found. '
                        err += 'Ensure that UAIDs were parsed by AbStar, '
                        err += 'use the -parse_uaids option to parse the UAIDs from the raw query input, '
                        err += 'or use the -u option for identity-based clustering.'
                        raise ValueError(err)
                results.append(d)
            except KeyError:
                continue
    # ...or get sequences from a MongoDB database
    else:
        logger.info('Getting sequences from MongoDB...')
        coll = db[collection]
        match = {'prod': 'yes'}
        project = {'_id': 0, 'seq_id': 1, args.clustering_field: 1, args.output_field: 1, 'raw_query': 1, 'raw_input': 1, 'v_gene.full': 1}
        if args.uaid is not None:
            project['uaid'] = 1
            project['uid'] = 1
        results = coll.find(match, project)
    if args.non_redundant:
        return results
    elif args.uaid:
        seqs = []
        for r in results:
            # raw_field and uid_field are necessary to maintain compatibility with legacy AbStar versions
            raw_field = 'raw_query' if 'raw_query' in r else 'raw_input'
            uid_field = 'uid' if 'uid' in r else 'uaid'
            if 'uaid' in r:
                seqs.append((r['seq_id'], r[uid_field], r[args.clustering_field], r[args.output_field], r[raw_field], r['v_gene']['full']))
            elif args.parse_uaids:
                if args.parse_uaids > 0:
                    seqs.append((r['seq_id'], r[raw_field][:args.parse_uaids], r[args.clustering_field], r[args.output_field], r[raw_field], r['v_gene']['full']))
                else:
                    seqs.append((r['seq_id'], r[raw_field][args.parse_uaids:], r[args.clustering_field], r[args.output_field], r[raw_field], r['v_gene']['full']))
            else:
                err = 'ERROR: UAID field was not found. '
                err += 'Ensure that UAIDs were parsed by AbStar, '
                err += 'use the -parse_uaids option to parse the UAIDs from the raw query input, '
                err += 'or use the -u option for identity-based clustering.'
                raise ValueError(err)
    else:
        seqs = [(r['seq_id'], r[seq_field], r[args.clustering_field], r[args.output_field], r[raw_field], r['v_gene']['full']) for r in results]
    logger.info('Found {} sequences\n'.format(len(seqs)))
    return seqs


def build_seq_db(seqs, args):
    logger.info('Building a SQLite database of sequences...')
    db_path = os.path.join(args.temp_dir, 'seq_db')
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    create_cmd = get_seq_db_creation_cmd(args)
    insert_cmd = get_seq_db_insert_cmd(args)
    c.execute('DROP TABLE IF EXISTS seqs')
    c.execute(create_cmd)
    c.executemany(insert_cmd, seqs)
    c.execute('CREATE INDEX seq_index ON seqs (seq_id)')
    conn.commit()
    conn.close()
    return db_path


def get_seq_db_creation_cmd(args):
    if args.uaid:
        return '''CREATE TABLE seqs (seq_id text, uaid text, clustering_seq text, output_seq text, raw text, v_gene text)'''
    return '''CREATE TABLE seqs (seq_id text, clustering_seq text, output_seq text, raw text, v_gene text)'''


def get_seq_db_insert_cmd(args):
    if args.uaid:
        return 'INSERT INTO seqs VALUES (?,?,?,?,?,?)'
    return 'INSERT INTO seqs VALUES (?,?,?,?,?)'


def remove_sqlite_db(args):
    db_path = os.path.join(args.temp_dir, 'seq_db')
    os.unlink(db_path)


def parse_germs(germ_file):
    germ_handle = open(germ_file, 'r')
    germs = {}
    for seq in SeqIO.parse(germ_handle, 'fasta'):
        germs[seq.id] = str(seq.seq).upper()
    return germs


def chunker(l, size=900):
    return (l[pos:pos + size] for pos in xrange(0, len(l), size))


def retrieve_clustering_seqs(seq_ids, seq_db_path):
    conn = sqlite3.connect(seq_db_path)
    seq_db = conn.cursor()
    seqs = []
    for chunk in chunker(seq_ids):
        seq_chunk = seq_db.execute('''SELECT seqs.seq_id, seqs.clustering_seq
                                   FROM seqs
                                   WHERE seqs.seq_id IN ({})'''.format(','.join('?' * len(chunk))), chunk)
        seqs.extend(seq_chunk)
    conn.commit()
    conn.close()
    return [Sequence(s[1], id=s[0]) for s in seqs]


def retrieve_output_seqs(seq_ids, seq_db_path):
    conn = sqlite3.connect(seq_db_path)
    seq_db = conn.cursor()
    seqs = []
    for chunk in chunker(seq_ids):
        seq_chunk = seq_db.execute('''SELECT seqs.seq_id, seqs.output_seq
                                   FROM seqs
                                   WHERE seqs.seq_id IN ({})'''.format(','.join('?' * len(chunk))), chunk)
        seqs.extend(seq_chunk)
    conn.commit()
    conn.close()
    return [Sequence(s[1], id=s[0]) for s in seqs]


# =========================================
#
#             NON-REDUNDANT
#
# =========================================


def unix_sort_unique(seqs, args):
    sort_file, input_count = make_sort_input(seqs, args)
    logger.info('Found {} sequences'.format(input_count))
    unique_file = tempfile.NamedTemporaryFile(dir=args.temp_dir, delete=False)
    sort_unique_cmd = 'sort -k2,2 -u -o {} {}'.format(unique_file.name, sort_file)
    logger.info('Running sort/uniq...')
    p = sp.Popen(sort_unique_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    p.communicate()
    os.unlink(sort_file)
    return unique_file.name



def make_sort_input(seqs, args):
    raw_field = 'raw_query' if 'raw_query' in seqs[0] else 'raw_input'
    seq_field = args.clustering_field
    sort_file = tempfile.NamedTemporaryFile(dir=args.temp_dir, delete=False)
    sort_handle = open(sort_file.name, 'a')
    count = 0
    for s in seqs:
        sort_handle.write('{} {} {}\n'.format(s['seq_id'], s[seq_field], s[raw_field]))
        count += 1
    sort_handle.close()
    return sort_file.name, count


# =================================
#
#           CLUSTERING
#
# =================================


def initial_clustering(seq_db_path, args):
    logger.info('\n{} clustering with CD-HIT...'.format('Initial UAID' if args.uaid else 'Identity-based'))
    start = time.time()
    conn = sqlite3.connect(seq_db_path)
    seq_db = conn.cursor()
    if args.uaid:
        seqs = seq_db.execute('''SELECT seqs.seq_id, seqs.uaid FROM seqs''')
    else:
        seqs = seq_db.execute('''SELECT seqs.seq_id, seqs.clustering_seq FROM seqs''')
    seqs = [Sequence(s[1], id=s[0]) for s in seqs]
    all_clusters = cluster(seqs, args.identity_threshold, temp_dir=args.temp_dir,
                           quiet=True, max_memory=0, debug=args.debug)
    passed_clusters = [c for c in all_clusters if c.size >= args.min_seqs]
    sizes = [c.size for c in passed_clusters]
    logger.info('{} total clusters identified'.format(len(all_clusters)))
    logger.info('{} clusters meet the minimum size cutoff ({} sequence{})'.format(len(passed_clusters), args.min_seqs, 's' if args.min_seqs > 1 else ''))
    logger.info('The average cluster contains {} sequences; the largest contains {} sequences'.format(round(1. * sum(sizes) / len(sizes), 2), max(sizes)))
    logger.info('Initial clustering took {} seconds\n'.format(round(time.time() - start, 2)))
    conn.commit()
    conn.close()
    return passed_clusters


def process_initial_clusters(initial_clusters, seq_db_path, args):
    logger.info('Calculating {} sequences...'.format('consensus' if args.consensus else 'centroid'))
    consentroids = []
    if args.debug:
        num_clusters = len(initial_clusters)
        update_progress(0, num_clusters, sys.stdout)
        for i, initial_cluster in enumerate(initial_clusters):
            clustering_seqs = retrieve_clustering_seqs(initial_cluster.ids, seq_db_path)
            output_seqs = retrieve_output_seqs(initial_cluster.ids, seq_db_path)
            _consentroids = process_initial_cluster(clustering_seqs, output_seqs, args)
            consentroids.extend(_consentroids)
            update_progress(i + 1, num_clusters, sys.stdout)
    else:
        async_results = []
        p = mp.Pool(maxtasksperchild=10)
        for initial_cluster in initial_clusters:
            clustering_seqs = retrieve_clustering_seqs(initial_cluster.ids, seq_db_path)
            output_seqs = retrieve_output_seqs(initial_cluster.ids, seq_db_path)
            async_results.append(p.apply_async(process_initial_cluster, (clustering_seqs, output_seqs, args)))
        monitor_mp_jobs(async_results)
        for ar in async_results:
            consentroids.extend(ar.get())
        p.close()
        p.join()
    return consentroids


def process_initial_cluster(clustering_seqs, output_seqs, args):
    consentroid_func = calculate_consensus if args.consensus else calculate_centroid
    consentroids = consentroid_func(clustering_seqs, output_seqs, args)
    return consentroids


def calculate_consensus(clustering_seqs, output_seqs, args):
    consensus_seqs = []
    clusters = cluster(clustering_seqs,
                       temp_dir=args.temp_dir,
                       threshold=args.identity_threshold,
                       return_just_seq_ids=True,
                       make_db=False,
                       threads=1,
                       quiet=True,
                       debug=args.debug)
    for clust_ids in clusters:
        _output_seqs = [s for s in output_seqs if s.id in clust_ids]
        consensus = make_consensus(_output_seqs)
        consensus_seqs.append(consensus)
    return consensus_seqs


def make_consensus(seqs):
    if len(seqs) == 1:
        return (Sequence(seqs[0].sequence), 1)
    _aln = muscle(seqs, as_file=True)
    aln = AlignIO.read(open(_aln, 'r'), 'fasta')
    summary_align = AlignInfo.SummaryInfo(aln)
    consensus = summary_align.gap_consensus(threshold=0.51, ambiguous='n')
    consensus_string = str(consensus).replace('-', '')
    consensus_seq = Sequence(consensus_string.upper())
    os.unlink(_aln)
    return (consensus_seq, len(seqs))


def calculate_centroid(clustering_seqs, output_seqs, args):
    threshold = args.identity_threshold if args.uaid else 0.7
    clustering_seqs = retrieve_clustering_seqs(seq_ids, seq_db_path)
    out_file, clust_file = cdhit(clustering_seqs, threads=1, quiet=True, temp_dir=args.temp_dir, debug=args.debug)
    centroids = parse_centroids(clust_file, output_seqs)
    os.unlink(out_file)
    os.unlink(clust_file)
    if args.only_largest_cluster:
        centroids = sorted(centroids, key=lambda x: x[1], reverse=True)[:1]
    return centroids


def parse_centroids(clust_file, output_seqs):
    centroid_ids = []
    sizes = []
    raw_clusters = [c.split('\n') for c in open(clust_file, 'r').read().split('\n>')]
    for rc in raw_clusters:
        size = 0
        for line in rc:
            if line.strip():
                size += 1
            if '*' in line:
                centroid_id = line.split()[2][1:-3]
                centroid_ids.append(centroid_id)
        sizes.append(size)
    centroids = []
    for c, s in zip(centroid_ids, sizes):
        c_seq = [_c for _c in output_seqs if _c.id == c][0]
        centroids.append((c_seq, size))
    return centroids


def process_singleton_clusters(singletons, seq_db_path, args):
    consentroids = []
    num_singletons = len(singletons)
    update_progress(0, num_singletons, sys.stdout)
    for i, singleton in enumerate(singletons):
        seq = retrieve_output_seqs(singleton.ids, seq_db_path)[0]
        seq.id = '{}_1'.format(uuid.uuid4()) if args.consensus else '{}_1'.format(seq.id)
        consentroids.append((seq, 1))
        update_progress(i + 1, num_singletons, sys.stdout)
    sys.stdout.flush()
    logger.info('')
    return consentroids



# =========================================
#
#         PROGRESS AND PRINTING
#
# =========================================



def _print_start_info(args):
    logger.info('')
    logger.info('')
    logger.info('')
    logger.info('-' * 25)
    logger.info('ERROR CORRECTION')
    logger.info('-' * 25)
    _log_params(args)


def _log_params(args):
    logger.info('')
    logger.info('INPUT DB: {}'.format(args.db))
    logger.info('OUTPUT DIR: {}'.format(args.output))
    logger.info('TEMP DIR: {}'.format(args.temp_dir))
    logger.info('UAIDs: {}'.format(args.uaid))
    logger.info('OUTPUT TYPE: {}'.format('Consensus' if args.consensus else 'Centroid'))
    logger.info('CLUSTERING TYPE: {}'.format('sort/uniq' if args.non_redundant else 'CD-HIT'))
    if args.non_redundant:
        logger.info('NON-REDUNDANT CLUSTERING FIELD: {}'.format(args.clustering_field))
        logger.info('IDENTITY THRESHOLD: 1.0')
        logger.info('MIN SEQS: 1')
        return
    elif not args.uaid:
        logger.info('CLUSTERING FIELD: {}'.format(args.clustering_field))
        logger.info('IDENTITY THRESHOLD: {}'.format(args.identity_threshold))
        logger.info('GERMLINES: {}'.format(args.germs))
    else:
        logger.info('PARSE UAIDS: {}'.format('False' if args.parse_uaids == 0 else args.parse_uaids))
    logger.info('MIN SEQS: {}'.format(args.min_seqs))
    logger.info('LARGEST CLUSTER ONLY: {}'.format(args.only_largest_cluster))
    if args.uaid:
        logger.info('INTRA-BIN CLUSTERING THRESHOLD: {}'.format(args.identity_threshold))
        logger.info('INTRA-BIN CLUSTERING FIELD: {}'.format(args.clustering_field))
        logger.info('{} FIELD: {}'.format('CONSENSUS' if args.consensus else 'CENTROID', args.output_field))



def monitor_mp_jobs(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        ready = [ar for ar in results if ar.ready()]
        finished = len(ready)
        update_progress(finished, jobs, sys.stdout)
    sys.stdout.write('\n')


def update_progress(finished, jobs, log, failed=None):
    pct = int(100. * finished / jobs)
    ticks = pct / 2
    spaces = 50 - ticks
    if failed:
        prog_bar = '\r({}/{}) |{}{}|  {}% ({}, {})'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct, finished - failed, failed)
    else:
        prog_bar = '\r({}/{}) |{}{}|  {}%'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct)
    sys.stdout.write(prog_bar)


def write_output(collection, sequences, sizes, collection_start_time, args):
    seq_type = 'consensus' if args.consensus else 'centroid'
    logger.info('Writing {} sequences to output file...'.format(seq_type))
    write_fasta_output(collection, sequences, args)
    if sizes is not None:
        write_stats_output(collection, sizes, args)
    logger.info('{} {} sequences were identified.'.format(len(sequences), seq_type))
    logger.info('{} was processed in {} seconds.'.format(collection,
        round(time.time() - collection_start_time, 2)))
    logger.info('')


def write_nr_output(collection, unique_file, collection_start_time, args):
    logger.info('Writing non-redundant sequences to output file...')
    oname = collection
    if os.path.isfile(collection):
        oname = os.path.basename(collection).rstrip('.json')
    outfile = os.path.join(args.output, '{}_nr.fasta'.format(oname))
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
    logger.info('{} was processed in {} seconds.'.format(collection,
        round(time.time() - collection_start_time, 2)))
    logger.info('')


def write_fasta_output(collection, sequences, args):
    seq_type = 'consensus' if args.consensus else 'centroids'
    oname = collection
    if os.path.isfile(collection):
        oname = os.path.basename(collection).rstrip('.json')
    outfile = os.path.join(args.output, '{}_{}.fasta'.format(oname, seq_type))
    out_handle = open(outfile, 'w')
    out_handle.write('\n'.join([s.fasta for s in sequences]))
    out_handle.close()


def write_stats_output(collection, sizes, args):
    sizes = [int(s) for s in sizes]
    bin_counts = np.bincount(sizes)[1:]
    bins = range(1, len(bin_counts))
    binned_data = zip(bins, bin_counts)
    bin_string = 'Cluster Size\tCount\n'
    bin_string += '\n'.join(['{}\t{}'.format(b[0], b[1]) for b in binned_data])
    oname = collection
    if os.path.isfile(collection):
        oname = os.path.basename(collection).rstrip('.json')
    outfile = os.path.join(args.output, '{}_cluster_sizes.txt'.format(oname))
    out_handle = open(outfile, 'w')
    out_handle.write(bin_string)
    out_handle.close()



def print_collection_info(collection, sample_name):
    logger.info('')
    logger.info('')
    logger.info(collection)
    logger.info('-' * len(collection))
    logger.info('SAMPLE NAME: {}'.format(sample_name))


def countdown(args):
    start = time.time()
    h = int(args.sleep / 60)
    m = int(args.sleep % 60)
    hz = '0' if len(str(h)) == 1 else ''
    mz = '0' if len(str(m)) == 1 else ''
    print('Countdown: {}{}:{}{}:{}'.format(hz, h, mz, m, '00'), end='')
    sys.stdout.flush()
    done = False
    while not done:
        time.sleep(0.25)
        elapsed = int(time.time() - start)
        if elapsed < args.sleep * 60:
            h = int((args.sleep - elapsed / 60.) / 60)
            m = int((args.sleep - elapsed / 60.) % 60)
            s = int((60 * args.sleep - elapsed) % 60)
            hz = '0' if len(str(h)) == 1 else ''
            mz = '0' if len(str(m)) == 1 else ''
            sz = '0' if len(str(s)) == 1 else ''
            print('\rCountdown: {}{}:{}{}:{}{}'.format(hz, h, mz, m, sz, s), end='')
            sys.stdout.flush()
        else:
            print('\rCountdown: 00:00:00')
            done = True


def run(**kwargs):
    '''
    Corrects antibody reads using UAIDs (molecular barcodes) or identity-based clustering.

    Either ``json`` or ``db`` is required.


    Args:

        db (str): Name of a MongoDB database to query.

        collection (str): Name of a MongoDB collection. If not provided, all
            collections in ``db`` will be iteratively processed.

        json: Can be one of two things:

            1. Path to a JSON file, containing sequence data annotated by AbStar.

            2. Path to a directory, containing one or more JSON files of
                AbStar-annotated data.

        output (str): Path to the output directory, into which corrected FASTA
            files will be deposited. If it does not exist, it will be created.

        log (str): Path to the log file. If the parent directory doesn't exist,
            it will be created.

        ip (str): IP address of the MongoDB server. Default is ``localhost``.

        port (int): Port of the MongoDB server. Default is ``27017``.

        user (str): Username with which to connect to the MongoDB database. If either
            of ``user`` or ``password`` is not provided, the connection to the MongoDB
            database will be attempted without authentication.

        password (str): Password with which to connect to the MongoDB database. If either
            of ``user`` or ``password`` is not provided, the connection to the MongoDB
            database will be attempted without authentication.

        min_seqs (int): Minimum number of sequences for a centroid/consensus sequence to be
            calculated. After clustering (either by identity or using UAIDs), clusters with
            at least ``min_seqs`` sequences will be retained for consensus/centroid calculation.
            Default is ``1``.

        uaid (bool): If ``True``, use Unique Antibody IDs (UAIDs) for error correction.
            Sequences will be binned by UAID and the sequences in each bin will be used to
            compute a centroid or consensus sequence. If ``False``, sequences will be clustered
            by identity and each cluster will be used for consensus/centroid determination.

        parse_uaids (int): If UAIDs haven't been pre-parsed by AbStar, indicate the length of the
            UAID sequence (in nucleotides) and the UAIDs will be parsed during correction. If
            ``parse_uaids`` is negative, the UAID will be parsed from the end of the sequence.
            Default is ``0``, which does not parse a UAID sequence.

        consensus (bool): If ``True``, consensus sequences are calculated. If ``False``, centroid
            sequences are calculated. Default is ``True``.

        identity_threshold (float): Identity threshold, if clustering by identity (not UAIDs).
            Must be a float between 0 and 1. Default is 0.975.

        only_largest_cluster (bool): When clustering using UAIDs, there is a some probability that
            different sequences get labeled with the same UAID. To limit incorrect consensus/centroid
            calculation, sequences in each UAID bin are clustered using ``identity_threshold`` before
            calculating consens/centroid sequences. By default, all UAID clusters that meet the
            ``min_seqs`` size threshold are used to generate consensus/centroid sequences. If that
            behavior is not desired, setting ``only_largest_cluster`` to ``True`` results in only
            the largest UAID cluster being used to generate centroid/consensus sequences.

        nr (bool): If ``True``, a non-redundant sequence dataset will be generated using ``sort | uniq``.
            This is much faster than normal sequence clustering with CD-HIT, but can only be performed at an
            identity threshold of 100%.

            .. note:

                Using ``nr`` may produce different results than clustering sequences with ``identity_threshold``
                set to ``1.0``. This is because sequences of different lengths that are otherwise identical
                will not be collapsed when using ``nr`` but will be collapsed using normal clustering.

        germs (str): Path to a file containing germline V-gene sequences. When clustering with ``min_seqs``
            equal to 2, the appropriate germline sequence will be added to the alignment to serve as a
            consensus tiebreaker.

        aa (bool): If ``True``, perform sequence clustering (either using ``identity_threshold`` or ``nr``)
            using amino acid sequences. Default is ``False``, which performs clustering using nucleotide
            sequences.

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
    _print_start_info(args)
    if args.sleep:
        countdown(args)
    for d in [args.output, args.temp_dir]:
        make_dir(d)
    if args.consensus and args.germs:
        germs = parse_germs(args.germs)
    else:
        germs = args.germs
    # check whether JSON files have been passed
    if args.json is not None and all([args.db is None, args.collection is None]):
    	if os.path.isfile(args.json) and args.json.endswith('.json'):
    		collections = [args.json, ]
    	else:
        	collections = list_files(args.json, extension='json')
        db = None
        sample_names = [os.path.basename(c).replace('.json', '') for c in collections]
    # otherwise, get sequences from MongoDB
    else:
        db = mongodb.get_db(args.db, args.ip, args.port, args.user, args.password)
        collections = mongodb.get_collections(db, collection=args.collection)
        sample_names = collections
    for collection, sample_name in zip(collections, sample_names):
        collection_start = time.time()
        print_collection_info(collection, sample_name)
        if args.non_redundant:
            seqs = get_seqs(db, collection, args, make_seq_db=False)
            unique_file = unix_sort_unique(seqs, args)
            write_nr_output(collection, unique_file, collection_start, args)
        else:
            seq_db_path = get_seqs(db, collection, args)
            initial_clusters = initial_clustering(seq_db_path, args)
            if args.min_seqs == 1:
                singletons = [ic for ic in initial_clusters if ic.size == 1]
                initial_clusters = [ic for ic in initial_clusters if ic.size > 1]
                logger.info('{} clusters contained only a single sequence. Processing singletons...'.format(len(singletons)))
                singleton_consentroids = process_singleton_clusters(singletons, seq_db_path, args)
                logger.info('')
            else:
                singleton_consentroids = []
            consentroids = process_initial_clusters(initial_clusters, seq_db_path, args)
            consentroids += singleton_consentroids
            sequences, sizes = zip(*consentroids)
            write_output(sample_name, sequences, sizes, collection_start, args)
            for ic in initial_clusters:
                ic.cleanup()
            remove_sqlite_db(args)


if __name__ == '__main__':
    parser = parse_args()
    args = parser.parse_args()
    logfile = args.log if args.log else os.path.join(args.output, 'abcorrect.log')
    log.setup_logging(logfile)
    logger = log.get_logger('abcorrect')
    main(args)
