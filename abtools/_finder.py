#!/usr/bin/python
# filename: _finder.py


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

import multiprocessing as mp
import platform
import os
import subprocess as sp
import sys
import tempfile
from threading import Thread
import time

import numpy as np
import pandas as pd

from pymongo import MongoClient

from Bio import SeqIO

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns

from abtools import color, log, mongodb
from abtools.utils import progbar


def parse_args():
    import argparse
    parser = argparse.ArgumentParser("For a MongoDB collection, plots the germline divergence against the sequence identity to a given 'subject' sequence.")
    parser.add_argument('-d', '--database', dest='db', required=True,
                        help="Name of the MongoDB database to query. Required.")
    parser.add_argument('-c', '--collection', dest='collection', default=None,
                        help="Name of the MongoDB collection to query. \
                        If not provided, all collections in the given database will be processed iteratively.")
    parser.add_argument('--collection-prefix', dest='collection_prefix', default=None,
                        help="If supplied, will iteratively process only collections beginning with <collection_prefix>.")
    parser.add_argument('-o', '--output', dest='output_dir', default=None,
                        help="Output directory figure files. If not provided, figures will not be generated. \
                        Directory will be created if it does not already exist.")
    parser.add_argument('-t', '--temp', dest='temp_dir', required=True,
                        help="Directory for temporary storage. \
                        Will be created if it does not already exist. Required.")
    parser.add_argument('-l', '--log', dest='log', default=None,
                        help="The log file, to which the blast_parse log info will be written. \
                        Default is <output>/abfinder.log.")
    parser.add_argument('-C', '--cluster', dest="cluster", default=False, action='store_true',
                        help="Use if performing computation on a Celery cluster. \
                        If set, input files will be split into many subfiles and passed to a Celery queue. \
                        If not set, input files will still be split, \
                        but will be distributed to local processors using multiprocessing.")
    parser.add_argument('-i', '--ip', dest='ip', default='localhost',
                        help="The IP address for the MongoDB server.  \
                        Defaults to 'localhost'.")
    parser.add_argument('--port', dest='port', default=27017,
                        help="The port for the MongoDB server. Defaults to '27017'.")
    parser.add_argument('-u', '--user', dest='user', default=None,
                        help="Username for the MongoDB server. Not used if not provided.")
    parser.add_argument('-p', '--password', dest='password', default=None,
                        help="Password for the MongoDB server. Not used if not provided.")
    parser.add_argument('-s', '--standard', dest='standard', required=True,
                        help='Path to a file containing the standard sequence(s) for which \
                        identity/divergence will be calculated, in FASTA format. \
                        All sequences in the standard file will iteratively processed. Required')
    parser.add_argument('-q', '--chain', dest='chain', default='heavy',
                        choices=['heavy', 'kappa', 'lambda', 'light'],
                        help="The chain type of the subject sequence. \
                        Options are 'heavy', 'kappa', 'lambda' and 'light'. \
                        Default is 'heavy'.")
    parser.add_argument('-n', '--no_update', dest='update', action='store_false', default=True,
                        help="Does not update the MongoDB with AbFinder info. \
                        Can save some time if the identity calculations aren't needed again.")
    parser.add_argument('--no_figure', dest='make_figure', action='store_false', default=True,
                        help="Does not make the identity/divergence figure. \
                        Useful if you don't want the figure, just the identity info written to the database.")
    parser.add_argument('--single-process-update', dest='single_process_update', action='store_true', default=False,
                        help="Perform the MongoDB update using a single process (without multiprocessing).")
    parser.add_argument('--update-threads', dest='update_threads', type=int, default=25,
                        help="Number of threads to use when update the MongoDB database. Default is 25.")
    parser.add_argument('-N', '--nucleotide', dest='is_aa', action='store_false', default=True,
                        help="Use nucleotide sequences for alignment. Default is amino acid sequences. \
                        Ensure standard format matches.")
    parser.add_argument('-x', '--xmin', dest='x_min', type=int, default=-1,
                        help="Minimum X-axis (germline divergence) value for the AbCompare plot. Default is -1.")
    parser.add_argument('-X', '--xmax', dest='x_max', type=int, default=35,
                        help="Maximum X-axis (germline divergence) value for the AbCompare plot. Default is 35.")
    parser.add_argument('-y', '--ymin', dest='y_min', type=int, default=65,
                        help="Minimum Y-axis (mAb identity) value for the AbCompare plot. Default is 65.")
    parser.add_argument('-Y', '--ymax', dest='y_max', type=int, default=101,
                        help="Maximum Y-axis (mAb identity) value for the AbCompare plot. Default is 101.")
    parser.add_argument('-g', '--gridsize', dest='gridsize', type=int, default=0,
                        help="Gridsize for the AbFinder hexbin plot. \
                        Default is 36 for amino acid sequences and 50 for nucleotide sequences.")
    parser.add_argument('--colormap', dest='colormap', default='Blues',
                        help="Colormap to be used in the AbFinder hexbin plots. \
                        Can accept a matplotlib cmap or the name of one of matplotlib's builtin cmaps. \
                        Default is 'Blues'.")
    parser.add_argument('--mincount', dest='mincount', default=3, type=int,
                        help="Minimum number of sequences in a hexbin for that hexbin to be colored. \
                        Default is 3.")
    parser.add_argument('--skip-padding', dest='remove_padding', default=True, action='store_false',
                        help="If set, will not remove padding field from MongoDB.")
    parser.add_argument('-D', '--debug', dest="debug", action='store_true', default=False,
                        help="If set, will write all failed/exception sequences to file \
                        and should give more informative errors.")
    return parser


class Args(object):
    def __init__(self, db=None, collection=None,
                 output=None, temp=None, log=None, cluster=False,
                 ip='localhost', port=27017, user=None, password=None, update=True,
                 standard=None, chain='heavy', is_aa=True,
                 x_min=-1, x_max=35, y_min=65, y_max=101, gridsize=0, mincount=3,
                 colormap='Blues', debug=False):
        super(Args, self).__init__()
        if not all([db, output, temp, standard]):
            err = 'You must provide a MongoDB database name, output and temp directories, \
                and a file containing one or more comparison (standard) sequences in FASTA format.'
            raise RuntimeError(err)
        self.db = db
        self.collection = collection
        self.output_dir = output
        self.temp_dir = temp
        self.log = log
        self.cluster = bool(cluster)
        self.ip = ip
        self.port = int(port)
        self.user = user
        self.password = password
        self.standard = standard
        if chain not in ['heavy', 'kappa', 'lambda', 'light']:
            err = 'Please select an appropriate chain. \
                Valid choices are: heavy, light, kappa and lambda.'
            raise RuntimeError(err)
        self.chain = chain
        self.update = bool(update)
        self.is_aa = bool(is_aa)
        self.x_min = int(x_min)
        self.x_max = int(x_max)
        self.y_min = int(y_min)
        self.y_max = int(y_max)
        self.gridsize = int(gridsize)
        mincount = int(mincount)
        self.colormap = colormap
        self.debug = bool(debug)



# ================================================
#
#             FILES AND DIRECTORIES
#
# ================================================



def make_directories(args):
    for d in [args.output_dir, args.temp_dir]:
        if d:
            _make_direc(d, args)


def _make_direc(d, args):
    if not os.path.exists(d):
        os.makedirs(d)
    if args.cluster:
        cmd = 'sudo chmod 777 {}'.format(d)
        p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()


def get_standards(args):
    standards = []
    for s in SeqIO.parse(open(args.standard, 'r'), 'fasta'):
        standards.append(s)
    return standards


def get_chain(args):
    if args.chain == 'light':
        return ['kappa', 'lambda']
    return [args.chain, ]


def get_sequences(db, collection, temp_dir, args):
    files = []
    fastas = []
    chunksize = 1000
    seq_counter = 0
    total_seq_counter = 0
    query_results = query(db, collection, args)
    iden_field = 'aa_identity' if args.is_aa else 'nt_identity'
    vdj_field = 'vdj_aa' if args.is_aa else 'vdj_nt'
    for seq in query_results:
        fastas.append('>{}_{}\n{}'.format(seq['seq_id'], seq[iden_field]['v'], seq[vdj_field]))
        seq_counter += 1
        total_seq_counter += 1
        if seq_counter == chunksize:
            files.append(write_to_temp_file(fastas, temp_dir))
            fastas = []
            seq_counter = 0
    if fastas:
        files.append(write_to_temp_file(fastas, temp_dir))
    return files


def write_to_temp_file(fastas, temp_dir):
    tfile = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
    tfile.write('\n'.join(fastas))
    tfile.close()
    return tfile.name


def clean_up(files):
    for f in files:
        os.unlink(f)



# ================================================
#
#                   MONGO
#
# ================================================



def query(db, collection, args):
    coll = db[collection]
    chain = get_chain(args)
    mongodb.index(db, collection, ['chain'])
    print_query_info()
    iden_field = 'aa_identity.v' if args.is_aa else 'nt_identity.v'
    vdj_field = 'vdj_aa' if args.is_aa else 'vdj_nt'
    return coll.find({'chain': {'$in': chain}, 'prod': 'yes'}, {'_id': 0, 'seq_id': 1, iden_field: 1, vdj_field: 1})


def chunker(l, n):
    'Generator that produces n-length chunks from iterable l.'
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def update_db(db, standard, scores, collection, args):
    db = mongodb.get_db(args.db, args.ip, args.port, args.user, args.password)
    print_index_info()
    mongodb.index(db, collection, ['seq_id'])
    print_update_info()
    start = time.time()
    conn = mongodb.get_connection(args.ip, args.port,
        args.user, args.password)
    mongo_version = conn.server_info()['version']
    standard = standard.replace('.', '_')
    g = scores.groupby('identity')
    groups = regroup(g.groups)


    for g in range(0, len(groups), args.update_threads):
        tlist = []
        for group in groups[g:g + args.update_threads]:
            t = Thread(target=update, args=(db, collection, group, standard, mongo_version, args))
            t.start()
            tlist.append(t)
        for t in tlist:
            t.join()
        progbar.progress_bar(g + args.update_threads, len(groups))


    # if platform.system().lower() == 'darwin' or args.debug or args.single_process_update:
    #     for i, group in enumerate(groups):
    #         update(db, collection, group, standard, mongo_version, args)
    #         progbar.progress_bar(i, len(groups))
    # else:
    #     p = mp.Pool(processes=25)
    #     async_results = []
    #     for group in groups:
    #         async_results.append(p.apply_async(update, args=(db, collection, group, standard, mongo_version, args)))
    #     monitor_update(async_results)
    #     p.close()
    #     p.join()
    print('')
    run_time = time.time() - start
    logger.info('Updating took {} seconds. ({} sequences per second)'.format(round(run_time, 2),
        round(len(scores) / run_time, 1)))


def update(db, collection, data, standard, version, args):
    db = mongodb.get_db(args.db, args.ip, args.port, args.user, args.password)
    coll = db[collection]
    score = data[0]
    ids = data[1]
    mab_id_field = 'mab_identity_aa' if args.is_aa else 'mab_identity_nt'
    if int(version.split('.')[0]) < 3:
        result = coll.update({'seq_id': {'$in': ids}},
                    {'$set': {'{}.{}'.format(mab_id_field, standard.lower()): float(score)}},
                    multi=True)
    else:
        result = coll.update_many({'seq_id': {'$in': ids}},
                         {'$set': {'{}.{}'.format(mab_id_field, standard.lower()): float(score)}})

        if args.debug:
            print('matched: {}'.format(result.matched_count))
            print('modified: {}'.format(result.modified_count))


def monitor_update(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        finished = len([r for r in results if r.ready()])
        progbar.progress_bar(finished, jobs)
    progbar.progress_bar(finished, jobs)


def regroup(oldgs):
    newgs = []
    for og in oldgs:
        if len(oldgs[og]) <= 500:
            newgs.append((og, oldgs[og]))
        else:
            for ng in chunker(oldgs[og], 500):
                newgs.append((og, ng))
    return newgs




# ================================================
#
#                  FIGURES
#
# ================================================



def make_figure(standard_id, scores, collection, args):
    print_fig_info()
    sns.set_style('white')
    fig_file = os.path.join(args.output_dir, '{0}_{1}_{2}.pdf'.format(args.db, collection, standard_id))
    x = list(scores['germ_divergence'].values)
    y = list(scores['identity'].values)
    xy_vals = zip(x, y)
    trunc_xy_vals = [v for v in xy_vals if v[0] <= args.x_max and v[1] >= args.y_min]
    x = [v[0] for v in trunc_xy_vals]
    y = [v[1] for v in trunc_xy_vals]
    # To make sure the gridsize is correct (since it's based on the actual values)
    # I need to add a single value near the max and min of each axis.
    # They're added just outside the visible plot, so there's no effect on the plot.
    x.extend([args.x_min - 1, args.x_max + 1])
    y.extend([args.y_min - 1, args.y_max + 1])
    # plot params
    cmap = color.get_cmap(args.colormap)
    plt.subplots_adjust(hspace=0.95)
    plt.subplot(111)
    plt.hexbin(x, y, bins='log', cmap=cmap, mincnt=3, gridsize=set_gridsize(args))
    plt.title(standard_id, fontsize=18)
    # set and label axes
    plt.axis([args.x_min, args.x_max, args.y_min, args.y_max])
    plt.xlabel('Germline divergence')
    plt.ylabel('{0} identity'.format(standard_id))
    # make and label the colorbar
    cb = plt.colorbar()
    cb.set_label('Sequence count (log10)', labelpad=10)
    # save figure and close
    plt.savefig(fig_file)
    plt.close()


def set_gridsize(args):
    if args.gridsize:
        return args.gridsize
    elif args.is_aa:
        return 36
    return 50




# ================================================
#
#                  PRINTING
#
# ================================================



def print_abfinder_start():
    logger.info('')
    logger.info('')
    logger.info('')
    logger.info('-' * 25)
    logger.info('ABFINDER')
    logger.info('-' * 25)


def print_standards_info(standards):
    logger.info('')
    logger.info('Found {} standard sequence(s):'.format(len(standards)))
    logger.info(', '.join([s.id for s in standards]))


def print_collections_info(collections):
    logger.info('')
    logger.info('Found {} collection(s):'.format(len(collections)))
    logger.info(', '.join(collections))


def print_single_standard(standard):
    standard_id_string = '{}'.format(standard.id)
    logger.info('')
    logger.info(standard_id_string)
    logger.info('-' * len(standard_id_string))


def print_single_collection(collection):
    logger.info('')
    logger.info('')
    logger.info(collection)
    logger.info('-' * len(collection))


def print_query_info():
    logger.info('Querying for comparison sequences...')


def print_remove_padding():
    logger.info('')
    logger.info('Removing MongoDB padding...')


def print_fig_info():
    logger.info('Making identity/divergence figure...')


def print_index_info():
    logger.info('Indexing the MongoDB collection...')


def print_update_info():
    logger.info('Updating the MongoDB database with identity scores:')




# ================================================
#
#                IDENTITY JOBS
#
# ================================================



def run_jobs(files, standard, args):
    logger.info('Running AbCompare...')
    if args.cluster:
        return _run_jobs_via_celery(files, standard, args)
    else:
        return _run_jobs_via_multiprocessing(files, standard, args)


def _run_jobs_via_multiprocessing(files, standard, args):
    from abtools.queue.tasks import identity
    results = []
    if args.debug:
        for f in files:
            results.extend(identity(f, standard, args.is_aa, args.debug))
    else:
        p = mp.Pool()
        async_results = []
        for f in files:
            async_results.append(p.apply_async(identity, (f, standard, args.is_aa)))
        monitor_mp_jobs(async_results)
        for a in async_results:
            results.extend(a.get())
        p.close()
        p.join()
    ids = [r[0] for r in results]
    identities = pd.Series([r[1] for r in results], index=ids)
    divergences = pd.Series([100. - r[2] for r in results], index=ids)
    d = {'identity': identities, 'germ_divergence': divergences}
    df = pd.DataFrame(d)
    return df


def monitor_mp_jobs(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        ready = [ar for ar in results if ar.ready()]
        finished = len(ready)
        update_progress(finished, jobs)
    print('')


def _run_jobs_via_celery(files, standard, args):
    from abtools.queue.tasks import identity
    async_results = []
    for f in files:
        async_results.append(identity.delay(f, standard, args.is_aa))
    succeeded, failed = monitor_celery_jobs(async_results)
    scores = []
    for s in succeeded:
        scores.extend(s.get())
    ids = [r[0] for r in scores]
    identities = pd.Series([r[1] for r in scores], index=ids)
    divergences = pd.Series([r[2] for r in scores], index=ids)
    d = {'identity': identities, 'germ_divergence': divergences}
    df = pd.DataFrame(d)
    return df


def monitor_celery_jobs(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        succeeded = [ar for ar in results if ar.successful()]
        failed = [ar for ar in results if ar.failed()]
        finished = len(succeeded) + len(failed)
        update_progress(finished, jobs, failed=len(failed))
    print('')
    return succeeded, failed


def update_progress(finished, jobs):
    pct = int(100. * finished / jobs)
    ticks = pct / 2
    spaces = 50 - ticks
    prog_bar = '\r({}/{}) |{}{}|  {}%'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct)
    sys.stdout.write(prog_bar)
    sys.stdout.flush()


def run(**kwargs):
    '''
    Mines NGS datasets for identity to known antibody sequences.

    All of ``db``, ``output``, ``temp`` and ``standard`` are required.


    Args:

        db (str): Name of a MongoDB database to query.

        collection (str): Name of a MongoDB collection. If not provided, all collections
            in ``db`` will be processed iteratively.

        output_dir (str): Path to the output directory, into which identity/divergence
            figures will be deposited.

        temp_dir (str): Path to a temporary directory.

        log (str): Path to a log file. If not provided, log information will not be retained.

        ip (str): IP address of the MongoDB server. Default is ``localhost``.

        port (str): Port of the MongoDB server. Default is ``27017``.

        user (str): Username with which to connect to the MongoDB database. If either
            of ``user`` or ``password`` is not provided, the connection to the MongoDB
            database will be attempted without authentication.

        password (str): Password with which to connect to the MongoDB database. If either
            of ``user`` or ``password`` is not provided, the connection to the MongoDB
            database will be attempted without authentication.

        standard (path): Path to a FASTA-formatted file containing one or more 'standard'
            sequences, against which the NGS sequences will be compared.

        chain (str): Antibody chain. Choices are 'heavy', 'kappa', 'lambda', and 'light'.
            Default is 'heavy'. Only NGS sequences matching ``chain`` (with 'light' covering
            both 'kappa' and 'lambda') will be compared to the ``standard`` sequences.

        update (bool): If ``True``, the MongoDB record for each NGS sequence will be updated
            with identity information for each standard. If ``False``, the updated is skipped.
            Default is ``True``.

        is_aa (bool): If ``True``, the ``standard`` sequences are amino acid sequences. If
            ``False``, they are nucleotide seqeunces. Default is ``False``.

        x_min (int): Minimum x-axis value on identity/divergence plots.

        x_max (int): Maximum x-axis value on identity/divergence plots.

        y_min (int): Minimum y-axis value on identity/divergence plots.

        y_max (int): Maximum y-axis value on identity/divergence plots.

        gridsize (int): Relative size of hexbin grids.

        mincount (int): Minimum number of sequences in a hexbin for the bin to be colored.
            Default is 3.

        colormap (str, colormap): Colormap to be used for identity/divergence plots.
            Default is ``Blues``.

        debug (bool): If ``True``, more verbose logging.
   '''
    args = Args(**kwargs)
    global logger
    logger = log.get_logger('abfinder')
    main(args)


def run_standalone(args):
    logfile = args.log if args.log else os.path.join(args.output_dir, 'abfinder.log')
    log.setup_logging(logfile)
    global logger
    logger = log.get_logger('abfinder')
    main(args)


def main(args):
    print_abfinder_start()
    db = mongodb.get_db(args.db, args.ip, args.port,
        args.user, args.password)
    make_directories(args)
    standards = get_standards(args)
    print_standards_info(standards)
    collections = mongodb.get_collections(db, args.collection, prefix=args.collection_prefix)
    print_collections_info(collections)
    for collection in collections:
        indexed = False
        print_single_collection(collection)
        if args.remove_padding:
            print_remove_padding()
            mongodb.remove_padding(db, collection)
        seq_files = get_sequences(db, collection, args.temp_dir, args)
        for standard in standards:
            print_single_standard(standard)
            scores = run_jobs(seq_files, standard, args)
            if args.output_dir:
                make_figure(standard.id, scores, collection, args)
            if args.update:
                if not indexed:
                    mongodb.index(db, collection, 'seq_id')
                    indexed = True
                update_db(db, standard.id, scores, collection, args)
        clean_up(seq_files)


if __name__ == '__main__':
    parser = parse_args()
    args = parser.parse_args()
    logfile = args.log if args.log else os.path.join(args.output_dir, 'abfinder.log')
    log.setup_logging(logfile)
    logger = log.get_logger('abfinder')
    main(args)
