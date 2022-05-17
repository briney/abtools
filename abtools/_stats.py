#!/usr/bin/env python
# filename: _stats.py


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

import argparse
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from abstar.core.germline import get_germlines
from abutils.io import read_sequences
from abutils.utils import log
from abutils.utils import mongodb
from abutils.utils.pipeline import make_dir, list_files


def parse_args():
    parser = argparse.ArgumentParser("Computes and plots basic repertoire information from one or more antibody NGS datasets.")
    parser.add_argument('-o', '--output', dest='output', required=True,
                        help="Output directory for figure and data files. Required.")
    parser.add_argument('-t', '--temp', dest='temp_dir', required=True,
                        help="Directory for temporary files. Required.")
    parser.add_argument('-l', '--log', dest='log',
                        help="Location for the AbStats log file. \
                        If not provided, will be written to <output>/abstats.log.")
    parser.add_argument('-d', '--database', dest='db', required=False,
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
    parser.add_argument('-i', '--ip', dest='ip', default='localhost',
                        help="IP address for the MongoDB server.  Defaults to 'localhost'.")
    parser.add_argument('--port', dest='port', default=27017, type=int,
                        help="Port for the MongoDB server.  Defaults to 27017.")
    parser.add_argument('-u', '--user', dest='user', default=None,
                        help="Username for the MongoDB server. Not used if not provided.")
    parser.add_argument('-p', '--password', dest='password', default=None,
                        help="Password for the MongoDB server. Not used if not provided.")
    parser.add_argument('-V', '--var', dest='var_plot', default=None,
                        help="Plot the distribution of variable families or genes. \
                        Options are 'fam', 'gene', or 'both'. \
                        If not provided, no plot will be made.")
    parser.add_argument('-D', '--div', dest='div_plot', default=None,
                        help="Plot the distribution of diversity families or genes. \
                        Options are 'fam', 'gene', or 'both'. \
                        If not provided, no plot will be made.")
    parser.add_argument('-J', '--join', dest='join_plot', default=None,
                        help="Plot the distribution of joining families or genes. \
                        Only option is 'gene'. \
                        If not provided, no plot will be made.")
    parser.add_argument('-3', '--cdr3', dest='cdr3_plot', default=None,
                        help="Plot the distribution of CDR3 lengths, as nucleotide or amino acid. \
                        Options are 'nt', 'aa', or 'both'. \
                        If not provided, no plot will be made.")
    parser.add_argument('-H', '--heatmap', dest='heatmap', default=False, action='store_true',
                        help="If set, generates a heatmap (or quiltplot, if you're a novice) of combined variable and joining gene use. \
                        If not provided, no plot will be made.")
    parser.add_argument('-I', '--isotype', dest='isotype_plot', default=False, action='store_true',
                        help="Plot the isotype frequency. \
                        Requires that isotypes have been identified using AbStar, or the 'isotype' field is present in the given MongoDB collection. \
                        If not provided, no plot will be made.")
    parser.add_argument('-C', '--chain', dest='chain', default='heavy', choices=['heavy', 'kappa', 'lambda'],
                        help="Select the antibody chain to analyze. \
                        Options are 'heavy', 'kappa', and 'lambda'. \
                        Default is 'heavy'.")
    parser.add_argument('-s', '--species', dest='species', default='human', choices=['human', 'mouse', 'macaque'],
                        help="Species of the sequence data. Choices are 'human', 'mouse', and 'macaque'. \
                        Default is 'human'.")
    parser.add_argument('--id_key', dest='id_key', default='sequence_id',
                        help="The field containing the sequence ID. Default is 'sequence_id'.")
    parser.add_argument('--raw_key', dest='raw_key', default='raw_input', required=True,
                        help="The field containing the raw input sequence, used to pass sequences. \
                        Default is 'raw_input'.")
    parser.add_argument('--debug', dest='debug', action='store_true', default=False,
                        help="If set, will run in debug mode.")
    return parser



class Args(object):
    """docstring for Args"""
    def __init__(self, output=None, temp=None, log=None,
        db=None, collection=None, json=None, tabular=None, sep='\t',
        ip='localhost', port=27017, user=None, password=None,
        var_plot=None, div_plot=None, join_plot=None,
        cdr3_plot=None, heatmap=False, isotype_plot=False, chain='heavy',
        id_key='sequence_id', raw_key='raw_input', debug=False):
        super(Args, self).__init__()
        if all([db is None, json is None, tabular is None]):
            print("\nERROR: One of 'db', 'json' or 'tabular' must be supplied.\n")
            sys.exit(1)
        if not all([output, temp]):
            print('\nERROR: --output and --temp are all required options.\n')
            sys.exit(1)

        self.output = output
        self.temp = temp
        self.log = log
        self.db = db
        self.collection = collection
        self.json = json
        self.tabular = tabular
        self.sep = sep
        self.ip = ip
        self.port = int(port)
        self.user = user
        self.password = password
        self.var_plot = var_plot
        self.div_plot = div_plot
        self.join_plot = join_plot
        self.cdr3_plot = cdr3_plot
        self.heatmap = heatmap
        self.isotype_plot = isotype_plot
        self.chain = chain
        self.id_key = id_key
        self.raw_key = raw_key
        self.debug = debug

def get_seqs(input, args):
    # read JSON data
    fields = [args.id_key, args.raw_key, "cdr3_len", "v_gene", "d_gene", "j_gene"]
    if args.json is not None:
        logger.info('reading input file...')
        seqs = read_sequences(file=input, format='json',
                              id_key=args.id_key,
                              sequence_key=args.raw_key,
                              fields=fields)
    # read tabular data
    elif args.tabular is not None:
        logger.info('reading input file...')
        seqs = read_sequences(file=input, format='tabular', sep=args.sep,
                              id_key=args.id_key,
                              sequence_key=args.raw_key,
                              fields=fields)
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
                              sequence_key=args.raw_key)
    else:
        seqs = []
    # make sure all of the necessary keys are present in the annotated data
    required_keys = [args.id_key, args.raw_key, "cdr3_len", "v_gene", "d_gene", "j_gene"]
    check_sequence_keys(seqs, required_keys)
    logger.info(f'found {len(seqs)} input sequences')
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

# def get_db(db, ip='localhost', port=27017, user=None, password=None):
#     if all([user is not None, password is not None]):
#         pwd = urllib.parse.quote_plus(password)
#         uri = 'mongodb://{}:{}@{}:{}'.format(user, pwd, ip, port)
#         conn = MongoClient(uri)
#     else:
#         conn = MongoClient(ip, port)
#     return conn[db]
#
#
# def get_collections(db):
#     if args.collection:
#         return [args.collection, ]
#     collections = db.collection_names(include_system_collections=False)
#     return sorted(collections)
#
#
#
# def query(db, collection, chain):
#     print('Getting sequences from MongoDB...', end='')
#     sys.stdout.flush()
#     c = db[collection]
#     results = c.find({'chain': chain, 'prod': 'yes'},
#                      {'_id': 0,
#                       'v_gene.gene': 1, 'd_gene.gene': 1, 'j_gene.gene': 1,
#                       'v_gene.fam': 1, 'd_gene.fam': 1, 'j_gene.fam': 1,
#                       'cdr3_len': 1, 'isotype': 1})
#     seqs = [r for r in results if '/OR' not in r['v_gene']['gene']]
#     print('Done.\nFound {} sequences\n'.format(len(seqs)))
#     return seqs


def aggregate(data, norm=True, sort_by='value', keys=None):
    '''
    Counts the number of occurrences of each item in 'data'.

    Inputs
    data: a list of values.
    norm: normalize the resulting counts (as percent)
    sort_by: how to sort the returned data. Options are 'value' and 'count'.

    Output
    a non-redundant list of values (from 'data') and a list of counts.
    '''
    if keys:
        vdict = {k: 0 for k in keys}
        for d in data:
            if d in keys:
                vdict[d] += 1
    else:
        vdict = {}
        for d in data:
            vdict[d] = vdict[d] + 1 if d in vdict else 1
    vals = [(k, v) for k, v in vdict.items()]
    if sort_by == 'value':
        vals.sort(key=lambda x: x[0])
    else:
        vals.sort(key=lambda x: x[1])
    xs = [v[0] for v in vals]
    if norm:
        raw_y = [v[1] for v in vals]
        total_y = sum(raw_y)
        ys = [100. * y / total_y for y in raw_y]
    else:
        ys = [v[1] for v in vals]
    return xs, ys




def cdr3_plot(seqs, collection, make_plot, chain, output_dir):
    if not make_plot:
        return None
    max_len = 40 if chain == 'heavy' else 20
    cdr3s = [s['cdr3_len'] for s in seqs if s['cdr3_len'] > 0 and s['cdr3_len'] <= max_len]
    x, y = aggregate(cdr3s, keys=list(range(1, max_len + 1)))
    color = sns.hls_palette(7)[4]
    plot_file = os.path.join(output_dir, '{}_{}_cdr3_lengths.pdf'.format(collection, chain))
    x_title = 'CDR3 Length (AA)'
    y_title = 'Frequency (%)'
    size = (9, 4) if chain == 'heavy' else (6, 4)
    make_barplot([str(i) for i in x], y,
                 color,
                 plot_file,
                 xlabel=x_title,
                 ylabel=y_title,
                 grid=True,
                 size=size,
                 xfontsize=7)


def germline_plot(seqs, gene, collection, output_dir, level, species, chain):
    germs = get_germlines(species, gene, chain=chain)
    if not level:
        return None
    level = ['fam', 'gene'] if level == 'both' else [level, ]
    for l in level:
        if l == 'fam':
            keys = list(set([g.name.split('-')[0] for g in germs]))
            size = (6, 4)
        else:
            keys = list(set([g.name.split('*')[0] for g in germs]))
            if gene == 'J':
                size = (6, 4)
            elif gene == 'D':
                size = (8, 4)
            else:
                size = (10, 4)
        gtype = '{}_gene'.format(gene.lower())
        short_chain = chain[0].upper()
        data = [s[gtype][l].rstrip('D') for s in seqs if gtype in s]
        x, y = aggregate(data, keys=keys)
        x = [i.replace('IG{}{}'.format(short_chain, gene), '{}{}'.format(gene, short_chain)) for i in x]
        plot_file = os.path.join(output_dir, '{}_{}{}_{}.pdf'.format(collection, gene, chain, l))
        colors = get_germline_plot_colors(x, l)
        make_barplot(x, y,
                     colors,
                     plot_file,
                     ylabel='Frequency (%)',
                     l=l,
                     size=size)


def get_germline_plot_colors(data, l):
    fams = [d.split('-')[0].replace('/OR', '') for d in data]
    nr_fams = sorted(list(set(fams)))
    num_colors = len(nr_fams)
    rgbs = sns.hls_palette(num_colors)
    rgb_dict = {i[0]: i[1] for i in zip(nr_fams, rgbs)}
    return [rgb_dict[f] for f in fams]



def isotypes():
    print('\n')
    print('Isotype plots not implemented ...')
    pass


def vj_heatmap(seqs, collection, make_plot, species, chain, output_dir):
    if not make_plot:
        return None
    plot_file = os.path.join(output_dir, '{}_VJheatmap_{}.pdf'.format(collection, chain))
    vj_data = group_by_vj(seqs, species, chain)
    vj_df = pd.DataFrame(vj_data)
    make_heatmap(vj_df.transpose(), plot_file)


def group_by_vj(data, species, chain):
    # from abtools.utils.germlines import germlines
    vs = list(set([g.name.split('*')[0] for g in get_germlines(species, 'V', chain=chain)]))
    js = list(set([g.name.split('*')[0] for g in get_germlines(species, 'J', chain=chain)]))
    vj = {}
    for v in vs:
        vj[v] = {}
        for j in js:
            vj[v][j] = 0
    for d in data:
        v = d['v_gene']['gene']
        j = d['j_gene']['gene'].rstrip('P')
        if v[:3] != j[:3]:
            continue
        if v not in vj:
            continue
        if j not in vj[v]:
            continue
        vj[v][j] = vj[v][j] + 1 if j in vj[v] else 1
    total = len(data)
    for v in list(vj.keys()):
        for j in list(vj[v].keys()):
            vj[v][j] = 100. * vj[v][j] / total
    return vj






def make_barplot(x, y, colors, ofile, xlabel=None, ylabel=None, l=None, grid=False, size=None, xfontsize=None):
    sns.set_style('ticks')
    # set bar locations and width
    ind = np.arange(len(x))
    width = 0.75
    # plot objects
    if size:
        fig = plt.figure(figsize=size)
    else:
        fig = plt.figure()
    ax = fig.add_subplot(111)
    # axis limits and ticks
    ax.set_ylim(0, 1.05 * max(y))
    ax.set_xlim(-width / 2, len(ind))
    ax.set_xticks(ind + width / 2)
    xtick_names = ax.set_xticklabels(x)
    if l == 'gene':
        plt.setp(xtick_names, rotation=90, fontsize=7)
    if grid:
        ax.yaxis.grid(True, alpha=0.5)
    # axis labels
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xfontsize:
        plt.setp(xtick_names, fontsize=xfontsize)
    ax.tick_params(axis='x', which='both', top='off', length=3, pad=1.5)
    ax.tick_params(axis='y', which='both', right='off', length=3, pad=1.5)
    # make the plot
    bar = ax.bar(ind, y, width, color=colors)
    fig.tight_layout()
    plt.savefig(ofile)
    plt.close()


def make_heatmap(df, ofile):
    sns.set()
    # set up plot, determine plot size
    h, w = df.shape
    f, ax = plt.subplots(figsize=(w / 1.75, h / 3))
    sns.heatmap(df,
                square=True,
                cmap='Blues',
                cbar=True,
                cbar_kws={'orientation': 'horizontal',
                          'fraction': 0.02,
                          'pad': 0.02,
                          'shrink': 0.675})
    # adjust labels
    ax.xaxis.tick_top()
    plt.xticks(rotation=90)
    f.tight_layout()
    # make the plot
    plt.savefig(ofile)
    plt.close()





def print_collection_info(collection):
    print('')
    print('')
    print('-' * 25)
    print(collection)
    print('-' * 25)
    print('')




def run(**kwargs):
    args = Args(**kwargs)
    main(args)


def run_standalone(args):
    logfile = args.log if args.log else os.path.join(args.output, 'abstats.log')
    log.setup_logging(logfile)
    global logger
    logger = log.get_logger('abstats')
    main(args)


def main(args):
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
        print_collection_info(group)
        seqs = get_seqs(group, args)
        name = group.split('.')[0]

        if len(seqs) == 0:
            print('No sequences to handle ...')
            continue

        if args.debug:
            print("Debugging ..., must show seq_id, vdj_nt, cdr3_len and v_gene if seqs are loaded:")
            print("Seq id :", seqs[0]['seq_id'])
            print("VDJ nt :", seqs[0]['vdj_nt'])
            print("CDR3 length :", seqs[0]['cdr3_len'])
            print("V gene :", seqs[0]['v_gene'])

        germline_plot(seqs, 'V', name, args.output, args.var_plot, args.species, args.chain)
        if args.chain == 'heavy':
            germline_plot(seqs, 'D', name, args.output, args.div_plot, args.species, args.chain)
        germline_plot(seqs, 'J', name, args.output, args.join_plot, args.species, args.chain)
        cdr3_plot(seqs, name, args.cdr3_plot, args.chain, args.output)
        vj_heatmap(seqs, name, args.heatmap, args.species, args.chain, args.output)
    print('\n')
    print('AbStats finished, successfully!.')



if __name__ == '__main__':
    args = parse_args()
    main()
