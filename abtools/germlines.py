#!/usr/bin/env python
# filename: germlines.py


###########################################################################
#
# Copyright (c) 2015 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################


import os
import sys

from Bio import SeqIO

from abtools.sequence import Sequence


def germline_names(species, segment, chain='all', resolution='allele'):
    chain_prefixes = _get_chain_prefixes(chain)
    direc = os.path.dirname(os.path.abspath(__file__))
    gg_file = os.path.join(direc, 'utils/germline_genes/{}_{}.fasta'.format(species, segment))
    gene_names = [s.id for s in SeqIO.parse(open(gg_file, 'r'), 'fasta')]
    gene_names = [g for g in gene_names if g[:3] in chain_prefixes]
    if resolution == 'family':
        gene_names = [g.split('-')[0] for g in gene_names]
    if resolution == 'gene':
        gene_names = [g.split('*')[0] for g in gene_names]
    return gene_names


def get_germline(gene, species):
    gene = gene.upper()
    segment = gene[3]
    direc = os.path.dirname(os.path.abspath(__file__))
    gg_file = os.path.join(direc, 'utils/germline_genes/{}_{}.fasta'.format(species, segment))
    gene = [s for s in SeqIO.parse(open(gg_file, 'r'), 'fasta') if s.id == gene][0]
    return Sequence(gene)


def _get_chain_prefixes(chain):
    prefixes = {'heavy': ['IGH'],
                'kappa': ['IGK'],
                'lambda': ['IKL'],
                'light': ['IGK', 'IGL'],
                'all': ['IGH', 'IGK', 'IGL']}
    return prefixes[chain]


def germlines(species, segment, chain):
    return germs[species][segment][chain]


germs = {'human':
            {'V':
                {'heavy': [
                    'IGHV1-18*01', 'IGHV1-18*02', 'IGHV1-18*03', 'IGHV1-18*04', 'IGHV1-2*01', 'IGHV1-2*02', 'IGHV1-2*03', 'IGHV1-2*04', 'IGHV1-2*05', 'IGHV1-24*01', 'IGHV1-3*01', 'IGHV1-3*02', 'IGHV1-45*01', 'IGHV1-45*02', 'IGHV1-45*03', 'IGHV1-46*01', 'IGHV1-46*02', 'IGHV1-46*03', 'IGHV1-58*01', 'IGHV1-58*02', 'IGHV1-69*01', 'IGHV1-69*02', 'IGHV1-69*03', 'IGHV1-69*04', 'IGHV1-69*05', 'IGHV1-69*06', 'IGHV1-69*07', 'IGHV1-69*08', 'IGHV1-69*09', 'IGHV1-69*10', 'IGHV1-69*11', 'IGHV1-69*12', 'IGHV1-69*13', 'IGHV1-69*14', 'IGHV1-69-2*01', 'IGHV1-69-2*02', 'IGHV1-8*01', 'IGHV1-8*02',
                    'IGHV2-26*01', 'IGHV2-5*01', 'IGHV2-5*02', 'IGHV2-5*03', 'IGHV2-5*04', 'IGHV2-5*05', 'IGHV2-5*06', 'IGHV2-5*08', 'IGHV2-5*09', 'IGHV2-70*01', 'IGHV2-70*02', 'IGHV2-70*03', 'IGHV2-70*04', 'IGHV2-70*05', 'IGHV2-70*06', 'IGHV2-70*07', 'IGHV2-70*08', 'IGHV2-70*09', 'IGHV2-70*10', 'IGHV2-70*11', 'IGHV2-70*12', 'IGHV2-70*13',
                    'IGHV3-11*01', 'IGHV3-11*03', 'IGHV3-11*04', 'IGHV3-11*05', 'IGHV3-11*06', 'IGHV3-13*01', 'IGHV3-13*02', 'IGHV3-13*03', 'IGHV3-13*04', 'IGHV3-13*05', 'IGHV3-15*01', 'IGHV3-15*02', 'IGHV3-15*03', 'IGHV3-15*04', 'IGHV3-15*05', 'IGHV3-15*06', 'IGHV3-15*07', 'IGHV3-15*08', 'IGHV3-20*01', 'IGHV3-20*02', 'IGHV3-21*01', 'IGHV3-21*02', 'IGHV3-21*03', 'IGHV3-21*04', 'IGHV3-23*01', 'IGHV3-23*02', 'IGHV3-23*03', 'IGHV3-23*04', 'IGHV3-23*05', 'IGHV3-30*01', 'IGHV3-30*02', 'IGHV3-30*03', 'IGHV3-30*04', 'IGHV3-30*05', 'IGHV3-30*06', 'IGHV3-30*07', 'IGHV3-30*08', 'IGHV3-30*09', 'IGHV3-30*10', 'IGHV3-30*11', 'IGHV3-30*12', 'IGHV3-30*13', 'IGHV3-30*14', 'IGHV3-30*15', 'IGHV3-30*16', 'IGHV3-30*17', 'IGHV3-30*18', 'IGHV3-30*19', 'IGHV3-30-3*01', 'IGHV3-30-3*02', 'IGHV3-30-3*03', 'IGHV3-33*01', 'IGHV3-33*02', 'IGHV3-33*03', 'IGHV3-33*04', 'IGHV3-33*05', 'IGHV3-33*06', 'IGHV3-43*01', 'IGHV3-43*02', 'IGHV3-48*01', 'IGHV3-48*02', 'IGHV3-48*03', 'IGHV3-48*04', 'IGHV3-49*01', 'IGHV3-49*02', 'IGHV3-49*03', 'IGHV3-49*04', 'IGHV3-49*05', 'IGHV3-53*01', 'IGHV3-53*02', 'IGHV3-53*03', 'IGHV3-53*04', 'IGHV3-64*01', 'IGHV3-64*02', 'IGHV3-64*03', 'IGHV3-64*04', 'IGHV3-64*05', 'IGHV3-66*01', 'IGHV3-66*02', 'IGHV3-66*03', 'IGHV3-66*04', 'IGHV3-7*01', 'IGHV3-7*02', 'IGHV3-7*03', 'IGHV3-72*01', 'IGHV3-72*02', 'IGHV3-73*01', 'IGHV3-73*02', 'IGHV3-74*01', 'IGHV3-74*02', 'IGHV3-74*03', 'IGHV3-9*01', 'IGHV3-9*02', 'IGHV3-9*03', 'IGHV3-NL1*01',
                    'IGHV4-28*01', 'IGHV4-28*02', 'IGHV4-28*03', 'IGHV4-28*04', 'IGHV4-28*05', 'IGHV4-28*06', 'IGHV4-28*07', 'IGHV4-30-2*01', 'IGHV4-30-2*02', 'IGHV4-30-2*03', 'IGHV4-30-2*04', 'IGHV4-30-2*05', 'IGHV4-30-2*06', 'IGHV4-30-4*01', 'IGHV4-30-4*02', 'IGHV4-30-4*03', 'IGHV4-30-4*04', 'IGHV4-30-4*05', 'IGHV4-30-4*06', 'IGHV4-30-4*07', 'IGHV4-31*01', 'IGHV4-31*02', 'IGHV4-31*03', 'IGHV4-31*04', 'IGHV4-31*05', 'IGHV4-31*06', 'IGHV4-31*07', 'IGHV4-31*08', 'IGHV4-31*09', 'IGHV4-31*10', 'IGHV4-34*01', 'IGHV4-34*02', 'IGHV4-34*03', 'IGHV4-34*04', 'IGHV4-34*05', 'IGHV4-34*06', 'IGHV4-34*07', 'IGHV4-34*08', 'IGHV4-34*09', 'IGHV4-34*10', 'IGHV4-34*11', 'IGHV4-34*12', 'IGHV4-34*13', 'IGHV4-38-2*01', 'IGHV4-38-2*02', 'IGHV4-39*01', 'IGHV4-39*02', 'IGHV4-39*03', 'IGHV4-39*04', 'IGHV4-39*05', 'IGHV4-39*06', 'IGHV4-39*07', 'IGHV4-4*01', 'IGHV4-4*02', 'IGHV4-4*03', 'IGHV4-4*04', 'IGHV4-4*05', 'IGHV4-4*06', 'IGHV4-4*07', 'IGHV4-4*08', 'IGHV4-59*01', 'IGHV4-59*02', 'IGHV4-59*03', 'IGHV4-59*04', 'IGHV4-59*05', 'IGHV4-59*06', 'IGHV4-59*07', 'IGHV4-59*08', 'IGHV4-59*09', 'IGHV4-59*10', 'IGHV4-61*01', 'IGHV4-61*02', 'IGHV4-61*03', 'IGHV4-61*04', 'IGHV4-61*05', 'IGHV4-61*06', 'IGHV4-61*07', 'IGHV4-61*08',
                    'IGHV5-10-1*01', 'IGHV5-10-1*02', 'IGHV5-10-1*03', 'IGHV5-10-1*04', 'IGHV5-51*01', 'IGHV5-51*02', 'IGHV5-51*03', 'IGHV5-51*04', 'IGHV5-51*05',
                    'IGHV6-1*01', 'IGHV6-1*02',
                    'IGHV7-4-1*01', 'IGHV7-4-1*02', 'IGHV7-4-1*03', 'IGHV7-4-1*04', 'IGHV7-4-1*05'
                    ],

                'kappa': ['IGKV1-12*01', 'IGKV1-12*02', 'IGKV1-13*01', 'IGKV1-13*02', 'IGKV1-16*01', 'IGKV1-16*02', 'IGKV1-17*01', 'IGKV1-17*02', 'IGKV1-17*03', 'IGKV1-27*01', 'IGKV1-33*01', 'IGKV1-39*01', 'IGKV1-39*02', 'IGKV1-5*01', 'IGKV1-5*02', 'IGKV1-5*03', 'IGKV1-6*01', 'IGKV1-6*02', 'IGKV1-8*01', 'IGKV1-9*01', 'IGKV1-NL1*01', 'IGKV1-12*01', 'IGKV1-12*02', 'IGKV1-13*01', 'IGKV1-13*02', 'IGKV1-16*01', 'IGKV1-16*02', 'IGKV1-17*01', 'IGKV1-33*01', 'IGKV1-39*01', 'IGKV1-43*01', 'IGKV1-8*01', 'IGKV1-8*02', 'IGKV1-8*03',
                    'IGKV2-24*01', 'IGKV2-28*01', 'IGKV2-29*01', 'IGKV2-29*02', 'IGKV2-29*03', 'IGKV2-30*01', 'IGKV2-30*02', 'IGKV2-40*01', 'IGKV2-40*02', 'IGKV2-26*01', 'IGKV2-26*02', 'IGKV2-26*03', 'IGKV2-28*01', 'IGKV2-29*01', 'IGKV2-29*02', 'IGKV2-30*01', 'IGKV2-40*01',
                    'IGKV3-11*01', 'IGKV3-11*02', 'IGKV3-15*01', 'IGKV3-20*01', 'IGKV3-20*02', 'IGKV3-NL1*01', 'IGKV3-NL2*01', 'IGKV3-NL3*01', 'IGKV3-NL4*01', 'IGKV3-NL5*01', 'IGKV3-11*01', 'IGKV3-11*02', 'IGKV3-15*01', 'IGKV3-15*02', 'IGKV3-15*03', 'IGKV3-20*01', 'IGKV3-7*01',
                    'IGKV4-1*01',
                    'IGKV5-2*01',
                    'IGKV6-21*01', 'IGKV6-21*02', 'IGKV6-21*01', 'IGKV6-21*02'
                    ],

                'lambda': [
                    'IGLV1-36*01', 'IGLV1-40*01', 'IGLV1-40*02', 'IGLV1-40*03', 'IGLV1-44*01', 'IGLV1-47*01', 'IGLV1-47*02', 'IGLV1-51*01', 'IGLV1-51*02', 'IGLV10-54*01', 'IGLV10-54*02', 'IGLV10-54*03',
                    'IGLV2-11*01', 'IGLV2-11*02', 'IGLV2-11*03', 'IGLV2-14*01', 'IGLV2-14*02', 'IGLV2-14*03', 'IGLV2-14*04', 'IGLV2-18*01', 'IGLV2-18*02', 'IGLV2-18*03', 'IGLV2-18*04', 'IGLV2-23*01', 'IGLV2-23*02', 'IGLV2-23*03', 'IGLV2-8*01', 'IGLV2-8*02', 'IGLV2-8*03',
                    'IGLV3-1*01', 'IGLV3-10*01', 'IGLV3-10*02', 'IGLV3-12*01', 'IGLV3-12*02', 'IGLV3-16*01', 'IGLV3-19*01', 'IGLV3-21*01', 'IGLV3-21*02', 'IGLV3-21*03', 'IGLV3-22*01', 'IGLV3-25*01', 'IGLV3-25*02', 'IGLV3-25*03', 'IGLV3-27*01', 'IGLV3-9*01', 'IGLV3-9*02',
                    'IGLV4-3*01', 'IGLV4-60*01', 'IGLV4-60*02', 'IGLV4-60*03', 'IGLV4-69*01', 'IGLV4-69*02',
                    'IGLV5-37*01', 'IGLV5-39*01', 'IGLV5-39*02', 'IGLV5-45*01', 'IGLV5-45*02', 'IGLV5-45*03', 'IGLV5-45*04', 'IGLV5-52*01',
                    'IGLV6-57*01', 'IGLV6-57*02',
                    'IGLV7-43*01', 'IGLV7-46*01', 'IGLV7-46*02',
                    'IGLV8-61*01', 'IGLV8-61*02', 'IGLV8-61*03',
                    'IGLV9-49*01', 'IGLV9-49*02', 'IGLV9-49*03'
                    ]
                },

            'D':
                {'heavy': [
                    'IGHD1-1*01', 'IGHD1-20*01', 'IGHD1-26*01', 'IGHD1-7*01',
                    'IGHD2-15*01', 'IGHD2-2*01', 'IGHD2-2*02', 'IGHD2-2*03', 'IGHD2-21*01', 'IGHD2-21*02', 'IGHD2-8*01', 'IGHD2-8*02',
                    'IGHD3-10*01', 'IGHD3-10*02', 'IGHD3-16*01', 'IGHD3-16*02', 'IGHD3-22*01', 'IGHD3-3*01', 'IGHD3-3*02', 'IGHD3-9*01',
                    'IGHD4-17*01', 'IGHD4-4*01',
                    'IGHD5-12*01', 'IGHD5-18*01', 'IGHD5-5*01',
                    'IGHD6-13*01', 'IGHD6-19*01', 'IGHD6-25*01', 'IGHD6-6*01',
                    'IGHD7-27*01'
                    ]
                },

            'J':
                {'heavy': [
                    'IGHJ1*01',
                    'IGHJ2*01',
                    'IGHJ3*01', 'IGHJ3*02',
                    'IGHJ4*01', 'IGHJ4*02', 'IGHJ4*03',
                    'IGHJ5*01', 'IGHJ5*02',
                    'IGHJ6*01', 'IGHJ6*02', 'IGHJ6*03', 'IGHJ6*04',
                    ],

                'kappa': [
                    'IGKJ1*01',
                    'IGKJ2*01', 'IGKJ2*02', 'IGKJ2*03', 'IGKJ2*04',
                    'IGKJ3*01',
                    'IGKJ4*01', 'IGKJ4*02',
                    'IGKJ5*01',
                    ],

                'lambda': [
                    'IGLJ1*01',
                    'IGLJ2*01',
                    'IGLJ3*01', 'IGLJ3*02',
                    'IGLJ6*01',
                    'IGLJ7*01', 'IGLJ7*02'
                    ]
                }
            }
        }
