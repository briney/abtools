#!/usr/bin/python
# filename: bnabs.py


#
# Copyright (c) 2018 Bryan Briney
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


from __future__ import print_function, absolute_import, division

import csv
from io import StringIO
import json
import os
import sys

import pandas as pd

from abutils.core.pair import assign_pairs
from abutils.core.sequence import Sequence, read_json

from . import CATNAP_PATH
from .metadata import Metadata
from .neut import get_neutralization
from .paper import Paper

if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
else:
    STR_TYPES = [str, unicode]


# CATNAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CATNAP_data')


# class bnAbPaper(Paper):
#     '''
#     Docstring for bnAbPaper.
#     '''
#     def __init__(self, paper):
#         self.raw = paper
#         pmid = paper[paper.find('(')+1:paper.rfind(')')]
#         super.__init__(pmid)


# class bnAbMetadata():
#     '''
#     Docstring for bnAbMetadata.
#     '''
#     def __init__(self, metadata):
#         self.raw = metadata
#         self._epitope = None
#         self._structures = None
#         self._donor = None
#         self._lineage = None
#         self._isolation_paper = None
#         self._neutralizing_antibody_feature_names = None
#         self._germline_paper = None

    
#     @property
#     def donor(self):
#         return self.raw.get('Donor', None) 
    
#     @property
#     def epitope(self):
#         return self.raw.get('Antibody binding type', None)
    
#     @property
#     def lineage(self):
#         return self.raw.get('Clonal lineage', None)
    
#     @property
#     def structures(self):
#         if self._structures is None:
#             raw_structures = self.raw['PDB or other structure'].split(';')
#             self._structures = [bnAbStructure(rs) for rs in raw_structures]
#         return self._structures
    
#     @property
#     def isolation_paper(self):
#         if self._isolation_paper is None:
#             if 'Isolation paper(Pubmed ID)' in self.raw:
#                 paper_string = self.raw['Isolation paper(Pubmed ID)'].strip()
#                 self._isolation_paper = bnAbPaper(paper_string)
#         return self._isolation_paper
            
#     @property
#     def germline_paper(self):
#         if self._germline_paper is None:
#             if 'Germline paper(Pubmed ID)' in self.raw:
#                 paper_string = self.raw['Germline paper(Pubmed ID)'].strip()
#                 self._germline_paper = bnAbPaper(paper_string)
#         return self._germline_paper


class bnAbStructure():
    '''
    Docstring for bnAbStructure.
    '''
    def __init__(self, structure):
        self._raw = structure
        self._url = None
        self.pdb_id, self.description = self._parse_structure()

    def __repr__(self):
        r = ['PDB ID: {}'.format(self.pdb_id)]
        r.append(self.description)
        r.append(self.url)
        return '\n'.join(r)
    
    
    @property
    def raw(self):
        return self._raw
    
    @property
    def url(self):
        if self._url is None:
            self._url = 'https://www.rcsb.org/structure/{}'.format(self.pdb_id.lower())
        return self._url
    
    
    def _parse_structure(self):
        pdb_id = self.raw.split('(')[0]
        description = self.raw[self.raw.find('(')+1:self.raw.rfind(')')]
        return pdb_id, description


def get_bnab(name):
    bnabs = get_bnabs([name, ])
    if bnabs:
        return bnabs[0]
    return bnabs


def get_bnabs(names=None):
    if type(names) in STR_TYPES:
        names = [names, ]
    # read sequence data from JSON file
    bnab_file = os.path.join(CATNAP_PATH, 'ab_seqs.json')
    seqs = read_json(bnab_file)
    pairs = assign_pairs(seqs=seqs, delim='_')
    if names is not None:
        pairs = [p for p in pairs if p.name in names]
    # read the metadata file
    metadata = {}
    with open(os.path.join(CATNAP_PATH, 'ab_info.txt')) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['Name'].strip():
                metadata[row['Name']] = row
    # read the neut file
    # neutralization = {}
    # neut_file = os.path.join(CATNAP_PATH, 'neut.txt')
    # neutralization = pd.read_csv(neut_file,
    #                              sep='\t',
    #                              dtype=str,
    #                              index_col='Virus name')
    # neutralization = neutralization.T
    # for col in neutralization.columns.values:
    #     if col['Name'].strip():
    #         neutralization[col['Name']] = col
    # attach metadata and neut info to each pair
    for p in pairs:
        p.metadata = Metadata(metadata.get(p.name, None))
        p.neut = get_neutralization(antibody=p.name)
    return pairs


def get_bnab_dict():
    bnabs = get_bnabs()
    return {bnab.name: bnab for bnab in bnabs}
