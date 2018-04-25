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

from collections import Counter
import csv
from io import StringIO
import json
import os
import requests
import sys

from Bio import SeqIO, Medline

from abstar import run as run_abstar

from abutils.core.pair import Pair, assign_pairs
from abutils.core.sequence import Sequence

if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
else:
    STR_TYPES = [str, unicode]



CATNAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CATNAP_data')



class bnAbMetadata():
    '''
    Docstring for bnAbMetadata.
    '''
    def __init__(self, metadata):
        self._raw = metadata
        self._epitope = None
        self._structures = None
        self._donor = None
        self._lineage = None
        self._isolation_paper = None
        self._neutralizing_antibody_feature_names = None
        self._germline_paper = None
        
    
    @property 
    def raw(self):
        return self._raw
    
    @property
    def donor(self):
        return self.raw.get('Donor', None) 
    
    @property
    def epitope(self):
        return self.raw.get('Antibody binding type', None)
    
    @property
    def lineage(self):
        return self.raw.get('Clonal lineage', None)
    
    @property
    def structures(self):
        if self._structures is None:
            raw_structures = self.raw['PDB or other structure'].split(';')
            self._structures = [bnAbStructure(rs) for rs in raw_structures]
        return self._structures
    
    @property
    def isolation_paper(self):
        if self._isolation_paper is None:
            if 'Isolation paper(Pubmed ID)' in self.raw:
                paper_string = self.raw['Isolation paper(Pubmed ID)'].strip()
                self._isolation_paper = bnAbPaper(paper_string)
        return self._isolation_paper
            
    @property
    def germline_paper(self):
        if self._germline_paper is None:
            if 'Germline paper(Pubmed ID)' in self.raw:
                paper_string = self.raw['Germline paper(Pubmed ID)'].strip()
                self._germline_paper = bnAbPaper(paper_string)
        return self._germline_paper


class bnAbPaper():
    '''
    Docstring for bnAbPaper.
    '''
    def __init__(self, paper):
        self._raw = paper
        self._pmid = None
        self._pmcid = None
        self._doi = None
        self._identifiers = self._retrieve_identifiers(paper)
        self._medline = self._retrieve_medline()
    
    def __str__(self):
        return self.citation
    
    def __repr__(self):
        r = [self.title]
        r.append(self.authors[0] + ' <...> ' + self.authors[-1])
        r.append(self.citation)
        return '\n'.join(r)
        
    @property
    def raw(self):
        return self._raw
    
    @property
    def medline(self):
        return self._medline
    
    @property
    def identifiers(self):
        return self._identifiers
    
    @property
    def pmid(self):
        if self.identifiers is not None:
            return self.identifiers['records'][0]['pmid']
        return None
    
    @property
    def pmcid(self):
        if self.identifiers is not None:
            return self.identifiers['records'][0]['pmcid']
        return None
    
    @property
    def doi(self):
        if self.identifiers is not None:
            return self.identifiers['records'][0]['doi']
        return None
    
    @property
    def abstract(self):
        return self._search_medline('AB')
    
    @property
    def authors(self):
        return self._search_medline('AU')
    
    @property
    def citation(self):
        return self._search_medline('SO')
    
    @property
    def date(self):
        return self._search_medline('DP')
    
    @property
    def journal(self):
        return self._search_medline('JT')
    
    @property
    def mesh(self):
        return self._search_medline('MH')
    
    @property
    def title(self):
        return self._search_medline('TI')
    

    def _retrieve_identifiers(self, paper):
        convert_url = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={}&format=json'
        pmid = paper[paper.find('(')+1:paper.rfind(')')]
        r = requests.get(convert_url.format(pmid))
        if r.status_code == requests.codes.ok:
            return r.json()
        else:
            return None
    
    
    def _retrieve_medline(self):
        medline_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={}&rettype=medline'
        r = requests.get(medline_url.format(self.pmid))
        if r.status_code == requests.codes.ok:
            return Medline.read(StringIO(r.text))
        else:
            return None
    
    
    def _search_medline(self, term):
        if self.medline is not None:
            return self.medline.get(term, None)
        return None


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


def get_bnabs(names=None):
    if type(names) in STR_TYPES:
        names = [names, ]
    # read sequence data from JSON file
    with open(os.path.join(CATNAP_PATH, 'ab_seqs.json')) as f:
        seqs = [Sequence(s) for s in json.load(f)]
    pairs = assign_pairs(seqs=seqs, delim='_')
    if names is not None:
        pairs = [p for p in pairs if p.name in names]
    # read the metadata file
    metadata = {}
    with open(os.path.join(CATNAP_PATH, 'ab_info.txt')) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['Antibody'].strip():
                metadata[row['Antibody']] = row
    # attach metadata to each pair
    for p in pairs:
        p.metadata = bnAbMetadata(metadata.get(p.name, None))
    return pairs


def get_bnab_dict():
    bnabs = get_bnabs()
    return {bnab.name: bnab for bnab in bnabs}
