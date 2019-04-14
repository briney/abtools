#!/usr/bin/python
# filename: neut.py


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

import os
import sys

import pandas as pd

from abutils.core.sequence import read_fasta
from abutils.utils.decorators import lazy_property

from . import CATNAP_PATH

if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
else:
    STR_TYPES = [str, unicode]


class Virus():
    '''
    Docstring for Virus
    '''
    def __init__(self, virus_data, nt_sequence=None, aa_sequence=None):
        self.raw = virus_data
        self.aligned_nt_sequence = nt_sequence
        self.aligned_aa_sequence = aa_sequence

    def __eq__(self, other):
        if other in [self.name] + self.aliases:
            return True
        return False
    
    @property
    def nt_sequence(self):
        if self.aligned_nt_sequence is not None:
            return self.aligned_nt_sequence.replace('-', '')
        return self.aligned_nt_sequence

    @property
    def aa_sequence(self):
        if self.aligned_aa_sequence is not None:
            return self.aligned_aa_sequence.replace('-', '')
        return self.aligned_aa_sequence
    
    @lazy_property
    def name(self):
        return self.raw['Virus name']
    
    @lazy_property
    def clade(self):
        return self.subtype
    
    @lazy_property
    def subtype(self):
        return self.raw['Subtype']
    
    @lazy_property
    def country(self):
        return self.raw['Country']
    
    @lazy_property
    def patient_health(self):
        return self.raw['Patient health']
    
    @lazy_property
    def days_post_infection(self):
        return self.raw['Days post infection']
    
    @lazy_property
    def days_from_seroconversion(self):
        return self.raw['Days from seroconversion']
    
    @lazy_property
    def fiebig(self):
        return self.raw['Fiebig']
    
    @lazy_property
    def risk_factor(self):
        return self.raw['Risk factor']
    
    @lazy_property
    def accession(self):
        return self.raw['Accession']
    
    @lazy_property
    def tier(self):
        return self.raw['Tier']
    
    @lazy_property
    def infection_stage(self):
        return self.raw['Infection stage']
    
    @lazy_property
    def aliases(self):
        return self.raw['Alias'].split(', ')
    
    @lazy_property
    def hiv_db_name(self):
        return self.raw['HIV DB name']
    
    @lazy_property
    def seq_data_exists(self):
        return self.raw['Seq data'] == 'Yes'
    
    @lazy_property
    def number_of_abs_tested(self):
        return self.raw['# of Abs tested']


def get_virus(virus):
    return get_viruses([virus, ])


def get_viruses(viruses=None):
    if type(viruses) in STR_TYPES:
        viruses = [viruses, ]
    # read the virus file
    virus_df = pd.read_csv(os.path.join(CATNAP_PATH, 'viruses.txt'), sep='\t')
    virus_aa = read_fasta(os.path.join(CATNAP_PATH, 'virus_seqs_aa.fasta'))
    virus_nt = read_fasta(os.path.join(CATNAP_PATH, 'virus_seqs_nt.fasta'))
    if viruses is None:
        viruses = virus_df['Virus name'].unique()
    virus_data = []
    for v in viruses:
        vdata = virus_df[virus_df['Virus name'] == v].to_dict()
        vaa = [s for s in virus_aa if s.id.split('.')[2] == v][0].sequence
        vnt = [s for s in virus_nt if s.id.split('.')[2] == v][0].sequence
        virus_data.append(Virus(vdata, nt_sequence=vnt, aa_sequence=vaa))
    if len(virus_data) == 1:
        return virus_data[0]
    return virus_data
