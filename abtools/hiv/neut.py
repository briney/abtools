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


from abc import abstractproperty
import csv
import itertools
import os
import pickle
import sys

import numpy as np
import pandas as pd

from abutils.utils.database import SQLiteDatabase
from abutils.utils.decorators import lazy_property

import abtools.hiv.bnabs
import abtools.hiv.virus

from . import CATNAP_PATH
from .paper import Paper

if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
else:
    STR_TYPES = [str, unicode]


class Neutralization():
    '''
    Docstring for Neutralization
    '''
    def __init__(self, neuts):
        self.neuts = neuts

    def __iter__(self):
        '''
        iter: Returns an iterator over the ``neuts`` attribute
        '''
        return iter(self.neuts)
    
    def __contains__(self, item):
        '''
        bool: Returns ``True`` if the provided item (either an Ab or virus name)
              is present in the ``Neutralization`` object.
        '''
        if any([item in self.antibodies, item in self.viruses]):
            return True
        return False
    
    def __getitem__(self, key):
        # if type(key) == slice:
        #     return self.sequence[key]
        if key in self.antibodies:
            return self._abdict[key]
        elif key in self.viruses:
            return self._virdict[key]
        return None


    @lazy_property
    def antibodies(self):
        return list(set([n.antibody for n in self.neuts]))

    @lazy_property
    def viruses(self):
        viruses = []
        for n in self.neuts:
            viruses.append(n.virus.name)
            viruses += n.virus.aliases
        return list(set(viruses))
    
    @lazy_property
    def _abdict(self):
        d = {}
        for a in self.antibodies:
            d[a] = {}
            for v in self.viruses:
                neuts = [n for n in self.neuts if all([n.antibody ==  a,
                                                       n.virus == v])]
                d[a][v] = neuts
        return d
    
    @lazy_property
    def _virdict(self):
        d = {}
        for v in self.viruses:
            d[v] = {}
            for a in self.antibodies:
                neuts = [n for n in self.neuts if all([n.antibody ==  a,
                                                       n.virus == v])]
                d[v][a] = neuts
        return d



class Neut():
    '''
    Docstring for Neut

    Contains neutralization data for a single Ab/virus pair
    '''
    def __init__(self, neut_data):
        self.raw = neut_data

    @property
    def antibody(self):
        return self.raw[0]['Antibody']

    @lazy_property
    def virus(self):
        return abtools.hiv.virus.get_virus(self.raw[0]['Virus'])

    @property
    def ic50_values(self):
        return [i.value for i in self.ic50s]

    @property
    def geomean_ic50(self):
        return self._geomean(self.ic50_values)

    @property
    def ic80_values(self):
        return [i.value for i in self.ic80s]

    @property
    def geomean_ic80(self):
        return self._geomean(self.ic80_values)

    @property
    def id50_values(self):
        return [i.value for i in self.id50s]

    @property
    def geomean_id50(self):
        return self._geomean(self.id50_values)

    @lazy_property
    def ic50s(self):
        ic50s = []
        for r in self.raw:
            if r['IC50'].strip():
                ic50s.append(IC50(r))
        return ic50s
    
    @lazy_property
    def ic80s(self):
        ic80s = []
        for r in self.raw:
            if r['IC80'].strip():
                ic80s.append(IC80(r))
        return ic80s
    
    @lazy_property
    def id50s(self):
        id50s = []
        for r in self.raw:
            if r['ID50'].strip():
                id50s.append(ID50(r))
        return id50s

    
    @staticmethod
    def _geomean(values):
        '''
        Returns the geometric mean of a list of values.

        Geomeans are calculated as described here: 
        https://www.hiv.lanl.gov/components/sequence/HIV/neutralization/help.html#geomean
        '''
        if all(['>' in str(v) for v in values]):
            vals = list(set(values))
            if len(vals) == 1:
                return vals[0]
            else:
                return vals
        else:
            vals = [v for v in values if '>' not in str(v)]
            vals = [float(v.replace('<', '') if type(v) in STR_TYPES else v for v in vals)]
            prod = np.prod(vals)
            return prod**(1 / len(vals))




class ICValue():
    '''
    Docstring for ICValue
    '''
    def __init__(self, ic_data):
        self.raw = ic_data
    
    @abstractproperty
    def value_key(self):
        return None
        
    @lazy_property
    def value(self):
        if '>' in self.raw[self.value_key]:
            return self.raw[self.value_key]
        elif '<' in self.raw[self.value_key]:
            return self.raw[self.value_key]
        else:
            return float(self.raw[self.value_key])
        
    @lazy_property
    def antibody(self):
        return abtools.hiv.bnabs.get_bnab(self.raw['Antibody'])
    
    @lazy_property
    def virus(self):
        return abtools.hiv.virus.get_virus(self.raw['Virus'])
    
    @lazy_property
    def pmid(self):
        return self.raw['Pubmed ID']
    
    @lazy_property
    def reference(self):
        return self.raw['Reference']
    
    @lazy_property
    def paper(self):
        return Paper(self.pmid)
    

class IC50(ICValue):
    '''
    Docstring for IC50.
    '''
    def __init__(self, ic_data):
        super.__init__(ic_data)
    
    @property
    def value_key(self):
        return 'IC50'
        

class IC80(ICValue):
    '''
    Docstring for IC80.
    '''
    def __init__(self, ic_data):
        super.__init__(ic_data)

    @property
    def value_key(self):
        return 'IC80'


class ID50(ICValue):
    '''
    Docstring for ID50.
    '''
    def __init__(self, ic_data):
        super.__init__(ic_data)

    @property
    def value_key(self):
        return 'ID50'


class NeutDB(SQLiteDatabase):
    '''

    '''
    def __init__(self, name=None, direc=None, in_memory=False, table_name=None):
        super(NeutDB, self).__init__(name=name, direc=direc,
                                     in_memory=in_memory, table_name=table_name)


    @property
    def structure(self):
        return [('antibody', 'text'), ('virus', 'text'), ('reference', 'text'), ('pmid', 'text'),
                ('ic50', 'text'), ('ic80', 'text'), ('id50', 'text'), ('raw', 'text')]
    
    
    @property
    def columns(self):
        return ['Antibody', 'Virus', 'Reference', 'Pubmed ID', 'IC50', 'IC80', 'ID50']
    
    
    def antibodies(self, virus=None):
        if virus is not None:
            if isinstance(virus, tuple(STR_TYPES)):
                virus = [virus, ]
            virus = [abtools.hiv.virus.get_standardized_name(v) for v in virus]
            query = '''SELECT {0}.antibody, {0}.virus
                       FROM {0}
                       WHERE {0}.virus in ({1})'''.format(self.table_name, ','.join('?' * len(virus)))
            results = self.cursor.execute(query, tuple(virus))
            return list(set(r[0] for r in results))
        else:
            query = 'SELECT DISTINCT {0}.antibody from {0}'.format(self.table_name)
            results = self.cursor.execute(query)
            return [r[0] for r in results]
        
    
    def viruses(self, antibody=None):
        if antibody is not None:
            if isinstance(antibody, tuple(STR_TYPES)):
                antibody = [antibody, ]
            query = '''SELECT {0}.virus, {0}.antibody
                       FROM {0}
                       WHERE {0}.antibody in ({1})'''.format(self.table_name, ','.join('?' * len(antibody)))
            results = self.cursor.execute(query, tuple(antibody))
            return list(set(r[0] for r in results))
        else:
            query = 'SELECT DISTINCT {0}.virus from {0}'.format(self.table_name)
            results = self.cursor.execute(query)
            return [r[0] for r in results]
        
    
    def insert_one(self, data):
        '''
        Inserts a single entry.

        Inputs:
          data: data to be inserted (as an iterable or dict)
        '''
        if isinstance(data, dict):
            data = [data[c] for c in self.columns]
        data = data[:-1] + [pickle.dumps(data[-1], protocol=0)]
        with self.connection as conn:
            conn.execute(self.insert_cmd, data)


    def insert_many(self, data):
        '''
        Inserts multiple entries.

        Inputs:
          data - a list/generator of iterables or dicts, each containing
                 data corresponding to a single entry.
        '''
        if isinstance(data[0], dict):
            keys = self.columns
            data = [[d[k] for k in keys] for d in data]
        data = [d[:-1] + [pickle.dumps(d[-1], protocol=0)] for d in data]
        with self.connection as conn:
            conn.executemany(self.insert_cmd, data)
    
    
    def find(self, antibody=None, virus=None):
        query = 'SELECT {0}.antibody, {0}.virus, {0}.raw FROM {0} '.format(self.table_name)
        query_args = []
        if any([antibody is not None, virus is not None]):
            query += 'WHERE '
        if antibody is not None:
            if isinstance(antibody, tuple(STR_TYPES)):
                antibody = [antibody, ]
            query += '{0}.antibody IN ({1})'.format(self.table_name, ','.join('?' * len(antibody)))
            query_args += antibody
        if all([antibody is not None, virus is not None]):
            query += 'AND '
        if virus is not None:
            if isinstance(virus, tuple(STR_TYPES)):
                virus = [virus, ]
            query += '{0}.virus IN ({1})'.format(self.table_name, ','.join('?' * len(virus)))
            query_args += virus
        results = self.cursor.execute(query, tuple(query_args))
        return [pickle.loads(r[-1]) for r in results]


def make_neut_db():
    db = NeutDB(name='neut', direc=CATNAP_PATH)
    df = pd.read_csv(os.path.join(CATNAP_PATH, 'neut.txt'), sep='\t').fillna('')
    data = []
    for row in df.iterrows():
        d = row[1].to_dict()
        db = NeutDB(name='neut', direc=CATNAP_PATH)
        data.append([d[c] if d[c] else None for c in db.columns] + [d])
    db.insert_many(data)
    db.index('Antibody')
    db.index('Virus')
    return db


def get_neut_db():
    db_path = os.path.join(CATNAP_PATH, 'neut')
    if not os.path.isfile(db_path):
        make_neut_db()
    return NeutDB(name='neut', direc=CATNAP_PATH)


def get_neutralization(antibodies=None, viruses=None):
    if type(antibodies) in STR_TYPES:
        antibodies = [antibodies, ]
    if type(viruses) in STR_TYPES:
        viruses = [viruses, ]
    viruses = [abtools.hiv.virus.get_standardized_name(v) for v in viruses]
    neut_db = get_neut_db()
    if antibodies is None:
        antibodies = neut_db.antibodies()
    neuts = []
    for a in antibodies:
        if viruses is None:
            viruses = neut_db.viruses(antibody=a)
        for v in viruses:
            d = neut_db.find(antibody=a, virus=v)
            neuts.append(Neut(d))
    return Neutralization(neuts)


    # # read the neut file
    # neut_df = pd.read_csv(os.path.join(CATNAP_PATH, 'neut.txt'), sep='\t')
    # ## Would pandas be faster here??

    # # with open(os.path.join(CATNAP_PATH, 'neut.txt')) as f:
    # #     reader = csv.DictReader(f, delimiter='\t')
    # #     data = [row for row in reader if row['Antibody'].strip()]
    # if antibodies is None:
    #     antibodies = neut_df['Antibody'].unique()
    # neuts = []
    # for a in antibodies:
    #     df = neut_df[neut_df['Antibody'] == a]
    #     if viruses is None:
    #         viruses = df['Virus'].unique()
    #     for v in viruses:
    #         d = df[df['Virus'] == v].to_dict(orient='records')
    #         neuts.append(Neut(d))
    # return Neutralization(neuts)

