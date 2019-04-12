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
import itertools

import numpy as np
import pandas as pd

from abutils.utils.decorators import lazy_property

from . import CATNAP_PATH
from .bnab import get_bnab
from .paper import Paper
from .virus import get_virus

if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
else:
    STR_TYPES = [str, unicode]


class Neutralization():
    '''
    Docstring for Neutralization
    '''
    def __init__(neuts):
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
        if key in self.antibodies:
            return self._abdict[key]
        elif key in self.viruses:
            return self._virdict[v]
        return None


    @lazy_property
    def antibodies(self):
        return list(set([n.raw['Antibody'] for n in self.neuts]))

    @lazy_property
    def viruses(self):
        return list(set([n.raw['Virus'] for n in self.neuts]))

    @lazy_property
    def _ic50s(self):
        ic50s = []
        for r in self.raw:
            if r['IC']
    
    @lazy_property
    def _abdict(self):
        d = {}
        for a in self.antibodies:
            d[a] = {}
            for v in self.viruses:
                neuts = [n for n in self.neuts if all([n.raw['Antibody'] ==  a,
                                                       n.raw['Virus'] == v])]
                d[a][v] = neuts
        return d
    
    @lazy_property
    def _virdict(self):
        d = {}
        for v in self.viruses:
            d[v] = {}
            for a in self.antibodies:
                neuts = [n for n in self.neuts if all([n.raw['Antibody'] ==  a,
                                                       n.raw['Virus'] == v])]
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

    @property
    def virus(self):
        return self.raw[0]['Virus']

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
        return get_bnab(self.raw['Antibody'])
    
    @lazy_property
    def virus(self):
        return get_virus(self.raw['Virus'])
    
    @lazy_property
    def pmic(self):
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



def get_neutralization(antibody=None, virus=None):
    # read the neut file
    with open(os.path.join(CATNAP_PATH, 'ab_info.txt')) as f:
        reader = csv.DictReader(f, delimiter='\t')
        data = [row for row in reader if row['Antibody'].strip()]
    antibodies = list(set([d['Antibody'] for d in data]))
    viruses = list(set([d['Virus'] for d in data]))
    neuts = []
    for a, v in itertools.product(antibodies, viruses):
        n = [d for d in data if all([d['Antibody'] == a, d['Virus'] == v])]
        neuts.append(Neut(n))
    # filter
    if type(antibody) in STR_TYPES:
        antibody = [antibody, ]
    if type(virus) in STR_TYPES:
        virus = [virus, ]
    if antibody is not None:
        neuts = [n for n in neuts if n.raw['Antibody'] in antibody]
    if virus is not None:
        neuts = [n for n in neuts if n.raw['Virus'] in virus]
    return Neutralization(neuts)

