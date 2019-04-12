#!/usr/bin/python
# filename: paper.py


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


from io import StringIO
import requests

from Bio import Medline


class Paper():
    '''
    Docstring for Paper.
    '''
    def __init__(self, pmid):
        self._pmid = pmid
        self._pmcid = None
        self._doi = None
        self._identifiers = self._retrieve_identifiers()
        self._medline = self._retrieve_medline()
    
    def __str__(self):
        return self.citation
    
    def __repr__(self):
        r = [self.title]
        r.append(self.authors[0] + ' <...> ' + self.authors[-1])
        r.append(self.citation)
        return '\n'.join(r)

    
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
    

    def _retrieve_identifiers(self):
        convert_url = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={}&format=json'
        r = requests.get(convert_url.format(self.pmid))
        if r.status_code == requests.codes.OK:
            return r.json()
        else:
            return None
    
    
    def _retrieve_medline(self):
        medline_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={}&rettype=medline'
        r = requests.get(medline_url.format(self.pmid))
        if r.status_code == requests.codes.OK:
            return Medline.read(StringIO(r.text))
        else:
            return None
    
    
    def _search_medline(self, term):
        if self.medline is not None:
            return self.medline.get(term, None)
        return None
