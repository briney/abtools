#!/usr/bin/python
# filename: metadata.py


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



class Metadata():
    '''
    Docstring for Metadata.
    '''
    def __init__(self, metadata):
        self.raw = metadata if metadata is not None else {}
        self._epitope = None
        self._structures = None
        self._donor = None
        self._lineage = None
        self._isolation_paper = None
        self._neutralizing_antibody_feature_names = None
        self._germline_paper = None

    
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
