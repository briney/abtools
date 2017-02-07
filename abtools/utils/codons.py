#!/usr/bin/python
# filename:  codons.py


#
# Copyright (c) 2015 Bryan Briney
# License:  The MIT license (http: //opensource.org/licenses/MIT)
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



codon_lookup = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
          "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
          "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
          "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
          "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
          "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
          "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
          "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
          "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
          "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
          "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
          "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
          "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
          "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
          "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
          "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", }
