#!/usr/bin/env python3

# bct and aed - 2018

'''
MIT License

Copyright (c) 2018 Banfield Lab - UC Berkeley

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

'''

import sys
import os
from collections import defaultdict
from itertools import chain
import argparse
import requests
from Bio import SeqIO


comp_tbl = str.maketrans('ACTGNactgnYRWSKMDVHBXyrwskmdvhbx',
                         'TGACNtgacnRYWSMKHBDVXrywsmkhbdvx')

bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]



def complement(s):
    return(s.translate(comp_tbl))

def extract_sequence(seq, s_start, s_end, comp):
    s_start -= 1
    subseq = seq[s_start:s_end]
    if comp == 'c':
        return(complement(subseq[::-1]))
    else:
        return(subseq)

def file_exists(f):
    return True if os.path.isfile(f) else False

def extract_sequence(seq, s_start, s_end, comp):
    s_start -= 1
    subseq = seq[s_start:s_end]
    if comp == 'c':
        return(complement(subseq[::-1]))
    else:
        return(subseq)

# count codons for sequence by window size (default whole sequence)
# use features to determine coding region
#
# param s str full length sequence
# param c dict dictionary of contig information
# param features list list of feature details
# param window_size int subseq size for counting
def count_codons(s, c, window_size=0):
    ws = len(s)
    #print(f'count_code: seq size is {ws}', file=sys.stderr)
    if window_size > 0:
        # window_size = window_size - window_size % 3  # make sure we count in 3's
        #print(f'window_size is {window_size}', file=sys.stderr)
        ws = window_size

    s = s.upper()  # uppercase sequence

    # compute all coding regions
    coding = dict()  # position -> codon
    for feat in c:
        seq = feat['seq']
        #print(seq)
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            coding[feat['start'] + i] = codon

    chunked = list()  # array of [from, to, chunk_counts_hash]
    for n in range(0, len(s), ws):
        #print(f'n={n} n+ws={n+ws-1}', file=sys.stderr)

        d = defaultdict(int)
        coding_length = 0
        thiswindow = {k: v for k, v in coding.items() if n <= k <= n+ws}
        for k, v in thiswindow.items():
            coding_length += 3
            d[v] += 1

        if (n+ws) > len(s):
            chunked.append([n, len(s)-1, coding_length, d])
        else:
            chunked.append([n, n+ws-1, coding_length, d])
    return(chunked)


# ctgname from to length TTT TTC ...
# % is calculated from number of codons in window divided by the number of bp in
# window in coding regions
def summarize_codons(chunks, ctg, name):
    hdr = 'contig\tstart\tend\tlength\t'
    for cdn in codons:
        hdr += f'{cdn}\t'
    hdr = hdr[:-1]
    print(hdr)

    for c in chunks:
        s = int(c[0])
        e = int(c[1])
        chunk_coding = c[2]
        d = c[3]
        l = e-s
        out = f'{name}\t{s}\t{e}\t{l}\t'

        for cdn in codons:
            if cdn in d:
                p = d[cdn]/float(chunk_coding)
            else:
                p = 0.0
            out += f'{p}\t'
        out = out[:-1]  # drop trailing tab
        print(out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate codon usage by window')
    parser.add_argument('-f', '--fasta', type=str,
                        help='path to fasta file of entire contig')
    parser.add_argument('-fna', '--features', type=str,
                        help='path to predicted protein file from contig')
    parser.add_argument('-w', '--window_size', type=int, default=0,
                        help='Window size (default: 0=use whole sequence')
    parser.add_argument('-v', '--verbose', action='store_true', default=False)

    args = parser.parse_args()
    fa = args.fasta
    print(fa)
    feats = args.features

    with open(fa) as f:
        for record in SeqIO.parse(f, 'fasta'):
            s = str(record.seq)
            name = str(record.id)

    c=[]
    with open(feats) as fea:
        for record in SeqIO.parse(fea, 'fasta'):
            rdes = str(record.description)
            start = rdes.split(' # ')[1]
            end = rdes.split(' # ')[2]
            #print(start)
            ide = str(record.id)
            atrs={}
            atrs['start']=int(start)
            atrs['end']=int(end)
            atrs['seq']=str(record.seq)
            atrs['name'] = ide

            c.append(atrs)

    chunks = count_codons(s, c, args.window_size)
    summarize_codons(chunks, c, name)
