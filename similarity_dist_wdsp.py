
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot similarity distribution of sequences
"""

import sys
import os
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from multiprocessing import Pool


def read_wdsp(wdsp_f):
    # readin seqs in fasta format
    # seqs foramt:[(pro,seq),...]
    with open(wdsp_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line.split() for line in lines]
        lines = [line for line in lines if len(line) > 1]
        begin = [i for i,line in enumerate(lines) if '>' in line]
        end = begin[1:] + [len(lines)]
        entries = [lines[b:e] for b,e in zip(begin,end)]
        seqs = [(e[0],''.join([''.join(ei[3:-1]) for ei in e[1:]])) for e in entries]
        seqs = [(seq[0][1], ''.join(seq[1:])) for seq in seqs]
        print seqs
        return seqs

def align(p):
    s1,s2,seq1,seq2 = p
    matrix = matlist.blosum62
    gap_open = -10  # usual value
    gap_extend = -0.5  # usual value

    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)

    seq1 = alns[0][0]
    seq2 = alns[0][1]
    identity = [1 for i, s in enumerate(seq1) if s == seq2[i]]
    identity = 1.0 * len(identity)/ len(seq1)

    return s1,s2,float('{0:<4.2f}'.format(identity))

def get_similarity(seqs):
    seq_pairs = []
    seq_num = len(seqs)
    for i in range(seq_num):
        for j in range(seq_num):
            if j > i:
                seq_pairs.append((i,j,seqs[i][1],seqs[j][1]))
    p = Pool(6)
    results = p.map(align,seq_pairs)
    p.close()

    results = sorted(results)
    scores = np.ones(shape=(seq_num,seq_num))
    for i,j,s in results:
        scores[i][j] = s
        scores[j][i] = s

    return [r[2] for r in results]


def main():
    seqs = read_wdsp(sys.argv[-1])
    scores = get_similarity(seqs)
    f,ax = plt.subplots()
    sns.distplot(scores)
    fname = os.path.split(sys.argv[-1])[1].split('.')[0]
    plt.savefig(fname+'_wdsp_similarity_dist_.png',dpi=300)

if __name__ == "__main__":
    main()


