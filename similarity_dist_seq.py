#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot similarity distribution of sequences
"""

import sys
import os
import numpy as np
import cPickle as pickle
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns


def read_fa(fa_f):
    # readin seqs in fasta format
    # seqs foramt:[(pro,seq),...]
    with open(fa_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        pro_line_num = [i for i, line in enumerate(
            lines) if '>' in line] + [len(lines)]
        seqs = [lines[n:pro_line_num[i + 1]]
                for i, n in enumerate(pro_line_num[:-1])]
        seqs = [(seq[0][1:], ''.join(seq[1:])) for seq in seqs]
        return seqs


def align(p):
    i1,i2,seq1,seq2,blosum = p
    gap_open = -10  # usual value
    gap_extend = -0.5  # usual value
    alns = pairwise2.align.globalds(seq1, seq2, blosum, gap_open, gap_extend)
    seq1 = alns[0][0]
    seq2 = alns[0][1]
    identical_res = [1 for i, s in enumerate(seq1) if s == seq2[i]]
    identity = 1.0 * len(identical_res)/ len(seq1)
    return i1,i2,round(identity,4)


def get_similarity(seqs):

    blosums = [matlist.blosum30,matlist.blosum35,matlist.blosum40,matlist.blosum45, \
               matlist.blosum50,matlist.blosum55,matlist.blosum60,matlist.blosum62, \
               matlist.blosum65,matlist.blosum70,matlist.blosum75,matlist.blosum80, \
               matlist.blosum85,matlist.blosum90,matlist.blosum95,matlist.blosum100]

    seq_pairs = []
    seq_num = len(seqs)
    for i in range(seq_num):
        for j in range(seq_num):
            if j > i:
                for blosum in blosums:
                    seq_pairs.append((i,j,seqs[i][1],seqs[j][1],blosum))
    p = Pool(6)
    results = p.map(align,seq_pairs)
    p.close()

    scores = np.ones(shape=(seq_num,seq_num))
    for i,j,s in results:
        if scores[i][j] == 1:
            scores[i][j] = s
        elif scores[i][j] < s:
            scores[i][j] = s
        scores[j][i] = scores[i][j]
    pair_scores = []
    for i in range(seq_num):
        for j in range(seq_num):
            if j > i:
                pair_scores.append(scores[i][j])

    return pair_scores

def main():

    seqs = read_fa(sys.argv[-1])
    scores = get_similarity(seqs)

    f,ax = plt.subplots()
    sns.distplot(scores,hist=False)
    fname = os.path.split(sys.argv[-1])[1].split('.')[0]
    plt.savefig(fname+'_seq_similarity_dist_.png',dpi=300)

if __name__ == "__main__":
    main()


