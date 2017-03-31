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

def align(seq1, seq2):
    matrix = matlist.blosum62
    gap_open = -10  # usual value
    gap_extend = -0.5  # usual value

    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)

    seq1 = alns[0][0]
    seq2 = alns[0][1]
    identity = [1 for i, s in enumerate(seq1) if s == seq2[i]]
    identity = 1.0 * len(identity)/ len(seq1)

    return float('{0:<4.2f}'.format(identity))

def get_similarity(seqs):
    scores = []
    seq_num = len(seqs)
    for i in range(seq_num):
        for j in range(seq_num):
            if j > i:
                scores.append(align(seqs[i][1],seqs[j][1]))
    return scores

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
    seqs = read_fa('Classified1.fasta')
    scores1 = get_similarity(seqs)

    seqs = read_fa('Classified2.fasta')
    scores2 = get_similarity(seqs)

    seqs = read_fa('Classified3.fasta')
    scores3 = get_similarity(seqs)

    seqs = read_fa('Classified4.fasta')
    scores4 = get_similarity(seqs)

    pickle.dump([scores1,scores2,scores3,scores4],open('scores.pickle','w'))
    # scores1,scores2,scores3,scores4 = pickle.load(open('scores.pickle'))

    f,ax = plt.subplots()
    sns.distplot(scores1,hist=False)
    sns.distplot(scores2,hist=False)
    sns.distplot(scores3,hist=False)
    sns.distplot(scores4,hist=False)
    fname = os.path.split(sys.argv[-1])[1].split('.')[0]
    plt.savefig(fname+'_seq_similarity_dist_.png',dpi=300)

if __name__ == "__main__":
    main()


