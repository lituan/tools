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


def read_msa(msa_f):
    with open(msa_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines if line]
        pro_line_num = [i for i, line in enumerate(
            lines) if '>' in line] + [len(lines)]
        seqs = [lines[n:pro_line_num[i + 1]]
                for i, n in enumerate(pro_line_num[:-1])]
        seqs = [(seq[0].split()[0][1:], ''.join(seq[1:])) for seq in seqs]

        return seqs

def get_similarity(seqs):

    def pim(seq1,seq2):
        identity = len([i for i,s in enumerate(seq1) if s == seq2[i]])
        return identity*1.0/len(seq1)

    scores = []
    seq_num = len(seqs)
    for i in range(seq_num):
        for j in range(seq_num):
            if j > i:
                scores.append(pim(seqs[i][1],seqs[j][1]))
    return scores

def main():
    seqs = read_msa(sys.argv[-1])
    scores = get_similarity(seqs)
    f,ax = plt.subplots()
    sns.distplot(scores)
    fname = os.path.split(sys.argv[-1])[1].split('.')[0]
    plt.savefig(fname+'_msa_similarity_dist_.png',dpi=300)

if __name__ == "__main__":
    main()


