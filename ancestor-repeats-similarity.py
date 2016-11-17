#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
calculate similarity between repeats
"""
import sys
import os
from collections import OrderedDict


def align_lis_lis(lis_lis):
    """align and  nested list to print a table"""
    lis_lis = [[str(l) for l in lis]
               for lis in lis_lis]  # trans every element to str
    # make all inner lists of the same length
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [lis + (inner_lis_max_len - len(lis)) * [''] for lis in lis_lis]
    # trans list, so that the elements of the same column are in one list
    lis_lis = [[lis[i] for lis in lis_lis] for i in range(inner_lis_max_len)]
    # make element in the same list have the same length
    aligned = []
    for lis in lis_lis:
        width = max([len(l) for l in lis])
        lis = [l + (width - len(l)) * ' ' for l in lis]
        aligned.append(lis)
    # trans list_list to the original list_list
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [[lis[i] for lis in aligned] for i in range(inner_lis_max_len)]
    return lis_lis


def align(seq1, seq2):
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62
    gap_open = -10  # usual value
    gap_extend = -0.5  # usual value

    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)

    seq1 = alns[0][0]
    seq2 = alns[0][1]
    print seq1
    print seq2
    identity = [1 for i, s in enumerate(seq1) if s == seq2[i]]
    identity = 1.0 * len(identity) / len(seq1)
    return float('{0:<4.2f}'.format(identity))


def read_repeats(repeat_f):
    with open(repeat_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        pro_repeats = []
        for line in lines:
            words = line.split(':')
            pro = words[0].split('.')[0]
            repeats = words[1].split()
            repeats = [''.join(r.split('-')) for r in repeats]
            pro_repeats.append((pro, repeats))
        return pro_repeats


def get_similarity(pro_repeats):
    pro_repeats_similarity = []

    len_repeats = len(pro_repeats[0][1])
    title = ['id']
    for i in range(len_repeats):
        for j in range(len_repeats):
            if j > i:
                title.append(str(i) + '_' + str(j))
    pro_repeats_similarity.append(title)

    for pro, repeats in pro_repeats:
        similarity = []
        for i in range(len_repeats):
            for j in range(len_repeats):
                if j > i:
                    sim = align(repeats[i], repeats[j])
                    similarity.append(sim)
        pro_repeats_similarity.append([pro] + similarity)
    return pro_repeats_similarity


def main():
    repeat_f = sys.argv[-1]
    pro_repeats = read_repeats(repeat_f)
    pro_repeats_similarity = get_similarity(pro_repeats)

    trans_pro_repeats_similarity = [
        [lis[i] for lis in pro_repeats_similarity] for i in range(len(pro_repeats_similarity[0]))]
    trans_pro_repeats_similarity = align_lis_lis(trans_pro_repeats_similarity)
    with open('similarity_trans.txt', 'w') as w_f:
        for pro_similarity in trans_pro_repeats_similarity:
            print >> w_f, ' '.join(pro_similarity)

    pro_repeats_similarity = align_lis_lis(pro_repeats_similarity)
    with open('similarity.txt', 'w') as w_f:
        for pro_similarity in pro_repeats_similarity:
            print >> w_f, ' '.join(pro_similarity)

if __name__ == "__main__":
    main()
