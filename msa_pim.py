#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
calculate pair-wise scores from multiple sequence alignment in clustal format
"""
import os
import sys
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


def read_msa(msa_f):
    with open(msa_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        msa = OrderedDict()
        for line in lines[1:]:
            words = line.split()
            if len(words) == 2:
                name,seq = words
                if not name in msa.keys():
                    msa[name] = seq
                else:
                    msa[name] = ''.join([msa[name],seq])

    return msa

def get_score(seq1,seq2):

    length = len([1 for i,j in zip(seq1,seq2) if (i!= '-' and j!= '-')])
    identity = [1 for i,j in zip(seq1,seq2) if (i==j and not '-' == i)]

    return '{0:4.2f}'.format(100.0*len(identity)/length)

def calculate_pim(msa):
    keys = msa.keys()
    length = len(keys)
    scores = []
    for i in xrange(length):
        score_i = []
        for j in xrange(length):
            if j < i:
                score_i.append(scores[j][i])
            elif j > i:
                score_i.append(get_score(msa[keys[i]],msa[keys[j]]))
            elif j == i:
                score_i.append(100.0)
        scores.append(score_i)
    return scores

def write_resutls(msa, scores, file_path, file_name):
    names = msa.keys()
    result = [[name] + score for name, score in zip(names, scores)]
    header = [['ID'] + names]
    result = header + result

    filename = os.path.join(file_path, file_name + '_scores_tab.txt')
    with open(filename, 'w') as w_f:
        for r in result:
            print >> w_f, '\t'.join([str(ri)for ri in r])

    result = align_lis_lis(result)
    filename = os.path.join(file_path, file_name + '_scores_align.txt')
    with open(filename, 'w') as w_f:
        for r in result:
            print >> w_f, '\t'.join([str(ri)for ri in r])

    pair_scores = []
    ids = names
    for i in range(len(ids)):
        for j in range(i):
            pair_scores.append([scores[i][j], (ids[i], ids[j])])

    pair_scores = sorted(pair_scores)
    filename = os.path.join(file_path, file_name + '_pair_scores.txt')
    with open(filename, 'w') as w_f:
        for score, pair in pair_scores:
            print >> w_f, score, '\t', '{0:<25}{1:<}'.format(pair[0], pair[1])


def main():
    msa = read_msa(sys.argv[-1])

    scores = calculate_pim(msa)

    file_path, file_name = os.path.split(sys.argv[-1])
    file_name, file_extention = os.path.splitext(file_name)
    write_resutls(msa, scores, file_path, file_name)

if __name__ == "__main__":
    main()








