#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
input seqs in fasta format, output sequence similarity matrix
usuage example
python pairwise_align example.fasta
"""
import os
import sys
# from wdsp import Wdsp


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
    identity = [1 for i, s in enumerate(seq1) if s == seq2[i]]
    identity = 1.0 * len(identity)/ len(seq1)


    return float('{0:<4.2f}'.format(identity))


def readfa(fa_f):
    # readin seqs in fasta format
    # seqs foramt:[(pro,seq),...]
    lines = fa_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    pro_line_num = [i for i, line in enumerate(
        lines) if '>' in line] + [len(lines)]
    seqs = [lines[n:pro_line_num[i + 1]]
            for i, n in enumerate(pro_line_num[:-1])]
    seqs = [(seq[0][1:], ''.join(seq[1:])) for seq in seqs]
    return seqs


def read_wdsp(wdsp_f):
    w = Wdsp(wdsp_f)
    seqs = w.seqs
    seqs = [(k, v) for k, v in seqs.iteritems()]
    return seqs


def align_seqs(seqs):
    lens = len(seqs)
    scores = []
    for i in xrange(lens):
        score_i = []
        for j in xrange(lens):
            # print i, '\t', j
            if j < i:
                score_i.append(scores[j][i])
            elif j >= i:
                score = align(seqs[i][1], seqs[j][1])
                score_i.append(score)
        scores.append(score_i)
    return scores


def write_resutls(seqs, scores, file_path, file_name):
    result = [[seq[0]] + score for seq, score in zip(seqs, scores)]
    header = [['ID'] + [seq[0] for seq in seqs]]
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
    ids = [seq[0] for seq in seqs]
    for i in range(len(ids)):
        for j in range(i):
            pair_scores.append([scores[i][j], (ids[i], ids[j])])

    pair_scores = sorted(pair_scores)
    filename = os.path.join(file_path, file_name + '_pair_scores.txt')
    with open(filename, 'w') as w_f:
        for score, pair in pair_scores:
            print >> w_f, score, '\t', '{0:<25}{1:<}'.format(pair[0], pair[1])

    # change similarity to dissimilarity
    length = len(scores)
    distances = []
    for i in xrange(length):
        distance = []
        for j in xrange(length):
            distance.append(1.0 - scores[i][j])
        distances.append(distance)

    result = [[seq[0]] + score for seq, score in zip(seqs, distances)]
    header = [['ID'] + [seq[0] for seq in seqs]]
    result = header + result

    filename = os.path.join(file_path, file_name + '_distances_tab.txt')
    with open(filename, 'w') as w_f:
        for r in result:
            print >> w_f, '\t'.join([str(ri)for ri in r])

    result = align_lis_lis(result)
    filename = os.path.join(file_path, file_name + '_distances_align.txt')
    with open(filename, 'w') as w_f:
        for r in result:
            print >> w_f, '\t'.join([str(ri)for ri in r])


def plot_heatmap(seqs, scores,file_name):
    import matplotlib.pyplot as plt
    from numpy import array

    column_labels = [s[0] for s in seqs]
    row_labels = column_labels
    scores = array(scores)

    fig, ax = plt.subplots()
    ax.axis('off')
    heatmap = ax.pcolor(scores, cmap=plt.cm.Blues)
    cb = plt.colorbar(heatmap)
    ax.set_xticklabels(row_labels,minor=False)
    ax.set_yticklabels(column_labels,minor=False)
    fig.savefig(file_name)


def main():
    with open(sys.argv[-1]) as fa_f:
        seqs = readfa(fa_f)
    # with open(sys.argv[-1]) as wdsp_f:
        # seqs = read_wdsp(wdsp_f)

    scores = align_seqs(seqs)

    file_path, file_name = os.path.split(sys.argv[-1])
    file_name, file_extention = os.path.splitext(file_name)
    write_resutls(seqs, scores, file_path, file_name)
    plot_heatmap(seqs, scores,file_name)

if __name__ == "__main__":
    main()
