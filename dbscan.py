#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
use dbscan to cluster sequences accordint to identity, and
get non-redundant seqs accordint to a identity threthold
not finished yet, dbscan cannot be used
"""

import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sys
import numpy as np
import matplotlib as mpl
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as spd

IDENTITY_CUTOFF = 0.9

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


    return float('{0:<4.3f}'.format(identity))


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

def dbscan(scores,eps=0.1):
    distance = [map(lambda x: 1.0-x,row) for row in scores]
    from sklearn.cluster import DBSCAN
    db = DBSCAN(eps=eps,min_samples=1,metric='precomputed').fit(distance)
    return db.labels_

def cluster_pros(pros,scores,labels):
    labels = list(labels)
    n = len(set(labels))
    clusters = []
    if -1 in labels:
        independent = [[i] for i,l in enumerate(labels) if l == -1]
        clusters += independent
        for label_n in range(n-1):
            c = [i for i,l in enumerate(labels) if l == label_n]
            clusters.append(c)
    else:
        for label_n in range(n):
            c = [i for i,l in enumerate(labels) if l == label_n]
            clusters.append(c)

    # sort pros according to cluster
    new_pros = [[pros[ci] for ci in c] for c in clusters]
    # sort scores according to cluster
    clu = []
    for c in clusters:
        clu += c
    new_scores = [[score_i[ci] for ci in clu] for score_i in scores]
    new_scores = [new_scores[ci] for ci in clu]

    # get non_redundant pros and scores
    non_redundant_pros = [pros[c[0]]  for c in clusters]
    non_redundant_scores = [[score_i[c[0]] for c in clusters] for score_i in scores]
    non_redundant_scores = [non_redundant_scores[c[0]] for c in clusters]

    return new_pros, new_scores, non_redundant_pros, non_redundant_scores

def write_results(nr_seqs):
    with open('nr_seqs.fa','w') as w_f:
        for pro,seq in nr_seqs:
            print >> w_f, '>{0}'.format(pro)
            for s in [seq[i:i+80] for i in xrange(0,len(seq),80)]:
                print >> w_f,s

def BlueGreenYellow():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                        (0.25, 0.0, 0.0),
                        (0.50, 0.5, 0.5),
                        (1.0, 1.0, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (0.25, 0.75, 0.75),
                       (0.50, 1.0, 1.0),
                       (1.0, 1.0, 1.0)),

             'blue':  ((0.0, 1.0, 1.0),
                         (0.25, 1.0, 1.0),
                         (0.75, 0.0, 0.0),
                         (1.0, 0.0, 0.0))
            }
    ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
    ### modulate between blue and cyan using the last y var in the first green tuple
    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def plot_heatmap(seqs, scores,file_name,method='average'):

    column_labels = [s[0] for s in seqs]
    row_labels = column_labels
    scores = [map(lambda x: float(x), row) for row in scores]
    scores = np.array(scores)
    distances = [map(lambda x: 1-x,row) for row in scores]
    linkage = sch.linkage(spd.squareform(distances),method=method)
    df = pd.DataFrame(scores,columns=column_labels, index=row_labels)

    if len(df.columns) > 20:
        sns_plot = sns.clustermap(df,row_linkage=linkage,col_linkage=linkage,xticklabels='',yticklabels='')
    else:
        sns_plot = sns.clustermap(df,figsize=figsize,row_linkage=linkage,col_linkage=linkage,annot=True,fmt='3.2f')
        plt.setp(sns_plot.ax_heatmap.yaxis.get_majorticklabels(), rotation=20)
        plt.setp(sns_plot.ax_heatmap.xaxis.get_majorticklabels(), rotation=70)
    # plt.yticks(rotation=90)
    sns_plot.savefig(file_name+'.png')
    plt.close('all')

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
    fa_f = open(sys.argv[-1])
    seqs = readfa(fa_f)
    pros = [seq[0] for seq in seqs]
    scores = align_seqs(seqs)
    eps = 1.0 - IDENTITY_CUTOFF
    labels = dbscan(scores,eps)
    new_pros,new_scores,nr_pros,nr_scores = cluster_pros(pros,scores,labels)
    nr_seqs = [seq for seq in seqs if seq[0] in nr_pros]
    write_results(nr_seqs)
    plot_heatmap(nr_seqs,nr_scores,'nr_seqs_identity_heatmap')



main()
