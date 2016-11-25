#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
remove redundant sequences at certain threthold
method description:
1. calculate all pairwise-identity from msa
2. find neighbors for all sequences
3. find the one with most neighbors
4. remove redundant sequences for this cluster
5. repeat procedures 2-4 until no removement is needed
"""

import sys
import os
import lt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

@lt.run_time
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

@lt.run_time
def get_pim(seqs):

    def pim(seq1,seq2):
        identity = len([i for i,s in enumerate(seq1) if s == seq2[i]])
        return identity*1.0/len(seq1)

    scores = []
    seqlen = len(seqs)
    for i in range(seqlen):
        score_i = []
        for j in range(seqlen):
            if j < i:
                score_i.append(scores[j][i])
            elif j > i:
                score_i.append(pim(seqs[i][1],seqs[j][1]))
            else:
                score_i.append(1.0)
        scores.append(score_i)
    return scores

@lt.run_time
def get_neighbors(seqs,scores,cutoff):
    neighbors = []
    for i,seq in enumerate(seqs):
        seq_neighbors = [seq[0]]
        for j in range(len(seqs)):
            if j != i:
                if scores[i][j] > cutoff:
                    seq_neighbors.append(seqs[j][0])
        neighbors.append(seq_neighbors)
    return neighbors

@lt.run_time
def remove_redundancy(seqs,scores,neighbors,good_list):
    while max([len(n) for n in neighbors]) > 1:
        neighbors = sorted(neighbors,key=lambda x: len(x),reverse=True)
        for seq in neighbors[0][1:]:
            if not seq in good_list:
                for n in neighbors:
                    if seq in n:
                        n.pop(n.index(seq))
                neighbors = [n for n in neighbors if len(n) > 0]
                neighbors = [neighbors[0]] + [n for n in neighbors[1:] if not n[0] == seq]
        if len(neighbors[0]) > 1:
            for seq in neighbors[0][:-1]:
                for n in neighbors:
                    if seq in n:
                        n.pop(n.index(seq))
                neighbors = [n for n in neighbors if len(n) > 0]
                neighbors = [n for n in neighbors if not n[0] == seq]

    nr_list = [n[0] for n in neighbors]
    nr_seqs = [seq for seq in seqs if seq[0] in nr_list ]
    nr_index = [i for i,seq in enumerate(seqs) if seq[0] in nr_list]
    nr_scores = [[score_i[i] for i in nr_index] for score_i in scores]
    nr_scores = [nr_scores[i] for i in nr_index]
    return nr_seqs,nr_scores


@lt.run_time
def plot_heatmap(seqs, scores,file_name):

    column_labels = [s[0] for s in seqs]
    row_labels = column_labels
    scores = [map(lambda x: float(x), row) for row in scores]
    scores = np.array(scores)
    df = pd.DataFrame(scores,columns=column_labels, index=row_labels)

    if len(df.columns) > 50:
        # sns_plot = sns.clustermap(df,annot=True,fmt='3.2f',xticklabels='',yticklabels='')
        sns_plot = sns.clustermap(df,xticklabels='',yticklabels='')
    else:
        # sns_plot = sns.clustermap(df,annot=True,fmt='3.2f')
        sns_plot = sns.clustermap(df)
        plt.setp(sns_plot.ax_heatmap.yaxis.get_majorticklabels(), rotation=20)
        plt.setp(sns_plot.ax_heatmap.xaxis.get_majorticklabels(), rotation=70)
    # plt.yticks(rotation=90)
    sns_plot.savefig(file_name+'.png')

@lt.run_time
def write_msa(seqs,filename):
    with open(filename+'.fasta','w') as w_f:
        for pro,seq in seqs:
            print >> w_f, '>{0}'.format(pro)
            for s in [seq[i:i+80] for i in range(0,len(seq),80)]:
                print >> w_f,s



@lt.run_time
def main():
    cutoff = 0.9
    good_list = []
    seqs = read_msa(sys.argv[-1])
    scores = get_pim(seqs)
    neighbors = get_neighbors(seqs,scores,cutoff)
    nr_seqs,nr_scores = remove_redundancy(seqs,scores,neighbors,good_list)
    plot_heatmap(nr_seqs,nr_scores,'nr_seqs_'+str(cutoff))
    write_msa(nr_seqs,'nr_seqs_'+str(cutoff))

main()




