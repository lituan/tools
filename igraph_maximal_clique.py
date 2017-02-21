#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
for a group of sequences,calculate pairwise identity
use igraph to find a group of sequences with pairwise identity in a certain range
"""

import sys
import os
import igraph
from multipleprocessing import Pool
from Bio import pairwise2
from Bio.SubMat import MatrixInfo as matlist


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

def readfa(fa_f):
    with open(fa_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        begin = [i for i,line in enumerate(lines) if '>' in line]
        seqs = [lines[b:e] for b,e in zip(begin,begin[1:]+[len(lines)])]
        seqs = [(seq[0][1:],''.join(seq[1:])) for seq in seqs]
        return seqs

def get_pim(seqs):
    scores = []
    seqnum = len(seqs)
    for i in range(seqnum):
        score_i = []
        for j in range(seqnum):
            if j < i:
                score_i.append(scores[j][i])
            elif j > i:
                score_i.append(align(seqs[i][1],seqs[j][1]))
            else:
                score_i.append(1.0)
        scores.append(score_i)
    return scores

def igraph_mc(labels,scores,cutoff1=0.6,cutoff2=0.9):

    adj_m = [map(lambda x: 1 if cutoff1 < x < cutoff2 else 0,row) for row in scores]
    for i in range(len(adj_m)):
        adj_m[i][i] = 0
    graph = igraph.Graph.Adjacency(adj_m,mode='undirected')
    indexes = igraph.maximal_cliques()
    mc_labels = [labels[i] for i in indexes[-1]]
    return mc_labels

def main():
    cutoff1 = 0.5
    cutoff2 = 0.9
    seqs = readfa(sys.argv[-1])
    scores = get_pim(seqs)
    mc_labels = igraph_mc([seq[0] for seq in seqs],scores,cutoff1,cutoff2)
    mc_seqs = [(pro,seq) for pro,seq in seqs if pro in mc_labels]

    fname = os.path.split(sys.argv[-1])[1].split('.')[0]
    fname = fname+'_mc_seqs_'+str(cutoff1)+'_'+str(cutoff2)
    with open(fname+'.fa','w') as w_f:
        for pro,seq in mc_seqs:
            print >> w_f,'>{0}'.format(pro)
            seqs = [seq[i:i+80] for i in range(0,len(seq),80)]
            for s in seqs:
                print >> w_f,s

if __name__ == "__main__":
    main()

