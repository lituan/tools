#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
try different combinations of gap_open and gap_extend to see the changes of score and identity
input seqs in fasta format, output sequence similarity matrix
usuage example
python pairwise_align_parameter example.fasta
"""
import os
import sys
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def align(seq1, seq2):
    matrix = matlist.blosum62
    opens = range(-2,-35,-1)
    extends = map(lambda x: x/100.0, range(-1,-200,-1))
    # gap_open = -10  # usual value
    # gap_extend = -0.5  # usual value

    scores = []
    identities = []
    for gap_open in opens:
        score_i = []
        identity_i = []
        for gap_extend in extends:
            alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
            seqq1 = alns[0][0]
            seqq2 = alns[0][1]
            identity = [1 for i, s in enumerate(seqq1) if s == seqq2[i]]
            identity = 1.0 * len(identity)/ len(seqq1)
            # grid_i.append((alns[0][2]))
            score_i.append(alns[0][2])
            identity_i.append(identity)

        scores.append(score_i)
        identities.append(identity_i)

    gx = [opens for i in range(len(extends))]
    gx = [[i for j in range(len(extends))] for i in opens]
    gy = [extends for i in range(len(opens)) ]
    # plot scores
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.plot_surface(gx,gy,scores,alpha=0.3)
    ax.contour(gx,gy,scores,zdir='z',offset=200,cmap=cm.coolwarm)
    ax.contour(gx,gy,scores,zdir='x',offset=-40,cmap=cm.coolwarm)
    ax.contour(gx,gy,scores,zdir='y',offset=0,cmap=cm.coolwarm)
    ax.set_xlim(-40,0)
    ax.set_ylim(-3,0)
    ax.set_zlim(200,900)
    # ax.plot_wireframe(gx,gy,scores)
    ax.set_xlabel('gap_open',labelpad=10)
    ax.set_ylabel('gap_extend',labelpad=10)
    ax.set_zlabel('score',labelpad=10)
    fig.savefig('align_score.png')

    # plot scores
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.plot_surface(gx,gy,identities,alpha=0.3)
    ax.contour(gx,gy,scores,zdir='z',offset=0.3,cmap=cm.coolwarm)
    ax.contour(gx,gy,scores,zdir='x',offset=-40,cmap=cm.coolwarm)
    ax.contour(gx,gy,scores,zdir='y',offset=0,cmap=cm.coolwarm)
    ax.set_xlim(-40,0)
    ax.set_ylim(-3,0)
    ax.set_zlim(0.3,0.38)
    ax.set_xlabel('gap_open',labelpad=10)
    ax.set_ylabel('gap_extend',labelpad=10)
    ax.set_zlabel('identity',labelpad=10)
    fig.savefig('align_identity.png')


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

def main():
    with open(sys.argv[-1]) as fa_f:
        seqs = readfa(fa_f)
        align(seqs[0][1],seqs[1][1])

if __name__ == "__main__":
    main()
