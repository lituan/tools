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
import numpy as np
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def align(seq1, seq2):
    matrix = matlist.blosum62
    # gap_open = -10  # usual value
    # gap_extend = -0.5  # usual value

    OPEN_BEGIN,OPEN_END,OPEN_STEP = -20.5,-2.0,0.1
    EXTEND_BEGIN,EXTEND_END,EXTEND_STEP = -2.0,-0.1,0.1

    opens = np.arange(OPEN_BEGIN,OPEN_END,OPEN_STEP)
    extends = np.arange(EXTEND_BEGIN,EXTEND_END,EXTEND_STEP)
    gx = [[opens[i] for j in range(len(extends))] for i in range(len(opens))]
    gy = [extends for i in range(len(opens))]
    scores = []
    identities = []
    for gap_open in opens:
        score_i = []
        identity_i = []
        for gap_extend in extends:
            if gap_open > gap_extend:
                score_i.append(0)
                identity_i.append(0)
            else:
                alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
                seqq1 = alns[0][0]
                seqq2 = alns[0][1]
                identity = [1 for i, s in enumerate(seqq1) if s == seqq2[i]]
                identity = 1.0 * len(identity)/ len(seqq1)
                score_i.append(alns[0][2])
                identity_i.append(identity)

        scores.append(score_i)
        identities.append(identity_i)

    # plot scores
    X_OFFSET,Y_OFFSET,Z_OFFSET = 0.5,0.1,20

    z_max = max([max([si for si in s if si != 0]) for s in scores])
    z_min = min([min([si for si in s if si != 0]) for s in scores])
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.plot_surface(gx,gy,scores,alpha=0.3)
    ax.contour(gx,gy,scores,zdir='z',offset=z_min-Z_OFFSET,cmap=cm.coolwarm)
    ax.contour(gx,gy,scores,zdir='x',offset=OPEN_BEGIN-X_OFFSET,cmap=cm.coolwarm)
    ax.contour(gx,gy,scores,zdir='y',offset=EXTEND_END+Y_OFFSET,cmap=cm.coolwarm)
    ax.set_xlim(OPEN_BEGIN-X_OFFSET,OPEN_END+X_OFFSET)
    ax.set_ylim(EXTEND_BEGIN-Y_OFFSET,EXTEND_END+Y_OFFSET)
    ax.set_zlim(z_min-Z_OFFSET,z_max+Z_OFFSET)
    # ax.plot_wireframe(gx,gy,scores)
    ax.set_xlabel('gap_open',labelpad=10)
    ax.set_ylabel('gap_extend',labelpad=10)
    ax.set_zlabel('score',labelpad=10)
    ax.set_title('Pairwise Alignment Score')
    fig.savefig('align_score.png')

    # plot identities
    X_OFFSET,Y_OFFSET,Z_OFFSET = 0.5,0.1,0.02

    z_max = max([max([si for si in s if si != 0]) for s in identities])
    z_min = min([min([si for si in s if si != 0]) for s in identities])
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.plot_surface(gx,gy,identities,alpha=0.3)
    ax.contour(gx,gy,identities,zdir='z',offset=z_min-Z_OFFSET,cmap=cm.coolwarm)
    ax.contour(gx,gy,identities,zdir='x',offset=OPEN_BEGIN-X_OFFSET,cmap=cm.coolwarm)
    ax.contour(gx,gy,identities,zdir='y',offset=EXTEND_END+Y_OFFSET,cmap=cm.coolwarm)
    ax.set_xlim(OPEN_BEGIN-X_OFFSET,OPEN_END+X_OFFSET)
    ax.set_ylim(EXTEND_BEGIN-Y_OFFSET,EXTEND_END+Y_OFFSET)
    ax.set_zlim(z_min-Z_OFFSET,z_max+Z_OFFSET)
    # ax.plot_wireframe(gx,gy,scores)
    ax.set_xlabel('gap_open',labelpad=10)
    ax.set_ylabel('gap_extend',labelpad=10)
    ax.set_zlabel('identity',labelpad=10)
    ax.set_title('Pairwise Alignment Identity')
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
