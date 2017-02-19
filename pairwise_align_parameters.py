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
from Bio.SubsMat.MatrixInfo import *
from multiprocessing import Pool
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm

def align(seq1, seq2, matrix,matrix_name):
    # matrix = matlist.blosum62
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
    best_score = 0
    best_score_c = ''
    best_identity = 0
    best_identity_c = ''
    para = []
    for gap_open in opens:
        score_i = []
        identity_i = []
        for gap_extend in extends:
            if gap_open > gap_extend:
                score_i.append(0)
                identity_i.append(0)
                para.append([gap_open,gap_extend,0,0])
            else:
                alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
                seqq1 = alns[0][0]
                seqq2 = alns[0][1]
                identity = [1 for i, s in enumerate(seqq1) if s == seqq2[i]]
                identity = 1.0 * len(identity)/ len(seqq1)
                score_i.append(alns[0][2])
                identity_i.append(identity)
                if alns[0][2] > best_score:
                    best_score = alns[0][2]
                    best_score_c = (gap_open,gap_extend)
                if identity > best_identity:
                    best_identity = identity
                    best_identity_c = (gap_open,gap_extend)
                para.append([gap_open,gap_extend,alns[0][2],identity])

        scores.append(score_i)
        identities.append(identity_i)

    # plot scores
    # X_OFFSET,Y_OFFSET,Z_OFFSET = 0.5,0.1,20

    z_max = max([max([si for si in s ]) for s in scores])
    z_min = min([min([si for si in s ]) for s in scores])
    # if z_max > 0:
        # fig = plt.figure()
        # ax = fig.add_subplot(111,projection='3d')
        # ax.plot_surface(gx,gy,scores,alpha=0.3)
        # ax.contour(gx,gy,scores,zdir='z',offset=z_min-Z_OFFSET,cmap=cm.coolwarm)
        # ax.contour(gx,gy,scores,zdir='x',offset=OPEN_BEGIN-X_OFFSET,cmap=cm.coolwarm)
        # ax.contour(gx,gy,scores,zdir='y',offset=EXTEND_END+Y_OFFSET,cmap=cm.coolwarm)
        # ax.set_xlim(OPEN_BEGIN-X_OFFSET,OPEN_END+X_OFFSET)
        # ax.set_ylim(EXTEND_BEGIN-Y_OFFSET,EXTEND_END+Y_OFFSET)
        # ax.set_zlim(z_min-Z_OFFSET,z_max+Z_OFFSET)
        # # ax.plot_wireframe(gx,gy,scores)
        # ax.set_xlabel('gap_open',labelpad=10)
        # ax.set_ylabel('gap_extend',labelpad=10)
        # ax.set_zlabel('score',labelpad=10)
        # ax.set_title('Pairwise Alignment Score')
        # fig.savefig(str(matrix_name)+'_align_score.png')
        # plt.close()

    # # plot identities
    # X_OFFSET,Y_OFFSET,Z_OFFSET = 0.5,0.1,0.02

    # z_max = max([max([si for si in s ]) for s in identities])
    # z_min = min([min([si for si in s ]) for s in identities])
    # if z_max > 0:
        # fig = plt.figure()
        # ax = fig.add_subplot(111,projection='3d')
        # ax.plot_surface(gx,gy,identities,alpha=0.3)
        # ax.contour(gx,gy,identities,zdir='z',offset=z_min-Z_OFFSET,cmap=cm.coolwarm)
        # ax.contour(gx,gy,identities,zdir='x',offset=OPEN_BEGIN-X_OFFSET,cmap=cm.coolwarm)
        # ax.contour(gx,gy,identities,zdir='y',offset=EXTEND_END+Y_OFFSET,cmap=cm.coolwarm)
        # ax.set_xlim(OPEN_BEGIN-X_OFFSET,OPEN_END+X_OFFSET)
        # ax.set_ylim(EXTEND_BEGIN-Y_OFFSET,EXTEND_END+Y_OFFSET)
        # ax.set_zlim(z_min-Z_OFFSET,z_max+Z_OFFSET)
        # # ax.plot_wireframe(gx,gy,scores)
        # ax.set_xlabel('gap_open',labelpad=10)
        # ax.set_ylabel('gap_extend',labelpad=10)
        # ax.set_zlabel('identity',labelpad=10)
        # ax.set_title('Pairwise Alignment Identity')
        # fig.savefig(str(matrix_name)+'_align_identity.png')
        # plt.close()

    # if z_max > 0:
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # ax.scatter(scores,identities)
        # ax.set_xlabel('Score',labelpad=10)
        # ax.set_ylabel('Identity',labelpad=10)
        # ax.set_title('Pairwise Alignment Identity and Score')
        # fig.savefig(str(matrix_name)+'_align_score_identity.png')
        # plt.close()

    return [best_score,best_score_c,best_identity,best_identity_c,para]

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

def single_align(p):
    seq1,seq2,m = p
    matrix = eval(m)
    matrix_name = str(m)
    try:
        result = align(seq1,seq2,matrix,matrix_name) + [matrix_name]
        return result
    except:
        return ''
    # pa.append(result)

def main():
    with open(sys.argv[-1]) as fa_f:
        seqs = readfa(fa_f)
        seqlen = len(seqs)
        for i in range(0,len(seqs),2):
            print seqs[i][1]
            print seqs[i+1][1]
            pa = []
            # for m in matlist.available_matrices:
            matrix_num = len(matlist.available_matrices)
            parameters = [[seqs[i][1],seqs[i+1][1],m] for m in matlist.available_matrices]
            p = Pool(16)
            pa = p.map(single_align,parameters)
            p.close()

            # for m in matlist.available_matrices:
                # matrix = eval(m)
                # matrix_name = str(m)
                # print i
                # result = align(seqs[i][1],seqs[i+1][1],matrix,matrix_name) + [matrix_name]
                # pa.append(result)
            with open(str(i)+'_pa.txt','w') as w_f:
                print >> w_f,'{0:<12}{1:<12}{2:<15}{3:<15}{4:<12}{5:<12}{6:<}'.format('best_score','gap_open','gap_extend','best_identity','gap_open','gap_extend','matrix')
                for p in pa:
                    if p:
                        bs,bsc,bi,bic,para,mn = p
                        print >> w_f,'{0:<12}{1:<12.2f}{2:<15.2f}{3:<15.2f}{4:<12.2f}{5:<12.2f}{6:<}'.format(bs,bsc[0],bsc[1],bi,bic[0],bic[1],mn)
                        with open(str(i)+'_'+mn+'_para.txt','w') as para_f:
                            print >> para_f,'{0:<12}{1:<12}{2:<12}{3:<12}'.format('gap_open','gap_extend','score','identity')
                            for go,ge,sc,ie in para:
                                print >> para_f,'{0:<12.2f}{1:<12.2f}{2:<12}{3:<12}'.format(go,ge,sc,ie)

if __name__ == "__main__":
    main()
