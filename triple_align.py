#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
align sequences as triples

contact, imlituan@gmail.com
"""
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2

gap_open = -10
gap_extend = -1.5

def gap_function(begin,length):
    if length == 0:
        return 0
    elif length == 1:
        return gap_open
    else:
        return gap_open + gap_extend*(length-1)

def wrap_match(blosum):
    def match_function(triple1,triple2):
        matrix = blosum
        score = []
        for a1,a2 in zip(triple1[0],triple2[0]):
            if matrix.has_key((a1,a2)):
                score.append(matrix[(a1,a2)])
            else:
                score.append(matrix[(a2,a1)])
        return sum(score)
    return match_function


def local_align_triple(p):
    blosums = [matlist.blosum30,matlist.blosum35,matlist.blosum40,matlist.blosum45, \
               matlist.blosum50,matlist.blosum55,matlist.blosum60,matlist.blosum62, \
               matlist.blosum65,matlist.blosum70,matlist.blosum75,matlist.blosum80, \
               matlist.blosum85,matlist.blosum90,matlist.blosum95,matlist.blosum100]

    i1,i2,seq1,seq2 = p
    seq1 = [[seqi] for seqi in seq1]
    seq2 = [[seqi] for seqi in seq2]

    max_identity = 0
    max_align = []
    for blosum in blosums:
        alns = pairwise2.align.localcc(seq1, seq2,wrap_match(blosum),gap_function,gap_function,gap_char=['---'])
        inner_identity = 0
        inner_align = []
        for aln in alns:
            align1 = [a if not isinstance(a,list) else a[0] for a in aln[0] ]
            align2 = [a if not isinstance(a,list) else a[0] for a in aln[1] ]
            align1 = ''.join(align1)
            align2 = ''.join(align2)
            identical_res = [1 for si,s in enumerate(align1) if s == align2[si]]
            identity = 1.0*len(identical_res)/len(align1)
            align1 = [align1[i:i+3] for i in range(0,len(align1),3) ]
            align2 = [align2[i:i+3] for i in range(0,len(align2),3) ]
            if identity > inner_identity:
                inner_identity = identity
                inner_align = [align1,align2]
        if inner_identity > max_identity:
            max_identity = inner_identity
            max_align = inner_align

    return i1,i2,max_identity,max_align[0],max_align[1]


p = [1,2,['ENK','AND','EHK','RND','ANC'],['STV','EHK','PSV','EHK','PSV','PWY','TYV']]
align_triple(p)
