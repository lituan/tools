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

def match_function(triple1,triple2):
    blosum = matlist.blosum62
    score = []
    for a1,a2 in zip(triple1[0],triple2[0]):
        if blosum.has_key((a1,a2)):
            score.append(blosum[(a1,a2)])
        else:
            score.append(blosum[(a2,a1)])
    return sum(score)

def align_triple(p):
    i1,i2,seq1,seq2 = p
    seq1 = [[seqi] for seqi in seq1]
    seq2 = [[seqi] for seqi in seq2]
    alns = pairwise2.align.localcc(seq1, seq2,match_function,gap_function,gap_function,gap_char=['---'])
    max_identity = 0
    max_align = []
    for aln in alns:
        align1 = [a if not isinstance(a,list) else a[0] for a in aln[0] ]
        align2 = [a if not isinstance(a,list) else a[0] for a in aln[1] ]
        align1 = ''.join(align1)
        align2 = ''.join(align2)
        identical_res = [1 for si,s in enumerate(align1) if s == align2[si]]
        identity = 1.0*len(identical_res)/len(align1)
        align1 = [align1[i:i+3] for i in range(0,len(align1),3) ]
        align2 = [align2[i:i+3] for i in range(0,len(align2),3) ]
        if identity > max_identity:
            max_identity = identity
            max_align = [align1,align2]

    print max_identity
    print max_align[0]
    print max_align[1]
    return i1,i2,max_identity,max_align

p = [1,2,['ENK','AND','EHK','RND','ANC'],['STV','EHK','PSV','EHK','PSV','PWY','TYV']]
align_triple(p)
