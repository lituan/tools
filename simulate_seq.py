#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
simulate a set of sequences
"""
import os
import sys
import numpy as np

AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
      'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

SEQ_LEN = 21
SEQ_NUM = 10
def random_seq(seqlen):
    seq = []
    for i in range(seqlen):
        seq.append(AA[np.random.randint(0,19)])
    return ''.join(seq)

with open('random_seq.fasta','w') as w_f:
    for i in range(SEQ_NUM):
        print >> w_f,'>seq'+str(i)
        seq = random_seq(SEQ_LEN)
        for s in [seq[i:i+80] for i in range(0,len(seq),80)]:
            print >> w_f,s

