#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
reform sequence so that each line contain 80 residues
and simplize uniprot-type id
"""

import os
import sys


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


def writefa(seqs):
    with open('seq80.fa', 'w') as w_f:
        seqs = [(pro.split()[0].split('|')[2], seq) for pro, seq in seqs]
        seqs = sorted(seqs)
        for pro, seq in seqs:
            print >> w_f, '>{}'.format(pro)
            for i in [seq[i:i + 80] for i in range(0, len(seq), 80)]:
                print >> w_f, i

with open(sys.argv[-1]) as fa_f:
    seqs = readfa(fa_f)
    writefa(seqs)
