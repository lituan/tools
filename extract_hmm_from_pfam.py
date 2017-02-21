#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extract hmm from pfam-A.hmm according to pfam acc
usage: python PF00400 Pfam-A.hmm
"""

import sys
import os

acc = sys.argv[-2]
pfam_hmm = sys.argv[-1]
got = []
with open(pfam_hmm) as o_f:
    for line in o_f:
        line = line.rstrip('\r\n')
        if 'HMMER3' in line:
            hmm = []
            hmm.append(line)
        elif '//' in line:
            hmm.append(line)
            if acc in hmm[2]:
                got.append(hmm)
        else:
            hmm.append(line)


if len(got) > 0:
    with open(acc+'.hmm','w') as w_f:
        for h in got:
            for line in h:
                print >> w_f,line



