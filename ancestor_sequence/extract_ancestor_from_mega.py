#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extract ancestor sequences from mega text output
"""
import os
import sys
from collections import OrderedDict

with open(sys.argv[-1]) as mega_f:
    lines = mega_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line for line in lines if line]
    seqs = OrderedDict()
    for line in lines:
        name,seq = line.split(':')
        seq = seq.strip(' ')
        if '.' in name:
            if not name in seqs.keys():
                seqs[name] = seq
            else:
                seqs[name] = seqs[name]+seq

    with open('mega_ancestor.fa','w') as w_f:
        for name,seq in seqs.iteritems():
            print >> w_f, '{0:<40}{1}'.format(name,seq)
        for name,seq in seqs.iteritems():
            print >> w_f,'>{0}'.format(name)
            for s in [seq[i:i+80] for i in range(0,len(seq),80)]:
                print >> w_f,s


