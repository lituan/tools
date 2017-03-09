#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
select wdsps with 8 or more blades
"""
import sys
import os
from collections import OrderedDict

with open(sys.argv[-1]) as o_f:
    lines = o_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line for line in lines if line]
    begin = [i for i,l in enumerate(lines) if '>' in l]
    wdsps = [lines[b:e] for b,e in zip(begin,begin[1:]+[len(lines)])]
    long_wdsps = OrderedDict()
    for w in wdsps:
        pro = w[0].split()[1]
        if not pro in long_wdsps.keys():
            long_wdsps[pro] = w
        else:
            if len(w) > len(long_wdsps[pro]):
                print pro
                print len(long_wdsps[pro])
                print len(w)
                long_wdsps[pro] = w
    with open('select_long.wdsp','w') as w_f:
        for k,v in long_wdsps.iteritems():
            if len(v) >= 9:
                for vi in v:
                    print >> w_f,vi

