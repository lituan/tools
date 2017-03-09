#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
classify wdsps by blades num
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
    wdsps_class = {}
    for w in wdsps:
        blade_num = len(w) - 1
        if not blade_num in wdsps_class.keys():
            wdsps_class[blade_num] = [w]
        else:
            wdsps_class[blade_num].append(w)

    for blade_num,pros in wdsps_class.iteritems():
        print 'blade_num: ',blade_num,' ',len(pros)
        with open(str(blade_num)+'_blades.wdsp','w') as w_f:
            for p in pros:
                for line in p:
                    print >> w_f,line

