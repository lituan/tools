#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
delete redundant lines
"""
import sys
import os

with open(sys.argv[-1]) as o_f:
    lines = o_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    print 'total lines: ',len(lines)
    lines_unique = [line for line in lines if not line in lines_unique]
    print 'total unique lines: ',len(lines_unique)

    with open('unique.txt','w') as w_f:
        for line in lines_unique:
            print >> w_f, line
