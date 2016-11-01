#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
get sequence id from a fasta file
"""
import os
import sys

with open(sys.argv[-1]) as w_f:
    lines = w_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line.split()[0] for line in lines if len(line.split()) > 0]
    ids = [line for line in lines if '>' in line]
    # ids = [i.split()[0].split('|')[2] for i in ids ] # for fa from uniprot
    ids = [i.split('>')[-1] for i in ids]
    print 'total num of id:',len(ids)

    with open('ids.txt','w') as w_f:
        for i in ids:
            print >> w_f,i

