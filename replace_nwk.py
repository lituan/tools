#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
change node names in a nwk file according to a dic file
"""

import os
import sys

def read_dic(dic_f):
    lines = dic_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line.split() for line in lines if len(line.split()) == 2]
    dic = dict()
    for line in lines:
        dic[line[1]]=line[0]
    return dic

def main():
    with open(sys.argv[-2]) as dic_f:
        dic = read_dic(dic_f)

    with open(sys.argv[-1]) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        line = lines[0]
        for k,v in dic.iteritems():
            print k,v
            if k in line:
                print k
                line = line.replace(k,v)

    with open('relaced.nwk','w') as w_f:
        print >> w_f,line


if __name__ == "__main__":
    main()
