#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
relace strings in a file according to a dic file
"""

import os
import sys

def read_dic(dic_f):
    lines = dic_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line.split() for line in lines if len(line.split()) == 2]
    dic = dict()
    for line in lines:
        dic[line[0]]=line[1]
    return dic

def main():
    with open(sys.argv[-2]) as dic_f:
        dic = read_dic(dic_f)

    with open(sys.argv[-1]) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        newlines = [lines[0]]
        for line in lines[1:]:
            for k,v in dic.iteritems():
                if k in line:
                    line = line.replace(k,v)
                    break
            newlines.append(line)

    with open('relaced.txt','w') as w_f:
        for line in newlines:
            print >> w_f,line


if __name__ == "__main__":
    main()
