#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
get names for files in a directory
"""
import os
import sys

def files_in_dir(directory):
    for root, dirs, files in os.walk(directory):
        for f in files:
            yield os.path.join(root, f)

names = []
for f in files_in_dir(sys.argv[-1]):
    f_path,f_name = os.path.split(f)
    f_name,f_extention = os.path.splitext(f_name)
    names.append(f_name)

names = sorted(names)
with open('file_name.txt','w') as w_f:
    for n in names:
        print >> w_f, n
