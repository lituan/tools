#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
copy a list of files to new directory
"""
import os
import sys
from shutil import copyfile

def files_in_dir(directory):
    for root, dirs, files in os.walk(directory):
        for f in files:
            yield os.path.join(root, f)

with open(sys.argv[-2]) as o_f:
    lines = o_f.readlines()
    ids = [line.split()[0] for line in lines if len(line) > 0]

os.makedirs('new_files')

for f in files_in_dir(sys.argv[-1]):
    f_path,f_full_name = os.path.split(f)
    f_short_name,f_extention = os.path.splitext(f_full_name)
    for i in ids:
        if i == f_short_name:
            dir_file = os.path.join("new_files",f_full_name)
            copyfile(f,dir_file)




