#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
change extention of files in a directory
for example, from '.out' to '.wdsp'
"""
import sys
import os
import re


for f in lt.files_in_dir(sys.argv[-1]):
    with open(f) as o_f:
        lines = o_f.readlines()
        new_name = o_f.name
        pattern = re.compile(r'(?<=\d{%).*(?=.out)')
        new_name = pattern.findall(new_name)[0]
        new_name = new_name + '.wdsp'
        with open(new_name,'w') as w_f:
            for line in lines:
                print >> w_f,line
