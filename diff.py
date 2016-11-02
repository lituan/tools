"""
find common and special lines in files
"""
import os
import sys

files = sys.argv[1:]
file_lines = []
for i,f in enumerate(files):
    f_path,f_full_name = os.path.split(f)
    f_name,f_exten = os.path.splitext(f_full_name)

    with open(f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\n\r') for line in lines]
        lines = [line for line in lines if line]
        lines = set(lines) # delete redundant lines
        lines = sorted(lines)
        file_lines.append([f_path,f_name,lines])

#find common lines
common = file_lines[0][2]
common = sorted(common)
for i in range(len(files)-1):
    common = set(common).intersection(file_lines[i+1][2])
with open('common.txt','w') as w_f:
    for i in common:
        print >> w_f,i

total = [line for _,_,lines in file_lines for line in lines]
total = set(total)
total = sorted(total)
with open('total.txt','w') as w_f:
    for i in total:
        print >> w_f,i

for f_path,f_name,lines in file_lines:
    special = set(lines).difference(common)
    special = sorted(special)
    with open(f_name+"_"+'special.txt','w') as w_f:
        for i in special:
            print >> w_f,i

