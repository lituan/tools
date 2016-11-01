"""
put ids in 'one id one line' format
"""
import sys

with open(sys.argv[-1]) as o_f:
    lines = o_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line.split(';') for line in lines]
    ids = [i for line in lines for i in line]
    ids = set(ids)
    ids = sorted(ids)

    with open('id_in_line.txt','w') as w_f:
        for i in ids:
            print >> w_f,i

