#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
given a matrix of pairwise similarity, output the average similarity
the matrix must have lables on the top and right,as the follow shows
  lables    id1 id2
  id1       1.0  0.2
  id2       0.2  1.0
"""
import sys
import os
import string
import operator
import lt


def readmatrix(m_f):

    f_path,f_name = os.path.split(m_f)

    with open(m_f) as o_f:
        lines = o_f.readlines()
        lines = [line.strip('\r\n') for line in lines]  # remove EOF
        column_header = lines[0].split()[1:]
        row_header = column_header
        matrix = [line.split()[1:] for line in lines[1:]]
        # print matrix
        matrix = [map(float, row) for row in matrix]


        # transform to flat dist array
        m_len = len(matrix)
        flat_matrix = []
        for i in range(m_len):
            for j in range(m_len):
                if j > i:
                    flat_matrix.append(matrix[i][j])

        average = sum(flat_matrix)*1.0/len(flat_matrix)

    words = f_name.split('_')
    number = int(words[0])
    shape = '_'.join(words[1:-3])
    seq = words[-3]

    return [number,shape,seq,average]

def main():

    average = []
    for f in lt.files_in_dir(sys.argv[-1]):
        if 'align' in f:
            result = readmatrix(f)
            average.append(result)

    with lt.open_file() as w_f:
        average = sorted(average,key=operator.itemgetter(0),reverse=True)
        for a in average:
            print >> w_f,'{0:<10}{1:<20}{2:<10}{3:<10.8f}'.format(a[0],a[1],a[2],a[3])

if __name__ == "__main__":
    main()
