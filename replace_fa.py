"""
relace strings in a file according to a dic file
"""

import os
import sys
from collections import OrderedDict

def read_dic(dic_f):
    with open(dic_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        dic = dict()
        for line in lines:
            words = line.split()
            if len(words[0].split(',')) > 1:
                for i in words[0].split(','):
                    dic[i] = words[1]
            else:
                dic[words[0]] = words[1]
        return dic

def readfa(fa_f):
    # readin seqs in fasta format
    # seqs foramt:[(pro,seq),...]
    with open(fa_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        pro_line_num = [i for i, line in enumerate(lines) if '>' in line] + [len(lines)]
        seqs = [lines[n:pro_line_num[i + 1]]
                for i, n in enumerate(pro_line_num[:-1])]
        seqs = [(seq[0][1:], ''.join(seq[1:])) for seq in seqs]
        return seqs

def main():
    dic = read_dic(sys.argv[-2])

    seqs = readfa(sys.argv[-1])

    new_seqs = []
    for pro,seq in seqs:
        if dic.has_key(pro):
            new_seqs.append((dic[pro],seq))
        else:
            print 'not found:',pro

    with open('replaced.fa','w') as w_f:
        for pro,seq in new_seqs:
            print >> w_f,'>{0}'.format(pro)
            print >> w_f,seq


if __name__ == "__main__":
    main()
