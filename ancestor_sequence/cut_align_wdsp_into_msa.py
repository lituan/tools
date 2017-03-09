#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
cut  a,b,c in each blade
"""

import os
import sys


def trans_lis_lis(lis_lis):
    """align and trans nested list to print a table"""
    lis_lis = [[str(l) for l in lis]
               for lis in lis_lis]  # trans every element to str
    lis_lis_len = len(lis_lis)
    lis_max_len = max(len(lis) for lis in lis_lis)
    # construct a matrix
    lis_lis = [lis + (lis_max_len - len(lis)) * [''] for lis in lis_lis]
    # transform the matrix to its T matrix
    trans = [[lis_lis[i][j]
              for i in range(lis_lis_len)] for j in range(lis_max_len)]
    # aligned = [[l + (max([len(l) for l in lis]) - len(l)) * " " for l in lis] for lis in trans]
    aligned = []
    for lis in trans:
        width = max([len(l) for l in lis])
        lis = [l + (width - len(l)) * '-' for l in lis]
        aligned.append(lis)
    # transform the T matrix to the original  matrix
    trans = [[aligned[i][j]
              for i in range(lis_max_len)] for j in range(lis_lis_len)]
    return trans


# check blade_num and completeness of each balde
def check_entry(entry,blade_num):
    if len(entry) != blade_num + 1:
        return 0
    for e in entry[1:-1]:
        if len(e) != 12:
            return 0
    if len(entry[-1]) != 11:
        return 0
    return 1

def main():

    with open(sys.argv[-1]) as o_f:
        lines = o_f.readlines()
        # strip end of line symbol and empty lines
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        # cut into entries
        begin = [i for i, line in enumerate(lines) if '>' in line]
        end = begin[1:] + [len(lines)]
        entries = [lines[begin[i]:end[i]] for i in range(len(begin))]
        entries = [[l.split() for l in entry] for entry in entries]
        # filter non-complete entries
        entries = [entry for entry in entries if check_entry(entry,8)]
        # cut a b c
        entries = [[entry[0]] + [[e[5],e[7],e[9]] for e in entry[1:]] for entry in entries]
        # align each entry
        entries_new = []
        for entry in entries:
            combine = ['{0:<40}'.format(entry[0][1])]
            for e in entry[1:]:
                for ei in e:
                    combine += ei
            entries_new.append(combine)
        entries_new = trans_lis_lis(entries_new)

        msa = [[s[0].strip(' '),''.join(s[1:])]for s in entries_new]
        fname = os.path.split(sys.argv[-1])[1].split('.')[0]
        with open(fname+'cut_align_wdsp_into_msa.fas', 'w') as w_f:
            for pro,seq in msa:
                print >> w_f,'>{0}'.format(pro)
                s = [seq[i:i+80] for i in range(0,len(seq),80)]
                for si in s:
                    print >> w_f,si

        with open(fname+'cut_align_wdsp_into_msa_seq.fas', 'w') as w_f:
            for pro,seq in msa:
                seq = seq.replace('-','')
                print >> w_f,'>{0}'.format(pro)
                seqs = [seq[i:i+80] for i in range(0,len(seq),80)]
                for s in seqs:
                    print >> w_f,s
        with open(fname+'_ids.txt','w') as w_f:
            for pro,seq in msa:
                print >> w_f,pro


if __name__ == "__main__":
    main()
