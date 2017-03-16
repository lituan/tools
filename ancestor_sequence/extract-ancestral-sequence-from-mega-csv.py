#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
input mega ancestral sequences estimation in csv format, output all probable ancestral sequences ranked according to probability
you can set CUTOFF for the probability to remove some AAs with very low probability
it can consume lots of memory for large datasets, be careful
"""

import sys
import os
from collections import OrderedDict
from itertools import product
from memprof import memprof

CUTOFF = 0.05
NODE_CUTOFF = 30

def align_lis_lis(lis_lis):
    """align and trans nested list to print a table"""
    lis_lis = [[str(l) for l in lis]
               for lis in lis_lis]  # trans every element to str
    #make all inner lists of the same length
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [lis + (inner_lis_max_len - len(lis)) * [''] for lis in lis_lis]
    #trans list, so that the elements of the same column are in one list
    lis_lis = [[lis[i] for lis in lis_lis] for i in range(inner_lis_max_len)]
    #make element in the same list have the same length
    aligned = []
    for lis in lis_lis:
        width = max([len(l) for l in lis])
        lis = [l + (width - len(l)) * ' ' for l in lis]
        aligned.append(lis)
    #trans list_list to the original list_list
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [[lis[i] for lis in aligned] for i in range(inner_lis_max_len)]
    return lis_lis


def read_csv(mega_f):
    with open(mega_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        lines = [line.split(',') for line in  lines]
        nodes = lines[0][2:]
        len_seq = len(lines[1:])/20
        aalist = [line[1] for line in lines[1:21]]

        nodes_aa = OrderedDict()
        for i,node in enumerate(nodes):
            aa = []
            for j in range(len_seq):
                aa_j = []
                aa_probability = map(float,[line[i+2] for line in lines[j*20+1:j*20+21]])
                for a,p in zip(aalist,aa_probability):
                    if p > CUTOFF:
                        aa_j.append((a,p))
                aa_j = sorted(aa_j, key= lambda x: x[1],reverse=True)
                aa.append(aa_j)
            nodes_aa[node] = aa
        # print nodes_aa['Node_56']
        return nodes_aa




def write_nodes(nodes_aa):
    most_probable = []
    for node, aa in nodes_aa.iteritems():
        new_aa = []
        for a in aa:
            if a:
                new_aa.append(a)
            else:
                new_aa.append([(' ', 1.0)])
        new_aa = [a[0] for a in new_aa]
        aa_a = [a[0] for a in new_aa]
        aa_p = ['{0:<4.2f}'.format(a[1]) for a in new_aa]
        most_probable.append([node]+aa_p)
        most_probable.append([node]+aa_a)
    most_probable = align_lis_lis(most_probable)
    with open('most_probable.txt', 'w') as w_f:
        for i, mp in enumerate(most_probable):
            print >> w_f, ' '.join(mp)
            if (i+1)%2 == 0:
                print >> w_f, ' '

    combinations = OrderedDict()
    with open('all_probable_ancestor_sequences.txt','w') as w_f:
        for node, aa in nodes_aa.iteritems():
            node_id = int(node.split('_')[1])
            if node_id <= NODE_CUTOFF:
                seq = []
                for a in aa:
                    if a:
                        seq.append(a[0][0])
                    else:
                        seq.append('-')
                print >> w_f,node
                print >> w_f,''.join(seq)


            if node_id > NODE_CUTOFF:
                new_aa = []
                for a in aa:
                    if a:
                        new_aa.append(a)
                    else:
                        new_aa.append([(' ', 1.0)])
                aa = new_aa
                aa_combinations = product(*aa)

                aa_p = []
                for ac in aa_combinations:
                    ac_a = ['{0:<5}'.format(a[0]) for a in ac]
                    ac_p = ['{0:<4.3f}'.format(a[1]) for a in ac]
                    probability = eval('*'.join(map(str,ac_p)))
                    aa_p.append((probability, ac_p, ac_a))
                aa_p = sorted(aa_p,reverse=True)

                print >> w_f, node
                for probability,ac_p,ac_a in aa_p:
                    probability = '{0:<4.3f}'.format(probability)
                    print >> w_f, '{0:<8}{1}'.format(probability,' '.join(ac_a))
                    print >> w_f, '{0:<8}{1}'.format(' ',' '.join(ac_p))


def main():
    nodes_aa = read_csv(sys.argv[-1])
    write_nodes(nodes_aa)



if __name__ == "__main__":
    main()
