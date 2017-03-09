#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
input wdsp out file, output msa by aligning  by structures
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

def align_score(msa,column_cutoff=0.95,score_cutoff=0.80,bad_columns_cutoff=5):
    """
    sequence align score = matched aa /all aa
    matched aa is in a column, 95% is not indel
    msa format [['pro','ACT-G...']...]
    """

    # filter bad sequence, which has many bad columns (have many '-')
    seq_len = len(msa[0][1])
    columns = [[s[1][i] for s in msa] for i in range(seq_len)]
    columns = [i for i,c in enumerate(columns) if len([ci for ci in c if ci == '-'])*1.0/len(c) > column_cutoff]
    msa_bad_columns = [[len([s[1][i] for i in columns if s[1][i] != '-']),s[0],s[1]] for s in msa ]
    msa = [(pro,seq) for bad_columns,pro,seq in msa_bad_columns if bad_columns < bad_columns_cutoff]
    # find columns that are all '-', and elimit them
    columns = [[s[1][i] for s in msa] for i in range(seq_len)]
    columns = [i for i,c in enumerate(columns) if len([ci for ci in c if ci == '-'])*1.0/len(c) < 1.0]
    msa = [(pro,''.join([seq[i] for i in columns])) for pro,seq in msa]

    # filter bad sequence, which has '-' in good columns (have little '-')
    seq_len = len(msa[0][1])
    columns = [[s[1][i] for s in msa] for i in range(seq_len)]
    columns = [i for i,c in enumerate(columns) if len([ci for ci in c if ci != '-'])*1.0/len(c) > column_cutoff]
    msa_bad_columns = [[len([s[1][i] for i in columns if s[1][i] == '-']),s[0],s[1]] for s in msa ]
    msa = [(pro,seq) for bad_columns,pro,seq in msa_bad_columns if bad_columns < bad_columns_cutoff]
    # find columns that are all '-', and elimit them
    columns = [[s[1][i] for s in msa] for i in range(seq_len)]
    columns = [i for i,c in enumerate(columns) if len([ci for ci in c if ci == '-'])*1.0/len(c) < 1.0]
    msa = [(pro,''.join([seq[i] for i in columns])) for pro,seq in msa]

    # calculate score
    seq_len = len(msa[0][1])
    columns = [[s[1][i] for s in msa] for i in range(seq_len)]
    columns = [i for i,c in enumerate(columns) if len([ci for ci in c if ci != '-'])*1.0/len(c) > column_cutoff]
    msa_scores = [[len([s[1][i] for i in columns if s[1][i] != '-'])*1.0/len([si for si in s[1] if si != '-']),s[0],s[1]] for s in msa ]
    msa = [(pro,seq) for score,pro,seq in msa_scores if score > score_cutoff]
    # find columns that are '-', and elimit them
    columns = [[s[1][i] for s in msa] for i in range(seq_len)]
    columns = [i for i,c in enumerate(columns) if len([ci for ci in c if ci == '-'])*1.0/len(c) < 1.0]
    msa = [(pro,''.join([seq[i] for i in columns])) for pro,seq in msa]
    return msa


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
        # check last line of each entry
        for entry in entries:
            if len(entry[-1]) < len(entry[1]):
                complement = len(entry[1]) - len(entry[-1])
                entry[-1] = entry[-1][0:-1] + \
                    complement * [' '] + [entry[-1][-1]]
        # align each entry
        entries_new = []
        for entry in entries:
            combine = ['{0:<40}'.format(entry[0][1])]
            for b in [blade[3:-1] for blade in entry[1:]]:
                combine +=  b
            entries_new.append(combine)
        entries_new = trans_lis_lis(entries_new)

        column_cutoff = 0.95
        score_cutoff = 0.80
        bad_columns_cutoff = 5
        for bad_columns_cutoff in [1,2,3,4,5]:
            msa = [[s[0].strip(' '),''.join(s[1:])]for s in entries_new]
            msa = align_score(msa,column_cutoff,score_cutoff,bad_columns_cutoff)

            print 'origin: ',len(entries_new)
            print 'now: ',len(msa)
            fname = os.path.split(sys.argv[-1])[1].split('.')[0] + '_' + str(column_cutoff) +'_'+str(score_cutoff)+'_'+str(bad_columns_cutoff)
            with open(fname+'_align_wdsp_into_msa.fas', 'w') as w_f:
                for pro,seq in msa:
                    print >> w_f,'>{0}'.format(pro)
                    s = [seq[i:i+80] for i in range(0,len(seq),80)]
                    for si in s:
                        print >> w_f,si

            with open(fname+'_align_wdsp_into_msa_seq.fas', 'w') as w_f:
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
