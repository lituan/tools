#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pipeline
construct phylogenetic tree and estimate ancestors using msa
"""

import os
import sys
import subprocess
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool
from collections import OrderedDict
from scipy.stats import ttest_ind


def read_mega(mega_f):
    with open(mega_f) as o_f:

        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        pro_line_num = [i for i, line in enumerate(
            lines) if '>' in line] + [len(lines)]
        seqs = [lines[n:pro_line_num[i + 1]]
                for i, n in enumerate(pro_line_num[:-1])]
        seqs = [(seq[0][1:], ''.join(seq[1:])) for seq in seqs]

        seq_num = len(seqs)
        sequences = []
        i = 0
        for name,seq in seqs:
            if i < (seq_num + 1)/2:
                sequences.append(('current',name,seq))
            else:
                sequences.append(('ancestor',name,seq))
            i += 1
        return sequences

def get_similarity(seqs):
    # seqs format, [['ancestor',pro,'ATR...']...]
    def get_sim(r1,r2):
        return len([1 for i,r1i in enumerate(r1) if r1i == r2[i]])*1.0/len(r1)

    similarity = []
    for state,pro,seq in seqs:
        repeats = [seq[i:i+18] for i in range(0,len(seq),18)]
        sims = []
        for i in range(8):
            for j in range(8):
                if j > i:
                    sims.append(get_sim(repeats[i],repeats[j]))
        similarity.append((state,pro,sims))

    return similarity

def strip_plot(similarity,fname):
    # similarity format [['ancestor',pro,[0.3,04...]]...]
    p_label = []
    for i in range(28):
        anc = [si[i] for s,name,si in similarity if s == 'ancestor']
        now = [si[i] for s,name,si in similarity if s == 'current']
        pvalue = ttest_ind(anc,now,equal_var=False)
        if pvalue[-1] < 0.0001:
            p_label.append('****')
        elif pvalue[-1] < 0.001:
            p_label.append('***')
        elif pvalue[-1] < 0.01:
            p_label.append('**')
        elif pvalue[-1] < 0.05:
            p_label.append('*')
        else:
            p_label.append('')

    state = [ ]
    repeat_pair = []
    repeat_pair_similarity = []
    keys = [str(i)+'_'+str(j) for i in range(8) for j in range(8) if j > i]
    for s,name,sims in similarity:
        state += [s] * 28
        repeat_pair += keys
        repeat_pair_similarity += sims

    df = pd.DataFrame({'state':state,'repeat_pair':repeat_pair,'repeat_pair_similarity':repeat_pair_similarity})

    sns.set(style='whitegrid', color_codes=True)
    f,ax = plt.subplots(figsize=(12,8))
    # sns.stripplot(x='repeat_pair',y='repeat_pair_similarity',hue='state',palette='Set2',data=df,jitter=True)
    # sns.violinplot(x='repeat_pair',y='repeat_pair_similarity',hue='state',palette='Set2',data=df,split=True)
    sns.violinplot(x='repeat_pair',y='repeat_pair_similarity',hue='state',palette='muted',data=df,split=True)

    # add pvalue label
    anno_y = [max([si[i] for _,_,si in similarity])+0.10 for i in range(28) ]
    for i in range(28):
        ax.annotate(p_label[i],(i,anno_y[i]))
    # remove duplicated legends
    handles,labels = ax.get_legend_handles_labels()
    plt.legend(handles[:2],labels[:2])
    plt.savefig(fname+'_repeats_similarity.png',dpi=300)

    return p_label.count('*') + p_label.count('**') + p_label.count('***') + p_label.count('****')

def single_fun(mega_f):
    fpath,fname = os.path.split(mega_f)
    fname = os.path.splitext(fname)[0]
    seqs = read_mega(mega_f)
    similarity = get_similarity(seqs)
    f_name = os.path.join(fpath,fname)
    significant = strip_plot(similarity,f_name)
    return fname,significant


def main():
    parameters = []
    for root,dirs,files in os.walk(sys.argv[-1]):
        for f in files:
            if f[-17:] == 'most_probable.fas':
                f = os.path.join(root,f)
                parameters.append(f)

    # single_fun(parameters[0])
    p = Pool(6)
    result = p.map(single_fun,parameters)
    p.close()

    with open('significant.txt','w') as w_f:
        for r in result:
            print >> w_f,r[0],'\t',r[1]

    print [r for r in result if r[1] > 14]

if __name__ == "__main__":
    main()












