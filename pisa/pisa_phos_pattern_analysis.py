#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
analysis phos_binding_pattern in pdbs
"""

import os
import sys
import cPickle as pickle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from collections import Counter

AA = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', 'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}

def pattern_plot(phos_patterns,fname):
    # patterns format, [['TPO',['ARG','LYS'...]]...]
    # filter hetam
    def check_aa(pattern):
        for res in pattern:
            if not res in AA.keys():
                return False
        return True
    patterns = [(pdb,lig,p) for pdb,lig,p in phos_patterns if check_aa(p)]

    # plot amino acid freq
    aminoacids = []
    for _,_,p in patterns:
        aminoacids += p
    aa_count = Counter(aminoacids)
    aa = []
    count = []
    for a,c in aa_count.iteritems():
        c = c*1.0/len(aminoacids)
        aa.append(a)
        count.append(c)

    df = pd.DataFrame({'AA':aa,'Freq':count})
    df = df.sort_values('Freq',ascending=False)
    f,ax = plt.subplots()
    sns.set_style('whitegrid')
    sns.set_palette('pastel')
    sns.barplot(x='AA',y='Freq',data=df)
    plt.savefig(fname+'_aa_freq.png',dpi=300)

    # plot amino acid combination freq
    patterns = [(AA[pi] for pi in p) for _,_,p in patterns]
    patterns = [''.join(sorted(p)) for p in patterns]
    patterns_count = Counter(patterns)
    pattern = []
    count = []
    for p,c in patterns_count.iteritems():
        pattern.append(p)
        count.append(c)

    pattern_num = len(pattern)
    df = pd.DataFrame({'Pattern':pattern,'Freq':count})
    df = df.sort_values('Freq',ascending=False)
    f,ax = plt.subplots(figsize=(8,pattern_num*0.5))
    sns.set_style('whitegrid')
    sns.set_palette('pastel')
    sns.barplot(y='Pattern',x='Freq',data=df)
    plt.savefig(fname+'_pattern_freq.png',dpi=300)


def write_phos_patterns(phos_patterns,fname):
    with open(fname+'_patterns.txt','w') as w_f:
        for pdb,lig,res in phos_patterns:
            print >> w_f,'{0:<8}{1:<20}{2:<}'.format(pdb,lig,res)


def main():

    fname = '.'.join(os.path.split(sys.argv[-1])[1].split('.')[:-1])
    pdb_phos_sites = pickle.load(open(sys.argv[-1]))

    phos_sites = [(p[0],p[6],set(['_'.join(pi[1].split('_')[:3]) for pi in p[8]])) for p in pdb_phos_sites]

    phos_patterns = [(pdb,lig,[r.split('_')[2] for r in res]) for pdb,lig,res in phos_sites]
    write_phos_patterns(phos_patterns,fname+'_phos_patterns')

    pattern_plot(phos_patterns,fname)

    # tpo_patterns = [(pdb,lig,p) for pdb,lig,p in phos_patterns if lig == 'TPO']
    # pattern_plot(tpo_patterns,fname+'_TPO')

    # sep_patterns = [(pdb,lig,p) for pdb,lig,p in phos_patterns if lig == 'SEP']
    # pattern_plot(sep_patterns,fname+'_SEP')

    # ptr_patterns = [(pdb,lig,p) for pdb,lig,p in phos_patterns if lig == 'PTR']
    # pattern_plot(ptr_patterns,fname+'_PTR')

    with open(fname+'_phos_sites.txt','w') as w_f:
        for pdb,lig,res in phos_sites:
            print >> w_f,'{0:<8}{1:<12}{2:<}'.format(pdb,lig,res)


if __name__ == "__main__":
    main()
