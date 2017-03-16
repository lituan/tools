#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
analys pdb_lig_interactions and make graphs
"""

import os
import sys
import cPickle as pickle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from collections import Counter

def get_pdb_chain_sequences(p):
    # p format: [(pdbid,chain)...]
    pdbids = ','.join([pdbid for pdbid,chain in p])
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbids,
        'customReportColumns':'structureId,uniprotAcc,entityId,resolution,chainLength,releaseDate,sequence',
        'service':'wsfile',
        'format':'csv',
    }
    data = urllib.urlencode(data)
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    lines = response.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line for line in lines if line]
    lines = [line.split(',') for line in lines]
    lines = [[w.strip('"') for w in line] for line in lines]
    lines = [(line[0],line[1],line[-1])for line in lines if (line[0],line[1]) in p]
    lines = [line for line in lines if line[-1]]
    return lines


def patter_analysis(pdb_lig_interactions):

    # salt_bridges and hydrogenbonds
    pdb_lig_interaction_res= []
    for pdb,lig,interactions in pdb_lig_interactions:
        inter_res = []
        for interaction_type,interaction_res in interactions:
            inter_res += interaction_res
        if inter_res:
            pdb_lig_interaction_res.append((pdb,lig,set(inter_res)))

    patterns = [(pdb,lig,[r.split('_')[-1] for r in res]) for pbd,lig,res in pdb_lig_interaction_res]
    patterns = [p[2] for p in patterns if len(p[2])]
    # patterns = [p[2] for p in patterns if len(p[2]) >= 3]

    aminoacids = []
    for p in patterns:
        aminoacids += p
    aa_count = Counter(aminoacids)
    aa = []
    count = []
    for a,c in aa_count.iteritems():
        aa.append(a)
        count.append(c)

    fname = os.path.split(sys.argv[-1])[1].split('.')[0]
    df = pd.DataFrame({'AA':aa,'Freq':count})
    df = df.sort_values('Freq',ascending=False)
    f,ax = plt.subplots()
    sns.set_style('whitegrid')
    sns.set_palette('pastel')
    sns.barplot(x='AA',y='Freq',data=df)
    plt.savefig(fname+'_aa_freq.png',dpi=300)

    aa = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', 'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}
    patterns = [(aa[pi] for pi in p) for p in patterns]
    patterns = [''.join(sorted(p)) for p in patterns]
    patterns_count = Counter(patterns)
    pattern = []
    count = []
    for p,c in patterns_count.iteritems():
        pattern.append(p)
        count.append(c)

    fname = os.path.split(sys.argv[-1])[1].split('.')[0]
    df = pd.DataFrame({'Pattern':pattern,'Freq':count})
    df = df.sort_values('Freq',ascending=False)
    f,ax = plt.subplots()
    sns.set_style('whitegrid')
    sns.set_palette('pastel')
    sns.barplot(y='Pattern',x='Freq',data=df)
    plt.savefig(fname+'_pattern_freq.png',dpi=300)


def main():
    # format ['1nex','C_101_PTR','A',[('hydrophobic',['A_206_LYS'])]]
    pdb_lig_interactions = pickle.load(open(sys.argv[-1],'r'))
    patter_analysis(pdb_lig_interactions)


if __name__ == "__main__":
    main()
