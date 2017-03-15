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

AA = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', 'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}

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

def pattern_plot(patterns,fname):
    # patterns format, [['TPO',['ARG','LYS'...]]...]
    # filter hetam
    def check_aa(pattern):
        for res in pattern:
            if not res in AA.keys():
                return False
        return True
    patterns = [(lig,p) for lig,p in patterns if check_aa(p)]

    # plot amino acid freq
    aminoacids = []
    for _,p in patterns:
        aminoacids += p
    aa_count = Counter(aminoacids)
    aa = []
    count = []
    for a,c in aa_count.iteritems():
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
    patterns = [(AA[pi] for pi in p) for _,p in patterns]
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


def patter_analysis(pdb_interfaces,fname):
    pdb_lig_interaction_res= []
    for interface in pdb_interfaces:
        pdbid,p1,interface_area,p2,p3,p4,bonds = interface[:7]
        chains = interface[-1]
        phos_res = []
        interaction_res = []
        for bond in bonds:
            bond_type,bond_info = bond
            for bondi in bond_info:
                res1,res2,dist = bondi
                if not ('TPO' in res1 or 'SEP' in res1 or 'PTR' in res1):
                    if res1.split('_')[-1] != ' N  ':
                        interaction_res.append('_'.join(res1.split('_')[:3]))
                    else:
                        interaction_res.append('_'.join(res1.split('_')[:2]+['GLY']))
                if not ('TPO' in res2 or 'SEP' in res2 or 'PTR' in res2):
                    if res2.split('_')[-1] != ' N  ':
                        interaction_res.append('_'.join(res2.split('_')[:3]))
                    else:
                        interaction_res.append('_'.join(res2.split('_')[:2]+['GLY']))
                if 'TPO' in res1 or 'SEP' in res1 or 'PTR' in res1:
                    phos_res.append('_'.join(res1.split('_')[:3]))
                if 'TPO' in res2 or 'SEP' in res2 or 'PTR' in res2:
                    phos_res.append('_'.join(res2.split('_')[:3]))

        if len(set(phos_res)) != 1:
            continue
        if len(set([r.split('_')[0] for r in interaction_res])) != 1:
            continue
        phos_res = phos_res[0]
        interaction_res = list(set(interaction_res))
        pdb_lig_interaction_res.append((pdbid,phos_res,interaction_res))

    with open(fname+'_pdb_lig_interaction_res.txt','w') as w_f:
        for pdbid,phos_res,interaction_res in pdb_lig_interaction_res:
            print >> w_f, '{0:<8}{1:<15}{2:<}'.format(pdbid,phos_res,' '.join(interaction_res))
            interaction_res = [AA[r.split('_')[-1]] for r in interaction_res]
            interaction_res = ''.join(sorted(interaction_res))
            print >> w_f, '{0:<8}{1:<15}{2:<}'.format(pdbid,phos_res,interaction_res)

    patterns = [(lig.split('_')[-1],[r.split('_')[-1] for r in res]) for pdb,lig,res in pdb_lig_interaction_res]

    pattern_plot(patterns,fname)

    tpo_patterns = [(lig,p) for lig,p in patterns if lig == 'TPO']
    pattern_plot(patterns,fname+'_TPO')

    tpo_patterns = [(lig,p) for lig,p in patterns if lig == 'SEP']
    pattern_plot(patterns,fname+'_SEP')

    tpo_patterns = [(lig,p) for lig,p in patterns if lig == 'PTR']
    pattern_plot(patterns,fname+'_PTR')



def main():
    # format ['1nex','C_101_PTR','A',[('hydrophobic',['A_206_LYS'])]]
    fname = '.'.join(os.path.split(sys.argv[-1])[-1].split('.')[:-1])
    pdb_interfaces = pickle.load(open(sys.argv[-1],'r'))
    patter_analysis(pdb_interfaces,fname)


if __name__ == "__main__":
    main()
