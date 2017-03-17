#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
filter phos pattterns in kinase
"""

import os
import sys
import urllib
import urllib2
import cPickle as pickle
from collections import Counter
from multiprocessing import Pool


def get_pdb_chain_pfam(p):
    pdbid,inter_id,inter_chain,inter_res = p
    inter_res_num = [int(r.split('_')[1]) for r in inter_res]

    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbid,
        'customReportColumns':'structureId,entityId,pfamAccession,pfamId',
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
    inter_chain_id = [line for line in lines if line[1] == inter_chain][0][2]

    pfam = [(line[3],line[4]) for line in lines if line[2] == inter_chain_id if line[3]]

    inter_pfam = []
    if pfam:
        pfam_accs,pfam_ids = pfam[0]
        for pfam_a,pfam_i in zip(pfam_accs.split('#'),pfam_ids.split('#')):
            acc,res_range = pfam_a.split()
            res_range = res_range.strip('[]').split('-')
            res_range = map(int,res_range)
            for res_num in inter_res_num:
                if res_num >= min(res_range) and res_num <= max(res_range):
                    inter_pfam.append((acc,pfam_i))
    return pdbid,inter_id,inter_pfam

def get_all_pfam(pdb_interfaces):

    phos_inter = []
    for interface in pdb_interfaces:
        pdbid,p1,interface_area,p2,p3,p4,bonds = interface[:7]
        chains = interface[-1]
        phos_interacting_residues = {}
        PHOS = ['TPO_ O1P','TPO_ O2P','TPO_ O3P','TPO_ OG1','SEP_ O1P','SEP_ O2P','SEP_ O3P','SEP_ OG ','PTR_ O1P','PTR_ O2P','PTR _O3P','PTR OH ']
        for bond in bonds:
            bond_type,bond_info = bond
            for bondi in bond_info:
                res1,res2,dist = bondi
                if [p for p in PHOS if res1[-8:] == p]:
                    res1 = '_'.join(res1.split('_')[:3])
                    if not res1 in phos_interacting_residues.keys():
                        phos_interacting_residues[res1] = [res2]
                    else:
                        phos_interacting_residues[res1].append(res2)
                elif [p for p in PHOS if res2[-8:] == p]:
                    res2 = '_'.join(res2.split('_')[:3])
                    if not res2 in phos_interacting_residues.keys():
                        phos_interacting_residues[res2] = [res1]
                    else:
                        phos_interacting_residues[res2].append(res1)

        for phos,interacting_residues in phos_interacting_residues.items():
            inter_chain = interacting_residues[0].split('_')[0]
            phos_inter.append([pdbid,p1,inter_chain,interacting_residues])

    # p = Pool(4)
    # result = p.map(get_pdb_chain_pfam,phos_inter)
    # p.close()

    # result = []
    # for p in phos_inter:
        # result.append(get_pdb_chain_pfam(p))

    # pickle.dump(result,open('pfam_result.pickle','w'))
    result = pickle.load(open('pfam_result.pickle'))

    return result

    all_pfams = []
    for r in result:
        if r[2]:
            for p in r[2]:
                all_pfams.append(p[1])
    all_pfams_counter =  Counter(all_pfams)
    print all_pfams_counter
    all_pfams_counter = [(k,v) for k,v in all_pfams_counter.iteritems()]
    all_pfams_counter = sorted(all_pfams_counter,key=lambda x: x[1],reverse=True)
    print all_pfams_counter
    with open('pfams.txt','w') as w_f:
        for c in all_pfams_counter:
            print >> w_f,'{0:<30}{1:<}'.format(c[0],c[1])

    with open('pfam_result.txt','w') as w_f:
        for r in result:
            print >> w_f,r


def main():

    fname = '.'.join(os.path.split(sys.argv[-1])[1].split('.')[:-1])
    pdb_interfaces = pickle.load(open(sys.argv[-1]))

    pfam_result = get_all_pfam(pdb_interfaces)

    all_pfams = []
    for r in pfam_result:
        if r[2]:
            for p in r[2]:
                all_pfams.append(p[1])
    all_pfams_counter =  Counter(all_pfams)

    all_pfams_counter = [(k,v) for k,v in all_pfams_counter.iteritems()]
    all_pfams_counter = sorted(all_pfams_counter,key=lambda x: x[1],reverse=True)

    with open(fname+'_pfam_table.txt','w') as w_f:
        for c in all_pfams_counter:
            print >> w_f,'{0:<30}{1:<}'.format(c[0],c[1])

    with open(fname+'_pfam_result.txt','w') as w_f:
        for r in pfam_result:
            print >> w_f,r

    with open(fname+'_different_pfam.txt','w') as w_f:
        for r in pfam_result:
            r_pfams = [p[1] for p in r[2]]
            if len(set(r_pfams)) > 1:
                print >> w_f,r


if __name__ == "__main__":
    main()



