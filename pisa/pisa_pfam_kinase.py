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

KINASE = ['PF00069']

def get_pdb_chain_pfam(p):
    pdbid,inter_id,inter_chain,inter_res = p
    inter_res_num = [int(r.split('_')[1]) for r in inter_res]

    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbid,
        'customReportColumns':'structureId,pfamAccession',
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
    pfam = [line[2] for line in lines if line[1] == inter_chain][0]

    inter_pfam = ''
    if pfam:
        pfams = pfam.split('#')
        for f in pfams:
            acc,res_range = f.split()
            res_range = res_range.strip('[]').split('-')
            res_range = map(int,res_range)
            if min(inter_res_num) > min(res_range) and max(inter_res_num) < max(res_range):
                inter_pfam = acc
    return pdbid,inter_id,inter_pfam

def filter_kinae(pdb_interfaces):

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

    p = Pool(4)
    result = p.map(get_pdb_chain_pfam,phos_inter)
    p.close()

    # result = []
    # for p in phos_inter:
        # result.append(get_pdb_chain_pfam(p))

    non_kinase_result = [(r[0],r[1]) for r in result if not r[2] in KINASE]
    kinase_result = [(r[0],r[1]) for r in result if r[2] in KINASE]


    non_kinae_pdb_interfaces = [p for p in pdb_interfaces if (p[0],p[1]) in non_kinase_result]
    kinae_pdb_interfaces = [p for p in pdb_interfaces if (p[0],p[1]) in kinase_result]

    with open('kinase_pdb_interface.txt','w') as w_f:
        for p in kinae_pdb_interfaces:
            print >> w_f,p

    print 'after filter kinae',len(non_kinae_pdb_interfaces)
    return non_kinae_pdb_interfaces

def main():

    fname = '.'.join(os.path.split(sys.argv[-1])[1].split('.')[:-1])
    pdb_interfaces = pickle.load(open(sys.argv[-1]))


    pdb_interfaces = filter_kinae(pdb_interfaces)
    pickle.dump(pdb_interfaces,open(fname+'_non_kinae.pickle','w'))

if __name__ == "__main__":
    main()


