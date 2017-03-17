#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
get phos_binding residues and its interacting residues
"""

import os
import sys
import cPickle as pickle
from multiprocessing import Pool
from collections import OrderedDict


def select_phos_site(pdb_interfaces):
    pdb_phos_sites= []
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
                        phos_interacting_residues[res1] = [(bondi[0],res2)]
                    else:
                        phos_interacting_residues[res1].append((bondi[0],res2))
                elif [p for p in PHOS if res2[-8:] == p]:
                    res2 = '_'.join(res2.split('_')[:3])
                    if not res2 in phos_interacting_residues.keys():
                        phos_interacting_residues[res2] = [(bondi[1],res1)]
                    else:
                        phos_interacting_residues[res2].append((bondi[1],res1))

        for phos,interacting_residues in phos_interacting_residues.items():
            try:
                phos_chain = [c for c in chains if c[0] == phos.split('_')[0]][0]
                interacting_chain = [c for c in chains if c[0] == interacting_residues[0][1].split('_')[0]][0]
                phos_residue_num = int(phos.split('_')[1])
                phos_chain_residues_num  = [int(r.split('_')[1]) for r in phos_chain[4]]
                phos_position = min(max(phos_chain_residues_num)-phos_residue_num,phos_residue_num-min(phos_chain_residues_num))
                pdb_phos_sites.append((pdbid,p1,interface_area,p2,p3,p4,phos,phos_position,interacting_residues,phos_chain,interacting_chain))
            except:
                print phos
                print interacting_residues
                continue

    print 'num of phos sites',len(pdb_phos_sites)
    return pdb_phos_sites

def main():

    fname = '.'.join(os.path.split(sys.argv[-1])[1].split('.')[:-1])
    pdb_interfaces = pickle.load(open(sys.argv[-1]))
    pdb_phos_sites = select_phos_site(pdb_interfaces)

    pickle.dump(pdb_phos_sites,open(fname+'_phos_sites.pickle','w'))
    with open(fname+'_phos_sites.txt','w') as w_f:
        for site in pdb_phos_sites:
            print >> w_f,'{0:<6}{1:<6}{2:<8.2f}{3:<8.2f}{4:<6.2f}{5:<8.2f}{6:<16}{7:<6}\t{8:<}\t{9:<}\t{10:<}'.format(site[0],site[1],site[2],site[3],site[4],site[5],site[6],site[7],site[8],site[9][:4],site[10][:4])

    # for i in range(1,20):
        # phos_sites = [p for p in pdb_phos_sites if len(p[8]) == i ]
        # if phos_sites:
            # pickle.dump(pdb_phos_sites,open(fname+'_phos_sites_'+str(i)+'.pickle','w'))
            # with open(fname+'_phos_sites_'+str(i)+'.txt','w') as w_f:
                # for site in phos_sites:
                    # print >> w_f,'{0:<6}{1:<6}{2:<8.2f}{3:<8.2f}{4:<6.2f}{5:<8.2f}{6:<16}{7:<6}\t{8:<}\t{9:<}\t{10:<}'.format(site[0],site[1],site[2],site[3],site[4],site[5],site[6],site[7],site[8],site[9],site[10])

if __name__ == "__main__":
    main()




