#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check same interface phos_binding patterns
"""

import os
import sys
import urllib
import urllib2
import cPickle as pickle
from multiprocessing import Pool


def get_entityid(p):
    pdbid,interface_id,chain1,chain2 = p
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbid,
        'customReportColumns':'structureId,entityId',
        'service':'wsfile',
        'format':'csv',
    }
    data = urllib.urlencode(data)
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    lines = response.readlines()
    lines = [line.rstrip('\r\n') for line in lines[1:]]
    lines = [line for line in lines if line]
    lines = [line.split(',') for line in lines]
    lines = [[w.strip('"') for w in line] for line in lines]
    chain1_id = [line for line in lines if line[1] == chain1][0][2]
    chain2_id = [line for line in lines if line[1] == chain1][0][2]
    return pdbid,interface_id,chain1_id,chain2_id

def filter_same_interface(pdb_interfaces):

    pdbid_chain = [(p[0],p[1],p[-1][0][0],p[-1][1][0]) for p in pdb_interfaces]

    p = Pool(4)
    result = p.map(get_entityid,pdbid_chain)
    p.close()

    pdb_chain_entity = {}
    for r in result:
        if not (r[0],r[2],r[3]) in pdb_chain_entity.keys():
            pdb_chain_entity[(r[0],r[2],r[3])] = [r]
        else:
            pdb_chain_entity[(r[0],r[2],r[3])].append(r)

    with open('same_interface.txt','w') as w_f:
        same = []
        different = []
        for k,v in pdb_chain_entity.iteritems():
            if len(v) > 1:
                print >> w_f,k
                cluster = [p for p in pdb_interfaces if (p[0],p[1]) in [(vi[0],vi[1]) for vi in v]]
                cluster_patterns = []
                for c in cluster:
                    bonds = c[6]
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
                        if interacting_residues:
                            interacting_residues = ['_'.join(r.split('_')[:3]) for r in interacting_residues]
                            interacting_residues = list(set(interacting_residues))
                            interacting_residues = [r.split('_')[2] for r in interacting_residues]
                            interacting_residues = sorted(interacting_residues)
                            interacting_residues = '_'.join(interacting_residues)
                            cluster_patterns.append(interacting_residues)
                            print >> w_f,c[0],c[1],interacting_residues

                print cluster_patterns
                if len(cluster_patterns) > 1 and len(set(cluster_patterns)) == 1:
                    same.append(1)
                else:
                    different.append(1)
        print 'same',len(same)
        print 'different',len(different)




    pdb_unique_interface = [(v[0][0],v[0][1]) for k,v in pdb_chain_entity.iteritems()]
    pdb_interfaces = [p for p in pdb_interfaces if (p[0],p[1]) in pdb_unique_interface]
    print 'after filter same entity',len(pdb_interfaces)
    return pdb_interfaces

def filter_non_one_phos(pdb_interfaces):
    zero_phos_interfaces = []
    one_phos_interfaces = []
    more_phos_interfaces = []
    for interface in pdb_interfaces:
        pdbid,p1,interface_area,p2,p3,p4,bonds = interface[:7]
        phos_res = []
        for bond in bonds:
            bond_type,bond_info = bond
            for bondi in bond_info:
                res1,res2,dist = bondi
                if 'TPO' in res1 or 'SEP' in res1 or 'PTR' in res1:
                    phos_res.append('_'.join(res1.split('_')[:3]))
                if 'TPO' in res2 or 'SEP' in res2 or 'PTR' in res2:
                    phos_res.append('_'.join(res2.split('_')[:3]))
        phos_res = set(phos_res)
        if len(phos_res) == 1:
            one_phos_interfaces.append(interface)
        elif len(phos_res) > 1:
            more_phos_interfaces.append(interface)
        else:
            zero_phos_interfaces.append(interface)

    print 'after filter non_one_phos_interfaces',len(one_phos_interfaces)
    return one_phos_interfaces

def main():
    pdb_interfaces = pickle.load(open(sys.argv[-1]))

    pdb_interfaces = [p for p in pdb_interfaces if p[7][0][2].lower() == 'x,y,z' and p[7][1][2].lower() == 'x,y,z']

    pdb_interfaces = [p for p in pdb_interfaces if p[7][0][1] == 'Protein' and p[7][1][1] == 'Protein']

    pdb_interfaces = filter_non_one_phos(pdb_interfaces)
    pdb_interfaces = filter_same_interface(pdb_interfaces)



if __name__ == "__main__":
    main()


