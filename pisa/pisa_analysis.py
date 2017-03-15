#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
analysis interface from pisa
"""
import os
import sys
import urllib
import urllib2
import numpy as np
from multiprocessing import Pool
import cPickle as pickle


def filter_resolution(pdb_interfaces,resolution_cutoff):
    pdbids = [p[0] for p in pdb_interfaces]
    pdbids = ','.join(pdbids)
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbids,
        'customReportColumns':'structureId,resolution',
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
    lines = [line for line in lines if line[-1]]
    lines = [line[0] for line in lines if float(line[-1]) <= resolution_cutoff]
    pdb_interfaces = [p for p in pdb_interfaces if p[0] in lines]
    print 'after filter resolution',len(pdb_interfaces),'is below',resolution_cutoff
    return pdb_interfaces

def filter_pisa_interface(pdb_interfaces):
    plip_pdb_interfaces = pickle.load(open('rcsb_phos_pdb_lig_interactions_interface.pickle','r'))
    plip_pdb_interfaces = [(pdbid,face) for pdbid,_,face,_ in plip_pdb_interfaces]
    real_pdb_interfaces = []
    for interface in pdb_interfaces:
        pdbid,_,interface_area,_,_,_,bonds = interface[:7]
        chains = interface[-1]
        face = chains[0][0] + '-' + chains[1][0]
        if (pdbid,face) in plip_pdb_interfaces:
            real_pdb_interfaces.append(interface)
    print 'after filter pisa',len(real_pdb_interfaces)
    return real_pdb_interfaces

def select_important(pdb_interfaces):
    detail_interfaces = pickle.load(open('detail_analyse.pickle'))
    important_interfaces = []
    for interface in pdb_interfaces:
        pdbid,interface_id,interface_area,_,_,_,bonds = interface[:7]
        pdbid = pdbid.lower()
        chains = interface[-1]
        key = pdbid+'_'+interface_id
        if key in detail_interfaces.keys():
            inter_chain,chain_info,inter_info = detail_interfaces[pdbid+'_'+interface_id]
            important_res = inter_info[0][1] + inter_info[1][1]
            important_res = ''.join(important_res)
            if 'TPO' in  important_res or 'SEP' in important_res or 'PTR' in important_res:
                important_interfaces.append(interface)
    print 'phos_res as stabilizing res',len(important_interfaces)
    return important_interfaces


def select_hetero_interface(pdb_interfaces):
    hetero_pdb_interfaces = []
    for interface in pdb_interfaces:
        chains = interface[-1]
        if len(chains) == 2:
            chain1 = chains[0][0]
            chain2 = chains[1][0]
            if chain1 != chain2 and len(chain1) == 1 and len(chain2) == 1:
                hetero_pdb_interfaces.append(interface)
    if hetero_pdb_interfaces:
        print 'hetero_phos_interfaces',len(hetero_pdb_interfaces)
        return hetero_pdb_interfaces


def select_phos(pdb_interfaces):
    # filter interface with no bonds infomation
    pdb_interfaces = [interface for interface in pdb_interfaces if interface[7]]
    phos_interfaces = []
    for interface in pdb_interfaces:
        pdbid,_,interface_area,_,_,_,bonds = interface[:7]
        chains = interface[-1]
        for bond in bonds:
            bond_type,bond_info = bond
            for bondi in bond_info:
                res1,res2,dist = bondi
                b = res1+'-'+res2
                if 'TPO' in b or 'SEP' in b or 'PTR' in b:
                    if not interface in phos_interfaces:
                        phos_interfaces.append(interface)
                        break
    if phos_interfaces:
        print 'phos_interfaces',len(phos_interfaces)
        return phos_interfaces

def select_phos_site(pdb_interfaces):
    phos_interfaces_site = []
    for interface in pdb_interfaces:
        pdbid,p1,interface_area,p2,p3,p4,bonds = interface[:7]
        chains = interface[-1]
        phos_bonds = []
        for bond in bonds:
            bond_type,bond_info = bond
            if bond_type == 'hydrogen_bonds':
                hydrogen_bonds = []
                for bondi in bond_info:
                    res1,res2,dist = bondi
                    b = res1+'-'+res2
                    if 'TPO_ O1P' in b or 'TPO_ O2P' in b or 'TPO_ O3P' in b or 'TPO_ OG1' in b:
                        hydrogen_bonds.append(bondi)
                    if 'SEP_ O1P' in b or 'SEP_ O2P' in b or 'SEP_ O3P' in b or 'SEP_ OG' in b:
                        hydrogen_bonds.append(bondi)
                    if 'PTR_ O1P' in b or 'PTR_ O2P' in b or 'PTR_ O3P' in b or 'PTR_ OH' in b:
                        hydrogen_bonds.append(bondi)
                if hydrogen_bonds:
                    phos_bonds.append((bond_type,hydrogen_bonds))
            elif bond_type == 'salt_bridges':
                salt_bridges = []
                for bondi in bond_info:
                    res1,res2,dist = bondi
                    b = res1+'-'+res2
                    if 'TPO_ O1P' in b or 'TPO_ O2P' in b or 'TPO_ O3P' in b or 'TPO_ OG1' in b:
                        salt_bridges.append(bondi)
                    if 'SEP_ O1P' in b or 'SEP_ O2P' in b or 'SEP_ O3P' in b or 'SEP_ OG' in b:
                        salt_bridges.append(bondi)
                    if 'PTR_ O1P' in b or 'PTR_ O2P' in b or 'PTR_ O3P' in b or 'PTR_ OH' in b:
                        salt_bridges.append(bondi)
                if salt_bridges:
                    phos_bonds.append((bond_type,salt_bridges))

        if phos_bonds:
            phos_interfaces_site.append([pdbid,p1,interface_area,p2,p3,p4,phos_bonds,chains])

    if phos_interfaces_site:
        print 'phos_interfaces_site',len(phos_interfaces_site)
        return phos_interfaces_site

def select_one_phos_site(pdb_interfaces):
    one_phos_interfaces = []
    for interface in pdb_interfaces:
        pdbid,p1,interface_area,p2,p3,p4,bonds = interface[:7]
        chains = interface[-1]
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
        else:
            print pdbid,p1,phos_res

    if one_phos_interfaces:
        print 'one_phos_interfaces',len(one_phos_interfaces)
        return one_phos_interfaces



def main():
    pdb_interfaces = pickle.load(open(sys.argv[-1]))

    resolution_cutoff = 3.0
    pdb_interfaces = filter_resolution(pdb_interfaces,resolution_cutoff)
    pdb_interfaces = filter_pisa_interface(pdb_interfaces)
    pdb_interfaces = select_important(pdb_interfaces)


    if pdb_interfaces:
        pdb_interfaces = select_phos(pdb_interfaces)
        if pdb_interfaces:
            with open('pdb_phos_interfaces.txt','w') as w_f:
                for i in pdb_interfaces:
                    print >> w_f, i
            pickle.dump(pdb_interfaces,open('pdb_phos_interfaces.pickle','w'))

    if pdb_interfaces:
        pdb_interfaces = select_phos_site(pdb_interfaces)
        if pdb_interfaces:
            with open('pdb_phos_interfaces_site.txt','w') as w_f:
                for i in pdb_interfaces:
                    print >> w_f, i
            pickle.dump(pdb_interfaces,open('pdb_phos_interfaces_site.pickle','w'))

    if pdb_interfaces:
        pdb_interfaces = select_hetero_interface(pdb_interfaces)
        if pdb_interfaces:
            with open('hetero_pdb_phos_interfaces_site.txt','w') as w_f:
                for i in pdb_interfaces:
                    print >> w_f, i
            pickle.dump(pdb_interfaces,open('hetero_pdb_phos_interfaces_site.pickle','w'))

    if pdb_interfaces:
        pdb_interfaces = select_one_phos_site(pdb_interfaces)
        if pdb_interfaces:
            with open('one_phos_interfaces_site.txt','w') as w_f:
                for i in pdb_interfaces:
                    print >> w_f, i
            pickle.dump(pdb_interfaces,open('one_phos_interfaces_site.pickle','w'))


if __name__ == "__main__":
    main()







