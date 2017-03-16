#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import urllib
import urllib2
import cPickle as pickle
from multiprocessing import Pool

"""
filter non_bond, low-resolution, same-chain, non-phos, two-phos interfaces
"""
def write_pdb_interfaces(pdb_interfaces,fname):
    with open(fname+'.txt','w') as w_f:
        for i in pdb_interfaces:
            pdbid,interface_id,interface_area,inter_solvent,pvalue,binding_energy,interaction_res,chains = i
            print >> w_f, '{0:<8}{1:<8}{2:<8.2f}{3:<8.2f}{4:<8.2f}{5:<8.2f}{6:<}\t{7:<}'.format(pdbid,interface_id,interface_area,inter_solvent,pvalue,binding_energy,interaction_res,chains)
        pickle.dump(pdb_interfaces,open(fname+'.pickle','w'))


def filter_resolution(pdb_interfaces,cutoff):

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
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line for line in lines if line]
    lines = [line.split(',') for line in lines]
    lines = [[w.strip('"') for w in line] for line in lines]

    # filter pdbs with lower resolution
    lines = [line for line in lines[1:] if line[1]]
    pdbids = [line[0] for line in lines if float(line[1]) <= cutoff]
    good_pdb_interfaces = [p for p in pdb_interfaces if p[0] in pdbids]

    print 'after filter low resolution',cutoff,len(good_pdb_interfaces)
    return good_pdb_interfaces


def filter_symmety(pdb_interfaces):
    xyz_interfaces = [p for p in pdb_interfaces if p[7][0][2].lower() == 'x,y,z' and p[7][1][2].lower() == 'x,y,z']
    print 'after filter symmetry',len(xyz_interfaces)
    return xyz_interfaces

def filter_non_pro_interface(pdb_interfaces):
    non_pro_pdb_interfaces = [p for p in pdb_interfaces if p[7][0][1] != 'Protein' or p[7][1][1] != 'Protein']
    pro_pdb_interfaces = [p for p in pdb_interfaces if p[7][0][1] == 'Protein' and p[7][1][1] == 'Protein']
    print 'after non protein',len(pro_pdb_interfaces)
    return pro_pdb_interfaces

def filter_cov_interface(pdb_interfaces):
    non_cov_pdb_interfaces = [p for p in pdb_interfaces if not 'cov_bonds' in [bt[0] for bt in p[6]] and not 'ss_bonds' in [bt[0] for bt in p[6]]]
    cov_pdb_interfaces = [p for p in pdb_interfaces if 'cov_bonds' in [bt[0] for bt in p[6]] or 'ss_bonds' in [bt[0] for bt in p[6]]]
    print 'after filter covalent and ss bonds',len(non_cov_pdb_interfaces)
    return non_cov_pdb_interfaces


def filter_non_phos(pdb_interfaces):
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

    write_pdb_interfaces(more_phos_interfaces,'more_phos')
    write_pdb_interfaces(zero_phos_interfaces,'zero_phos')
    write_pdb_interfaces(one_phos_interfaces,'one_phos')

    print 'zero phos interfaces',len(zero_phos_interfaces)
    print 'more phos interfaces',len(more_phos_interfaces)
    print 'after filter non_one_phos_interfaces',len(one_phos_interfaces+more_phos_interfaces)
    return one_phos_interfaces + more_phos_interfaces


def filter_small_interface(pdb_interfaces,cutoff):
    big_pdb_interfaces = [p for p in pdb_interfaces if p[2] >= cutoff]
    small_pdb_interfaces = [p for p in pdb_interfaces if p[2] < cutoff]
    write_pdb_interfaces(small_pdb_interfaces,'small')
    print 'after filter small interface',len(big_pdb_interfaces)
    return big_pdb_interfaces

def filter_unstable_assembly(pdb_interfaces):
    stable_interfaces = pickle.load(open('stable_pdb_interfaces.pickle'))
    stable_pdb_if = []
    for pdbid,interfaces in stable_interfaces.iteritems():
        for i in interfaces:
            stable_pdb_if.append(pdbid.upper()+i)
    stable_pdb_interfaces = [p for p in pdb_interfaces if p[0].upper()+p[1] in stable_pdb_if]
    unstable_pdb_interfaces = [p for p in pdb_interfaces if not p[0].upper()+p[1] in stable_pdb_if]

    write_pdb_interfaces(unstable_pdb_interfaces,'unstable')

    print 'after filter unstable assembly',len(stable_pdb_interfaces)
    return stable_pdb_interfaces

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

def filter_kinase(pdb_interfaces):

    KINASE = ['PF00069']
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

    non_kinase_result = [(r[0],r[1]) for r in result if not r[2] in KINASE]
    kinase_result = [(r[0],r[1]) for r in result if r[2] in KINASE]

    non_kinae_pdb_interfaces = [p for p in pdb_interfaces if (p[0],p[1]) in non_kinase_result]
    kinae_pdb_interfaces = [p for p in pdb_interfaces if (p[0],p[1]) in kinase_result]

    write_pdb_interfaces(kinae_pdb_interfaces,'kinase')
    print 'after filter kinae',len(non_kinae_pdb_interfaces)
    return non_kinae_pdb_interfaces




def main():
    interface_cutoff = 400
    resolution_cutoff = 3.0
    pdb_interfaces = pickle.load(open(sys.argv[-1]))

    pdb_interfaces = filter_resolution(pdb_interfaces,resolution_cutoff)
    # write_pdb_interfaces(pdb_interfaces,'0filter_resolution')

    pdb_interfaces = filter_cov_interface(pdb_interfaces)
    # write_pdb_interfaces(pdb_interfaces,'2filter_cov_interface')

    pdb_interfaces = filter_non_pro_interface(pdb_interfaces)
    # write_pdb_interfaces(pdb_interfaces,'1filter_non_pro')

    pdb_interfaces = filter_non_phos(pdb_interfaces)
    # write_pdb_interfaces(pdb_interfaces,'3filter_non_phos')

    pdb_interfaces = filter_small_interface(pdb_interfaces,interface_cutoff)
    # write_pdb_interfaces(pdb_interfaces,'5filter_small')

    pdb_interfaces = filter_unstable_assembly(pdb_interfaces)
    # write_pdb_interfaces(pdb_interfaces,'4filter_unstable')

    pdb_interfaces = filter_kinase(pdb_interfaces)
    # write_pdb_interfaces(pdb_interfaces,'6filter_kinase')

    write_pdb_interfaces(pdb_interfaces,'base_filter')

if __name__ == "__main__":
    main()

