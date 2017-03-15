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

def filter_symmety(pdb_interfaces):
    xyz_interfaces = [p for p in pdb_interfaces if p[7][0][2].lower() == 'x,y,z' and p[7][1][2].lower() == 'x,y,z']
    print 'after filter symmetry',len(xyz_interfaces)
    return xyz_interfaces

def filter_non_pro_interface(pdb_interfaces):
    non_pro_pdb_interfaces = [p for p in pdb_interfaces if p[7][0][1] != 'Protein' or p[7][1][1] != 'Protein']
    pro_pdb_interfaces = [p for p in pdb_interfaces if p[7][0][1] == 'Protein' and p[7][1][1] == 'Protein']
    write_pdb_interfaces(non_pro_pdb_interfaces,'non_pro')
    print 'after non protein',len(pro_pdb_interfaces)
    return pro_pdb_interfaces

def filter_cov_interface(pdb_interfaces):
    non_cov_pdb_interfaces = [p for p in pdb_interfaces if not 'cov_bonds' in [bt[0] for bt in p[6]] and not 'ss_bonds' in [bt[0] for bt in p[6]]]
    cov_pdb_interfaces = [p for p in pdb_interfaces if 'cov_bonds' in [bt[0] for bt in p[6]] or 'ss_bonds' in [bt[0] for bt in p[6]]]
    write_pdb_interfaces(cov_pdb_interfaces,'cov')
    print 'after filter covalent and ss bonds',len(non_cov_pdb_interfaces)
    return non_cov_pdb_interfaces

def filter_short_peptide(pdb_interfaces,num_cutoff):
    long_pdb_interfaces = [p for p in pdb_interfaces if p[7][0][3] >= 3 and p[7][1][3] >= 3]
    short_pdb_interfaces = [p for p in pdb_interfaces if p[7][0][3] < 3 or p[7][1][3] < 3]
    short_pdb_inter_num = [p[7] for p in short_pdb_interfaces]
    with open('short_num.txt','w') as w_f:
        for s in short_pdb_inter_num:
            print >> w_f,s
    write_pdb_interfaces(short_pdb_interfaces,'short_'+str(num_cutoff))
    print 'after filter short peptide',len(long_pdb_interfaces)
    return long_pdb_interfaces

def filter_same_chain(pdb_interfaces):
    hetero_pdb_interfaces = [p for p in pdb_interfaces if p[7][0][0] != p[7][1][0]]
    same_pdb_interfaces = [p for p in pdb_interfaces if p[7][0][0] == p[7][1][0]]
    write_pdb_interfaces(same_pdb_interfaces,'same_chain')
    print 'after filter same chain',len(hetero_pdb_interfaces)
    return hetero_pdb_interfaces


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

    write_pdb_interfaces(more_phos_interfaces,'more_phos')
    write_pdb_interfaces(zero_phos_interfaces,'zero_phos')

    print 'zero phos interfaces',len(zero_phos_interfaces)
    print 'more phos interfaces',len(more_phos_interfaces)
    print 'after filter non_one_phos_interfaces',len(one_phos_interfaces)
    return one_phos_interfaces


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

    p = Pool(8)
    result = p.map(get_entityid,pdbid_chain)
    p.close()

    pdb_chain_entity = []
    pdb_unique_interface = []
    for r in result:
        if (r[0],r[2],r[3]) in pdb_chain_entity:
            pass
        else:
            pdb_unique_interface.append((r[0],r[1]))
            pdb_chain_entity.append((r[0],r[2],r[3]))

    pdb_interfaces = [p for p in pdb_interfaces if (p[0],p[1]) in pdb_unique_interface]
    print 'after filter same entity',len(pdb_interfaces)
    return pdb_interfaces

# def filter_non_standard(pdb_interfaces):
    # AA = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', 'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}
    # stanard_interfaces = []
    # non_standard_interfaces = []
    # for interface in pdb_interfaces:
        # pdbid,p1,interface_area,p2,p3,p4,bonds = interface[:7]
        # for bond in bonds:
            # bond_type,bond_info = bond
            # for bondi in bond_info:
                # res1,res2,dist = bondi
                # if 'TPO' in res1 or 'SEP' in res1 or 'PTR' in res1:
                    # if res2.split('_')[2] not in AA.keys():
                        # non_standard_interfaces.append(interface)
                        # break
                # if 'TPO' in res2 or 'SEP' in res2 or 'PTR' in res2:
                    # if res1.split('_')[2] not in AA.keys():
                        # non_standard_interfaces.append(interface)
                        # break
        # stanard_interfaces.append(interface)


    # print 'after non_standard interaction residues',len(stanard_interfaces)
    # return stanard_interfaces




def main():
    pdb_interfaces = pickle.load(open(sys.argv[-1]))

    pdb_interfaces = filter_symmety(pdb_interfaces)
    write_pdb_interfaces(pdb_interfaces,'filter_symmety')

    pdb_interfaces = filter_non_pro_interface(pdb_interfaces)
    write_pdb_interfaces(pdb_interfaces,'filter_non_pro')

    pdb_interfaces = filter_cov_interface(pdb_interfaces)
    write_pdb_interfaces(pdb_interfaces,'filter_cov_interface')

    pdb_interfaces = filter_short_peptide(pdb_interfaces,6)
    write_pdb_interfaces(pdb_interfaces,'filter_short_peptide')

    pdb_interfaces = filter_same_chain(pdb_interfaces)
    write_pdb_interfaces(pdb_interfaces,'filter_same_chain')

    pdb_interfaces = filter_non_one_phos(pdb_interfaces)
    write_pdb_interfaces(pdb_interfaces,'filter_non_one_phos')

    # resolution_cutoff = 3.0
    # pdb_interfaces = filter_resolution(pdb_interfaces,resolution_cutoff)
    # write_pdb_interfaces(pdb_interfaces,'filter_low_resolution')


    # pdb_interfaces = filter_non_standard(pdb_interfaces)
    # write_pdb_interfaces(pdb_interfaces,'filter_non_standard')


    # pdb_interfaces = filter_same_interface(pdb_interfaces)
    # write_pdb_interfaces(pdb_interfaces,'filter_same_entity')

    write_pdb_interfaces(pdb_interfaces,'final')


if __name__ == "__main__":
    main()

