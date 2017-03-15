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



def main():
    pdb_interfaces = pickle.load(open(sys.argv[-1]))

    pdb_interfaces = filter_symmety(pdb_interfaces)
    write_pdb_interfaces(pdb_interfaces,'filter_symmety')

    pdb_interfaces = filter_non_pro_interface(pdb_interfaces)
    write_pdb_interfaces(pdb_interfaces,'filter_non_pro')

    pdb_interfaces = filter_cov_interface(pdb_interfaces)
    write_pdb_interfaces(pdb_interfaces,'filter_cov_interface')

    pdb_interfaces = filter_non_phos(pdb_interfaces)
    write_pdb_interfaces(pdb_interfaces,'filter_non_phos')



if __name__ == "__main__":
    main()

