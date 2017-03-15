#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
parse pisa interface xml file
"""

import os
import sys
from bs4 import BeautifulSoup
from collections import OrderedDict
from collections import Counter
import cPickle as pickle
from multiprocessing import Pool

def parse_interface(pisa_f):
    soup = BeautifulSoup(open(pisa_f),'lxml')

    pdb_interfaces = []
    pdb_entries = soup.find_all('pdb_entry')
    if pdb_entries:
        for pdb_entry in pdb_entries:
            pdbid = pdb_entry.find('pdb_code').get_text()
            n_interfaces = int(pdb_entry.find('n_interfaces').get_text())
            if n_interfaces > 0:
                interfaces = pdb_entry.find_all('interface')
                for interface in interfaces:
                    interface_id = interface.find('id').get_text()
                    single_interface = [pdbid,interface_id]
                    int_area = float(interface.find('int_area').get_text())
                    int_solv_en = float(interface.find('int_solv_en').get_text())
                    pvalue = float(interface.find('pvalue').get_text())
                    stab_en = float(interface.find('stab_en').get_text())

                    single_interface += [int_area,int_solv_en,pvalue,stab_en]
                    all_bonds = []
                    hydrogen_bonds = []
                    h_bonds = interface.find('h-bonds')
                    bond_num = int(h_bonds.find('n_bonds').get_text())
                    if bond_num > 0:
                        bonds = h_bonds.find_all('bond')
                        for bond in bonds:
                            chain1 = bond.find('chain-1').get_text()
                            res1 = bond.find('res-1').get_text()
                            seqnum1 = bond.find('seqnum-1').get_text()
                            atname1 = bond.find('atname-1').get_text()

                            chain2 = bond.find('chain-2').get_text()
                            res2 = bond.find('res-2').get_text()
                            seqnum2 = bond.find('seqnum-2').get_text()
                            atname2 = bond.find('atname-2').get_text()

                            dist = bond.find('dist').get_text()

                            bond_info = [chain1+'_'+seqnum1+'_'+res1+'_'+atname1,chain2+'_'+seqnum2+'_'+res2+'_'+atname2,dist]

                            hydrogen_bonds.append(bond_info)
                    if hydrogen_bonds:
                        all_bonds.append(('hydrogen_bonds',hydrogen_bonds))

                    salt_bridges = []
                    s_bridges = interface.find('salt-bridges')
                    bond_num = int(s_bridges.find('n_bonds').get_text())
                    if bond_num > 0:
                        bonds = s_bridges.find_all('bond')
                        for bond in bonds:
                            chain1 = bond.find('chain-1').get_text()
                            res1 = bond.find('res-1').get_text()
                            seqnum1 = bond.find('seqnum-1').get_text()
                            atname1 = bond.find('atname-1').get_text()

                            chain2 = bond.find('chain-2').get_text()
                            res2 = bond.find('res-2').get_text()
                            seqnum2 = bond.find('seqnum-2').get_text()
                            atname2 = bond.find('atname-2').get_text()

                            dist = bond.find('dist').get_text()

                            bond_info = [chain1+'_'+seqnum1+'_'+res1+'_'+atname1,chain2+'_'+seqnum2+'_'+res2+'_'+atname2,dist]

                            salt_bridges.append(bond_info)
                    if salt_bridges:
                        all_bonds.append(('salt_bridges',salt_bridges))

                    ss_bonds = []
                    s_bonds = interface.find('ss-bonds')
                    bond_num = int(s_bonds.find('n_bonds').get_text())
                    if bond_num > 0:
                        bonds = s_bonds.find_all('bond')
                        for bond in bonds:
                            chain1 = bond.find('chain-1').get_text()
                            res1 = bond.find('res-1').get_text()
                            seqnum1 = bond.find('seqnum-1').get_text()
                            atname1 = bond.find('atname-1').get_text()

                            chain2 = bond.find('chain-2').get_text()
                            res2 = bond.find('res-2').get_text()
                            seqnum2 = bond.find('seqnum-2').get_text()
                            atname2 = bond.find('atname-2').get_text()

                            dist = bond.find('dist').get_text()

                            bond_info = [chain1+'_'+seqnum1+'_'+res1+'_'+atname1,chain2+'_'+seqnum2+'_'+res2+'_'+atname2,dist]

                            ss_bonds.append(bond_info)

                    if ss_bonds:
                        all_bonds.append(('ss_bonds',ss_bonds))

                    cov_bonds = []
                    c_bonds = interface.find('cov-bonds')
                    bond_num = int(c_bonds.find('n_bonds').get_text())
                    if bond_num > 0:
                        bonds = c_bonds.find_all('bond')
                        for bond in bonds:
                            chain1 = bond.find('chain-1').get_text()
                            res1 = bond.find('res-1').get_text()
                            seqnum1 = bond.find('seqnum-1').get_text()
                            atname1 = bond.find('atname-1').get_text()

                            chain2 = bond.find('chain-2').get_text()
                            res2 = bond.find('res-2').get_text()
                            seqnum2 = bond.find('seqnum-2').get_text()
                            atname2 = bond.find('atname-2').get_text()

                            dist = bond.find('dist').get_text()

                            bond_info = [chain1+'_'+seqnum1+'_'+res1+'_'+atname1,chain2+'_'+seqnum2+'_'+res2+'_'+atname2,dist]

                            cov_bonds.append(bond_info)
                    if cov_bonds:
                        all_bonds.append(('cov_bonds',cov_bonds))

                    single_interface.append(all_bonds)

                    molecules = interface.find_all('molecule')
                    mols = []
                    for mol in molecules:
                        mol_chain = mol.find('chain_id').get_text()
                        mol_class = mol.find('class').get_text()
                        mol_symop = mol.find('symop').get_text()
                        mol_int_nres = int(mol.find('int_nres').get_text())
                        mol_res = [mol_chain+'_'+res.find('seq_num').get_text()+'_'+res.find('name').get_text() for res in mol.find_all('residue')]
                        mols.append((mol_chain,mol_class,mol_symop,mol_int_nres,mol_res))
                    single_interface.append(mols)

                    if single_interface:
                        pdb_interfaces.append(single_interface)



            else:
                print pdbid+' is empty'

    if pdb_interfaces:
        print pisa_f,'parse complete'
        return pdb_interfaces
    else:
        print  pisa_f,'got nothing'



def main():
    parameters = []
    for root,dirs,files in os.walk(sys.argv[-1]):
        for f in files:
            if f[-4:] == '.xml':
                f = os.path.join(root,f)
                parameters.append(f)

    # result = []
    # for p in parameters:
        # r =parse_interface(p)
        # result.append(r)
        # break

    p = Pool(6)
    result = p.map(parse_interface,parameters)
    p.close()

    pdb_interfaces =[]
    for r in result:
        if r:
            pdb_interfaces += r

    with open('pdb_interfaces.txt','w') as w_f:
        for i in pdb_interfaces:
            pdbid,interface_id,interface_area,inter_solvent,pvalue,binding_energy,interaction_res,chains = i
            print >> w_f, '{0:<8}{1:<8}{2:<8.2f}{3:<8.2f}{4:<8.2f}{5:<8.2f}\t{6:<}\t{7:<}'.format(pdbid,interface_id,interface_area,inter_solvent,pvalue,binding_energy,interaction_res,chains)
    pickle.dump(pdb_interfaces,open('pdb_interfaces.pickle','w'))


if __name__ == "__main__":
    main()




